"""Compute quantum chemistry using Mainz-Austin-Budapest-Gainesville's CFOUR executable."""

import copy
import pprint
from decimal import Decimal
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
from qcelemental.models import AtomicInput, AtomicResult, Provenance
from qcelemental.util import safe_version, unnp, which

from ...util import execute
from ..model import ProgramHarness
from ..qcvar_identities_resources import build_atomicproperties, build_out
from .germinate import muster_modelchem
from .harvester import harvest
from .keywords import format_keywords

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class CFOURHarness(ProgramHarness):
    """

    Notes
    -----
    * Looks for basis set file ``../basis/GENBAS`` from ``xcfour`` executable. If this doesn't work, file an issue.

    """

    _defaults = {
        "name": "CFOUR",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "xcfour", return_bool=True, raise_error=raise_error, raise_msg="Please install via http://cfour.de/"
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("xcfour")
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "ZMAT"], {"ZMAT": "\nHe\n\n"})

            if success:
                for line in output["stdout"].splitlines():
                    if "Version" in line:
                        branch = " ".join(line.strip().split()[1:])
                self.version_cache[which_prog] = safe_version(branch)

        return self.version_cache[which_prog]

    def compute(self, input_model: AtomicInput, config: "TaskConfig") -> AtomicResult:
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_model)

    def build_input(
        self, input_model: AtomicInput, config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        cfourrec = {"infiles": {}, "scratch_directory": config.scratch_directory}

        opts = copy.deepcopy(input_model.keywords)

        # Handle memory
        # for cfour, [GiB] --> [MB]
        opts["memory_size"] = int(config.memory * (1024 ** 3) / 1e6)
        opts["mem_unit"] = "mb"

        # Handle molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="cfour", units="Bohr", return_data=True)
        opts.update(moldata["keywords"])

        # Handle calc type and quantum chemical method
        mdcopts = muster_modelchem(input_model.model.method, input_model.driver.derivative_int())
        opts.update(mdcopts)

        # Handle basis set
        # * why, yes, this is highly questionable
        #   * assuming relative file location between xcfour exe and GENBAS file
        #   * reading a multi MB file into the inputs dict
        opts["basis"] = input_model.model.basis

        # Handle conversion from schema (flat key/value) keywords into local format
        optcmd = format_keywords(opts)

        xcfour = which("xcfour")
        genbas = Path(xcfour).parent.parent / "basis" / "GENBAS"
        cfourrec["infiles"]["ZMAT"] = molcmd + optcmd
        cfourrec["infiles"]["GENBAS"] = genbas.read_text()
        cfourrec["command"] = [xcfour]

        return cfourrec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            ["GRD", "FCMFINAL", "DIPOL"],
            scratch_messy=False,
            scratch_directory=inputs["scratch_directory"],
        )
        return success, dexe

    def parse_output(
        self, outfiles: Dict[str, str], input_model: AtomicInput
    ) -> AtomicResult:  # lgtm: [py/similar-function]

        stdout = outfiles.pop("stdout")
        stderr = outfiles.pop("stderr")

        # c4mol, if it exists, is dinky, just a clue to geometry of cfour results
        qcvars, c4hess, c4grad, c4mol, version, errorTMP = harvest(input_model.molecule, stdout, **outfiles)

        if c4grad is not None:
            qcvars["CURRENT GRADIENT"] = c4grad
            qcvars[f"{input_model.model.method.upper()[3:]} TOTAL GRADIENT"] = c4grad

        if c4hess is not None:
            qcvars[f"{input_model.model.method.upper()[3:]} TOTAL HESSIAN"] = c4hess
            qcvars["CURRENT HESSIAN"] = c4hess

        if input_model.driver.upper() == "PROPERTIES":
            retres = qcvars[f"CURRENT ENERGY"]
        else:
            retres = qcvars[f"CURRENT {input_model.driver.upper()}"]

        if isinstance(retres, Decimal):
            retres = float(retres)
        elif isinstance(retres, np.ndarray):
            retres = retres.ravel().tolist()

        build_out(qcvars)
        atprop = build_atomicproperties(qcvars)

        output_data = {
            "schema_version": 1,
            "extras": {"outfiles": outfiles, **input_model.extras},
            "properties": atprop,
            "provenance": Provenance(creator="CFOUR", version=self.get_version(), routine="xcfour"),
            "return_result": retres,
            "stderr": stderr,
            "stdout": stdout,
            "success": True,
        }

        # got to even out who needs plump/flat/Decimal/float/ndarray/list
        # Decimal --> str preserves precision
        output_data["extras"]["qcvars"] = {
            k.upper(): str(v) if isinstance(v, Decimal) else v for k, v in unnp(qcvars, flat=True).items()
        }

        return AtomicResult(**{**input_model.dict(), **output_data})
