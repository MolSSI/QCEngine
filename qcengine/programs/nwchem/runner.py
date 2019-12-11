"""
Calls the NWChem executable.
"""
import copy
import logging
import pprint
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple

import numpy as np
import qcelemental as qcel
from qcelemental.models import AtomicResult, Provenance
from qcelemental.util import safe_version, which

from qcengine.config import TaskConfig, get_config
from qcengine.exceptions import UnknownError

from ...exceptions import InputError
from ...util import execute, create_mpi_invocation
from ..model import ProgramHarness
from .germinate import muster_modelchem
from .harvester import extract_formatted_properties, harvest
from .keywords import format_keywords

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
logger = logging.getLogger(__name__)


class NWChemHarness(ProgramHarness):
    """

    Notes
    -----
    * To use the TCE, specify ``AtomicInput.model.method`` as usual, then also include ``qc_module = True`` in ``AtomicInput.keywords``.

    """

    _defaults = {
        "name": "NWChem",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": False,
        "node_parallel": True,
        "managed_memory": True,
    }
    # ATL: OpenMP only >=6.6 and only for Phi; potential for Mac using MKL and Intel compilers
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "nwchem",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via http://www.nwchem-sw.org/index.php/Download",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        # Get the node configuration
        config = get_config()

        # Run NWChem
        which_prog = which("nwchem")
        if config.use_mpiexec:
            command = create_mpi_invocation(which_prog, config)
        else:
            command = [which_prog]
        command.append("v.nw")

        if which_prog not in self.version_cache:
            success, output = execute(command, {"v.nw": ""}, scratch_directory=config.scratch_directory)

            if success:
                for line in output["stdout"].splitlines():
                    if "nwchem branch" in line:
                        branch = line.strip().split()[-1]
                    if "nwchem revision" in line:
                        revision = line.strip().split()[-1]
                self.version_cache[which_prog] = safe_version(branch + "+" + revision)
            else:
                raise UnknownError(output["stderr"])

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs NWChem in executable mode
        """
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        if "There is an error in the input file" in dexe["stdout"]:
            raise InputError(dexe["stdout"])
        if "not compiled" in dexe["stdout"]:
            # recoverable with a different compilation with optional modules
            raise InputError(dexe["stdout"])

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_model)
        else:
            raise UnknownError(dexe["stderr"])

    def build_input(
        self, input_model: "AtomicInput", config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        nwchemrec = {"infiles": {}, "scratch_directory": config.scratch_directory}

        opts = copy.deepcopy(input_model.keywords)
        opts = {k.lower(): v for k, v in opts.items()}

        # Handle memory
        # for nwchem, [GiB] --> [B]
        # someday, replace with this: opts['memory'] = str(int(config.memory * (1024**3) / 1e6)) + ' mb'
        memory_size = int(config.memory * (1024 ** 3))
        if config.use_mpiexec:  # It is the memory per MPI rank
            memory_size //= config.nnodes * config.ncores
        opts["memory"] = memory_size

        # Handle molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="nwchem", units="Bohr", return_data=True)
        opts.update(moldata["keywords"])

        # Handle calc type and quantum chemical method
        mdccmd, mdcopts = muster_modelchem(
            input_model.model.method, input_model.driver.derivative_int(), opts.pop("qc_module", False)
        )
        opts.update(mdcopts)

        # Handle basis set
        # * for nwchem, still needs sph and ghost
        for el in set(input_model.molecule.symbols):
            opts[f"basis__{el}"] = f"library {input_model.model.basis}"

        # Log the job settings
        logger.debug("JOB_OPTS")
        logger.debug(pp.pformat(opts))

        # Handle conversion from schema (flat key/value) keywords into local format
        optcmd = format_keywords(opts)

        nwchemrec["infiles"]["nwchem.nw"] = "echo\n" + molcmd + optcmd + mdccmd

        # Determine the command
        if config.nnodes > 1:
            nwchemrec["command"] = create_mpi_invocation(which("nwchem"), config)
        else:
            nwchemrec["command"] = [which("nwchem")]

        return nwchemrec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            ["job.movecs", "job.hess", "job.db", "job.zmat", "job.fd_dipole"],
            scratch_messy=False,
            scratch_directory=inputs["scratch_directory"],
        )
        return success, dexe

    def parse_output(
        self, outfiles: Dict[str, str], input_model: "AtomicInput"
    ) -> "AtomicResult":  # lgtm: [py/similar-function]

        stdout = outfiles.pop("stdout")

        # nwmol, if it exists, is dinky, just a clue to geometry of nwchem results
        qcvars, nwhess, nwgrad, nwmol, version, errorTMP = harvest(input_model.molecule, stdout, **outfiles)

        if nwgrad is not None:
            qcvars["CURRENT GRADIENT"] = nwgrad

        if nwhess is not None:
            qcvars["CURRENT HESSIAN"] = nwhess

        retres = qcvars[f"CURRENT {input_model.driver.upper()}"]
        if isinstance(retres, Decimal):
            retres = float(retres)
        elif isinstance(retres, np.ndarray):
            retres = retres.ravel().tolist()

        # Get the formatted properties
        qcprops = extract_formatted_properties(qcvars)

        # Format them inout an output
        output_data = {
            "schema_name": "qcschema_output",
            "schema_version": 1,
            "extras": {"outfiles": outfiles},
            "properties": qcprops,
            "provenance": Provenance(creator="NWChem", version=self.get_version(), routine="nwchem"),
            "return_result": retres,
            "stdout": stdout,
        }

        # got to even out who needs plump/flat/Decimal/float/ndarray/list
        # Decimal --> str preserves precision
        output_data["extras"]["qcvars"] = {
            k.upper(): str(v) if isinstance(v, Decimal) else v for k, v in qcel.util.unnp(qcvars, flat=True).items()
        }

        output_data["success"] = True
        return AtomicResult(**{**input_model.dict(), **output_data})
