"""
Calls the Madness moldft executable.
"""
import re
import copy
import logging
import pprint
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple

import numpy as np
import qcelemental as qcel
from qcelemental.models import AtomicResult, Provenance, AtomicInput
from qcelemental.util import safe_version, which, which_import

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


class MadnessHarness(ProgramHarness):
    """
     Notes
     -----
     * To use the TCE, specify ``AtomicInput.model.method`` as usual, then also include ``qc_module = True`` in ``AtomicInput.keywords``.
     """

    _defaults = {
        "name": "Madness",
        "scratch": True,
        "thread_safe": True,
        "thread_parallel": True,
        "node_parallel": True,
        "managed_memory": True,
    }
    # ATL: OpenMP only >=6.6 and only for Phi; potential for Mac using MKL and Intel compilers
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Whether Madness harness is ready for operation, with both the QC program and any particular dependencies found.

        Parameters
        ----------
        raise_error: bool
            Passed on to control negative return between False and ModuleNotFoundError raised.

         Returns
         -------
         bool
             If both nwchem and its harness dependency networkx are found, returns True.
             If raise_error is False and nwchem or networkx are missing, returns False.
             If raise_error is True and nwchem or networkx are missing, the error message for the first missing one is raised.

       """
        qc = which(
            "madness",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via https://github.com/m-a-d-n-e-s-s/madness",
        )
        #         dep = which_import(
        #             "networkx",
        #             return_bool=True,
        #             raise_error=raise_error,
        #             raise_msg="For NWChem harness, please install via `conda install networkx -c conda-forge`.",
        #         )
        return qc  # and dep

    ## gotta figure out which input file and from where
    def get_version(self) -> str:
        self.found(raise_error=True)

        # Get the node configuration
        config = get_config()

        # Run MADNESS
        which_prog = which("madness")
        if config.use_mpiexec:
            command = create_mpi_invocation(which_prog, config)
        else:
            command = [which_prog]
        command.append("v.moldft")

        if which_prog not in self.version_cache:
            success, output = execute(
                command,
                {
                    "v.moldft": "dft\nxc hf\nend\n\ngeometry\n He    0.00000000      0.00000000     0.00000000    \n end\n"
                },
                scratch_directory=config.scratch_directory,
            )

            if success:
                for line in output["stdout"].splitlines():
                    if "multiresolution suite" in line:
                        version = line.strip().split()[1]
                self.version_cache[which_prog] = safe_version(version)
            else:
                raise UnknownError(output["stderr"])

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
         Runs madness in executable mode
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
        self, input_model: AtomicInput, config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        #
        madnessrec = {"infiles": {}, "scratch_directory": config.scratch_directory}

        opts = copy.deepcopy(input_model.keywords)
        opts = {k.lower(): v for k, v in opts.items()}

        # Handle Molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="madness", units="bohr", return_data=True)
        molData = {}
        for k, v in moldata["keywords"].items():
            molData["dft__" + k] = v
        opts.update(molData)

        ## Handle Calc Type (ROBERT)
        mdccmd, mdcopts = muster_modelchem(input_model.model.method, input_model.driver, opts.pop("qc_module", False))
        opts.update(mdcopts)

        ## Handle the basis set (ROBERT) the question is what value of k

        # Log the job settings (LORI)  Not sure if i need this
        logger.debug("JOB_OPTS")
        logger.debug(pp.pformat(opts))

        # Handle conversion from schema (flat key/value) keywords into local format
        optcmd = format_keywords(opts)

        # optcmd="dft\n xc hf \nend\n"

        madnessrec["infiles"]["input"] = optcmd + molcmd
        ## Determine the command
        # Determine the command
        madnessrec["command"] = [which("madness")]
        print(madnessrec["infiles"]["input"])
        return madnessrec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:
        success, dexe = execute(inputs["command"], inputs["infiles"],)
        return success, dexe

    def parse_output(
        self, outfiles: Dict[str, str], input_model: "AtomicInput"
    ) -> "AtomicResult":  # lgtm: [py/similar-function]

        # Get the stdout from the calculation (required)
        stdout = outfiles.pop("stdout")

        # Read the NWChem stdout file and, if needed, the hess or grad files
        qcvars, nwhess, nwgrad, nwmol, version, errorTMP = harvest(input_model.molecule, stdout, **outfiles)

        if nwgrad is not None:
            qcvars["CURRENT GRADIENT"] = nwgrad

        if nwhess is not None:
            qcvars["CURRENT HESSIAN"] = nwhess

        # Normalize the output as a float or list of floats
        if input_model.driver.upper() == "PROPERTIES":
            retres = qcvars[f"CURRENT ENERGY"]
        else:
            print(qcvars)
            retres = qcvars[f"CURRENT {input_model.driver.upper()}"]

        if isinstance(retres, Decimal):
            retres = float(retres)
        elif isinstance(retres, np.ndarray):
            retres = retres.tolist()

        # Get the formatted properties
        qcprops = extract_formatted_properties(qcvars)

        # Format them inout an output
        output_data = {
            "schema_name": "qcschema_output",
            "schema_version": 1,
            "extras": {"outfiles": outfiles, **input_model.extras},
            "properties": qcprops,
            "provenance": Provenance(creator="MADNESS", version=self.get_version(), routine="madness"),
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
