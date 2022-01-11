"""
Calls the Madness moldft executable.
"""
# import re
import copy
import logging
import pprint
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple
from pathlib import Path

import numpy as np
import qcelemental as qcel
from qcelemental.models import AtomicResult, Provenance, AtomicInput
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


class MadnessHarness(ProgramHarness):
    """
    Notes
    -----
    * To use the TCE, specify ``AtomicInput.model.method`` as usual, then also include ``qc_module = True`` in ``AtomicInput.keywords``.
    """

    _defaults = {
        "name": "madness",
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
                    "v.moldft": "dft\nxc lda\nend\ngeometry\nO  0.0    0.0 0.0\nH  1.4375 0.0 1.15\nH - 1.4375 0.0 1.15\nend\n"
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
        if "There is an error in the input file" in dexe["moldft"]["stdout"]:
            raise InputError(dexe["moldft"]["stdout"])
        if "not compiled" in dexe["moldft"]["stdout"]:
            # recoverable with a different compilation with optional modules
            raise InputError(dexe["moldft"]["stdout"])
        if success:
            num_commands = len(dexe)
            print(num_commands)
            if num_commands == 2:
                dexe["moldft"]["outfiles"]["stdout"] = dexe["moldft"]["stdout"]
                dexe["moldft"]["outfiles"]["stderr"] = dexe["moldft"]["stderr"]
                dexe["molresponse"]["outfiles"]["stdout"] = dexe["molresponse"]["stdout"]
                dexe["molresponse"]["outfiles"]["stderr"] = dexe["molresponse"]["stderr"]
            else:
                dexe["moldft"]["outfiles"]["stdout"] = dexe["moldft"]["stdout"]
                dexe["moldft"]["outfiles"]["stderr"] = dexe["moldft"]["stderr"]
            return self.parse_output(dexe, input_model)
        else:
            print(dexe["stdout"])

            raise UnknownError(dexe["stderr"])

    def build_input(
            self, input_model: AtomicInput, config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        #
        madnessrec = {
            "infiles": {},
            "scratch_directory": config.scratch_directory,
            "scratch_messy": config.scratch_messy,
        }

        ## These are the madness keywords
        opts = copy.deepcopy(input_model.keywords)
        opts = {k.lower(): v for k, v in opts.items()}

        # Handle Molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="madness", units="bohr", return_data=True)
        molData = {}
        for k, v in moldata["keywords"].items():
            molData["dft__" + k] = v
        opts.update(molData)

        ## Handle Calc Type (ROBERT)
        ## now returns respnse options as well
        mdccmd, mdcopts = muster_modelchem(input_model.model.method, input_model.driver)

        opts.update(mdcopts)

        ## Handle the basis set (ROBERT) the question is what value of k

        # Log the job settings (LORI)  Not sure if i need this
        logger.debug("JOB_OPTS")
        logger.debug(pp.pformat(opts))

        # Handle conversion from schema (flat key/value) keywords into local format
        optcmd = format_keywords(opts)
        madnessrec["commands"] = {}
        if mdccmd == "response":
            dft_cmds = optcmd.split(mdccmd)
            dft_cmds[1] = "response\n" + dft_cmds[1]

            madnessrec["infiles"]["moldft"] = {}
            madnessrec["infiles"]["moldft"]["input"] = dft_cmds[0] + molcmd
            madnessrec["infiles"]["molresponse"] = {}
            madnessrec["infiles"]["molresponse"]["rinput"] = dft_cmds[1]
            madnessrec["commands"]["moldft"] = [which("moldft")]
            madnessrec["commands"]["molresponse"] = [which("molresponse")]
        else:
            dft_cmds = optcmd
            madnessrec["infiles"]["moldft"] = {}
            madnessrec["infiles"]["moldft"]["input"] = dft_cmds + molcmd
            madnessrec["commands"]["moldft"] = [which("moldft")]

        print(dft_cmds)
        # optcmd="dft\n xc hf \nend\n"
        # print(madnessrec["infiles"]["input"])
        return madnessrec

    def execute(
            self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:
        num_commands = len(inputs["commands"])
        oexe = {}
        if num_commands == 2:
            success, dexe = execute(
                inputs["commands"]["moldft"],
                inputs["infiles"]["moldft"],
                scratch_exist_ok=True,
                scratch_name=inputs.get("scratch_name", None),
                scratch_directory=inputs["scratch_directory"],
                scratch_messy=True,
            )
            oexe["moldft"] = dexe
            success, dexe_response = execute(
                inputs["commands"]["molresponse"],
                inputs["infiles"]["molresponse"],
                scratch_messy=True,
                scratch_name=Path(dexe["scratch_directory"]).name,
                scratch_exist_ok=True,
            )
            oexe["molresponse"] = dexe_response
            print(dexe)
            print(dexe_response)
            return success, oexe
        else:
            print(inputs["commands"]["moldft"])
            success, dexe = execute(
                inputs["commands"]["moldft"],
                inputs["infiles"]["moldft"],
                scratch_exist_ok=True,
                scratch_name=inputs.get("scratch_name", None),
                scratch_directory=inputs["scratch_directory"],
                scratch_messy=True,
            )
            oexe["moldft"] = dexe
            return success, oexe

    def parse_output(
            self, outfiles: Dict[str, str], input_model: "AtomicInput"
    ) -> "AtomicResult":  # lgtm: [py/similar-function]

        # Get the stdout from the calculation (required)
        stdout = outfiles["moldft"]["stdout"]
        if "molresponse" in outfiles.keys():
            stdout += outfiles["molresponse"]["stdout"]

        # Read the MADNESj stdout file and, if needed, the hess or grad files
        qcvars, madhess, madgrad, madmol, version, errorTMP = harvest(input_model.molecule, **outfiles)
        ## pop the files because I think I need to
        outfiles.pop("moldft")
        if "molresponse" in outfiles.keys():
            outfiles.pop("molresponse")

        if madgrad is not None:
            qcvars["CURRENT GRADIENT"] = madgrad

        if madhess is not None:
            qcvars["CURRENT HESSIAN"] = madhess
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
