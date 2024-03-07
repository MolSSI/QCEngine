"""
Calls the Madness moldft executable.
"""

# import re
import copy
import json
import logging
import pprint
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple

import numpy as np

import qcelemental as qcel
from qcelemental.models import AtomicResult, Provenance, AtomicInput
from qcelemental.util import safe_version, which
from qcengine.config import TaskConfig, get_config
from qcengine.exceptions import UnknownError
from .germinate import muster_modelchem
from .harvester import extract_formatted_properties, harvest
from .keywords import format_keywords
from ..model import ProgramHarness
from ...exceptions import InputError
from ...util import execute, create_mpi_invocation, popen

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
logger = logging.getLogger(__name__)


class MadnessHarness(ProgramHarness):
    """
    Notes
    -----
    """

    _defaults = {
        "name": "madness",
        "scratch": True,
        "thread_safe": False,  # TODO: Check that MADNESS is thread safe?
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
             If both m-a-d-n-e-s-s and its harness dependency networkx are found, returns True.
             If raise_error is False and MADNESS, returns False.

        """

        # We only look for moldft but madness provides other QC functionality through more applications

        qc = which(
            "moldft",
            return_bool=False,
            raise_error=raise_error,
            raise_msg="Please install via https://github.com/m-a-d-n-e-s-s/madness",
        )
        return qc

    def get_version(self) -> str:
        self.found(raise_error=True)

        # Get the node configuration
        config = get_config()

        # Run MADNESS
        which_prog = which("moldft", return_bool=False)

        command = str(which_prog)

        if which_prog not in self.version_cache:
            with popen([which_prog, "--help"]) as exc:
                exc["proc"].wait(timeout=30)

            if exc["proc"].returncode != 0:
                raise UnknownError(exc["stderr"])
            else:
                for line in exc["stdout"].splitlines():
                    if "multiresolution suite" in line:
                        version = line.strip().split()[1]
                        self.version_cache[which_prog] = safe_version(version)
        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs madness in executable mode
        """
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)
        print("job_inputs", job_inputs)
        success, output = self.execute(job_inputs, extra_outfiles=["mad.calc_info.json"])

        print("output", output["stdout"])
        if "There is an error in the input file" in output["stdout"]:
            raise InputError(output["stdout"])
        if "not compiled" in output["stdout"]:
            # recoverable with a different compilation with optional modules
            raise InputError(output["stdout"])
        if success:
            num_commands = len(output)
            stdin = job_inputs["infiles"]["input"]
            output["outfiles"]["stdout"] = output["stdout"]
            output["outfiles"]["stderr"] = output["stderr"]
            output["outfiles"]["input"] = stdin
            return self.parse_output(output, input_model)
        else:

            raise UnknownError(output["stderr"])

    def build_input(
        self, input_model: AtomicInput, config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        #
        madnessrec = {
            "infiles": {},
            "scratch_directory": config.scratch_directory,
            "scratch_messy": config.scratch_messy,
        }

        # Prepare to write out the options
        opts = copy.deepcopy(input_model.keywords)
        opts = {k.lower(): v for k, v in opts.items()}
        print("opts: ", opts)

        # Determine the command to use to launch the code
        if config.use_mpiexec:
            madnessrec["command"] = create_mpi_invocation(which("moldft"), config)
            logger.info(f"Launching with mpiexec: {' '.join(madnessrec['command'])}")
        else:
            madnessrec["command"] = [which("nwchem")]

        # Handle Molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="madness", units="Bohr", return_data=True)
        print("moldata: ", moldata)
        print("molcmd: ", molcmd)
        molData = {}
        for k, v in moldata["keywords"].items():
            molData["dft__" + k] = v
        opts.update(molData)
        print("Method", input_model.model.method)
        mdccmd, mdcopts = muster_modelchem(input_model.model.method, input_model.driver)
        opts.update(mdcopts)

        logger.debug("JOB_OPTS")
        logger.debug(pp.pformat(opts))

        # Handle conversion from schema (flat key/value) keywords into local format
        optcmd = format_keywords(opts)
        # I need to split to geometry keywords and add it to the end of the geometry command in molcommand
        # if the geometry keyword exits
        if optcmd.find("geometry") != -1:
            geo_index = optcmd.find("geometry")  # find first occurrence of geometry
            end_index = optcmd[geo_index:].find("end")  # find first occurrence of end after geometry
            geometry_input = optcmd[
                geo_index + 8 : end_index + geo_index
            ]  # grab everything in between geometry and end

            optcmd = optcmd[0:geo_index] + optcmd[geo_index + end_index + 4 :]  # optcmd becomes everything else
            molcmd = molcmd.replace(
                "end", geometry_input.strip() + "\nend"
            )  # replace end with the added geometry input

        madnessrec["command"] = {}
        dft_cmds = optcmd
        madnessrec["infiles"] = {}
        madnessrec["infiles"]["input"] = dft_cmds + molcmd
        madnessrec["command"] = [which("moldft")]

        print("madness_rec", madnessrec)
        return madnessrec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:
        oexe = {}
        success, dexe = execute(
            command=inputs["command"],
            infiles=inputs["infiles"],
            outfiles=["mad.calc_info.json", "mad.scf_info.json"],
            scratch_exist_ok=True,
            scratch_name=inputs.get("scratch_name", None),
            scratch_directory=inputs["scratch_directory"],
            scratch_messy=True,
        )
        print("success", success)
        print("dexe", dexe)

        return success, dexe

    def parse_output(self, outfiles, input_model: "AtomicInput") -> "AtomicResult":  # lgtm: [py/similar-function]

        qcvars, madhess, madgrad, madmol, version, errorTMP = harvest(input_model.molecule, outfiles)

        moldft_out = outfiles

        m_stdout = moldft_out.pop("stdout")
        m_stderr = moldft_out.pop("stderr")
        print("after pop", moldft_out["outfiles"].keys())

        response_out = None
        r_stdout = None
        r_stderr = None
        if "molresponse" in outfiles.keys():
            response_out = outfiles.pop("molresponse")
            r_stdout = response_out.pop("stdout")
            r_stderr = response_out.pop("stderr")

        stdout = m_stdout
        if r_stdout is not None:
            stdout.update(r_stdout)

        if madgrad is not None:
            qcvars["CURRENT GRADIENT"] = madgrad

        if madhess is not None:
            qcvars["CURRENT HESSIAN"] = madhess
        # Normalize the output as a float or list of floats
        if input_model.driver.upper() == "PROPERTIES":
            retres = qcvars[f"RETURN_ENERGY"]
        else:
            retres = qcvars["RETURN_ENERGY"]

        if isinstance(retres, Decimal):
            retres = float(retres)
        elif isinstance(retres, np.ndarray):
            retres = retres.tolist()

        # Get the formatted properties
        qcprops = extract_formatted_properties(qcvars)
        # Format them inout an output
        m_native_files = {k: v for k, v in moldft_out["outfiles"].items() if v is not None}
        native_files = {
            "input": m_native_files["input"],
            "calc_info": m_native_files["mad.calc_info.json"],
            "scf_info": m_native_files["mad.scf_info.json"],
        }
        output_data = {
            "schema_name": "qcschema_output",
            "schema_version": 1,
            "extras": {"outfiles": outfiles, **input_model.extras},
            "native_files": native_files,
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
        output_data["extras"]["outfiles"] = {
            "input": native_files["input"],
            "calc_info": json.loads(native_files["calc_info"]),
            "scf_info": json.loads(native_files["scf_info"]),
        }

        output_data["success"] = True
        return AtomicResult(**{**input_model.dict(), **output_data})
