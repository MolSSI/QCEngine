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
from .harvester import extract_formatted_properties, tensor_to_numpy
from .keywords import format_keywords
from ..model import ProgramHarness
from ...exceptions import InputError
from ...util import execute, create_mpi_invocation, popen
from ..util import error_stamp


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
        which_prog = which("moldft")
        if config.use_mpiexec:
            which_prog = create_mpi_invocation(which_prog, config)
        else:
            command = [which_prog]

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

    def compute(self, input_model: AtomicInput, config: TaskConfig) -> "AtomicResult":
        """
        Runs madness in executable mode
        """
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)

        # Location resolution order config.scratch_dir, /tmp
        parent = config.scratch_directory
        scratch_messy = config.scratch_messy

        error_message = None
        compute_success = False
        # Now here we may have mulitple applications to run defined by the keys in the madnessrec dictionary

        extra_outfiles = {"moldft": ["mad.calc_info.json", "mad.scf_info.json"], "molresponse": ["response_base.json"]}

        all_output = {}
        app_succeed = []

        for app, job_inputs in job_inputs.items():
            output = {}

            print("app:", app)
            print("job_inputs:", job_inputs)
            print("____________________________________")

            outfiles = extra_outfiles[app]
            command = job_inputs["command"]
            infiles = job_inputs["infiles"]

            # print evertything
            print("command", command)
            print("infiles", infiles)
            print("outfiles", outfiles)
            print("parent", parent)
            print("scratch_name", app)
            print("scratch_suffix", "_madness")

            print("____________________________________")

            if app == "moldft":
                success, dexe = execute(
                    command=command,
                    infiles=infiles,
                    outfiles=outfiles,
                    scratch_directory=parent,
                    scratch_name=app,
                    scratch_messy=True,
                    scratch_exist_ok=True,
                )
            else:
                tempdir = dexe["scratch_directory"]
                success, dexe = execute(
                    command=command,
                    infiles=infiles,
                    outfiles=outfiles,
                    scratch_directory=tempdir,
                    scratch_name=app,
                    scratch_messy=True,
                    scratch_exist_ok=True,
                )
            print("____________________________________")

            print("success", success)
            print("dexe", dexe)

            stdin = job_inputs["infiles"]["input"]
            if "There is an error in the input file" in output["stdout"]:
                raise InputError(error_stamp(stdin, output["stdout"], output["stderr"]))
            if "not compiled" in output["stdout"]:
                # recoverable with a different compilation with optional modules
                raise InputError(error_stamp(stdin, output["stdout"], output["stderr"]))

            if success:
                output["outfiles"]["stdout"] = output["stdout"]
                output["outfiles"]["stderr"] = output["stderr"]
                output["outfiles"]["input"] = stdin
                app_succeed.append(app)
                all_output[app] = output
            else:
                raise UnknownError(error_stamp(stdin, output["stdout"], output["stderr"]))
        return self.parse_output(all_output, input_model)

    def build_input(
        self, input_model: AtomicInput, config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:

        method = input_model.model.method.lower()
        run_type = input_model.driver
        #
        madnessrec = {}

        for exec, exec_keywords in input_model.keywords.items():

            exec_rec = {}

            # Prepare to write out the options
            opts = copy.deepcopy(exec_keywords)
            opts = {k.lower(): v for k, v in opts.items()}

            # Determine the command to use to launch the code
            if config.use_mpiexec:
                exec_rec["command"] = create_mpi_invocation(which(exec), config)
                logger.info(f"Launching with mpiexec: {' '.join(exec_rec['command'])}")
            else:
                exec_rec["command"] = [which(exec)]

            if exec == "moldft":
                # Handle Molecule
                molcmd, moldata = input_model.molecule.to_string(dtype="madness", units="Bohr", return_data=True)
                molData = {}
                for k, v in moldata["keywords"].items():
                    molData["dft__" + k] = v
                opts.update(molData)
                if run_type == "gradient":
                    opts["dft__derivatives"] = True
                elif run_type == "hessian":
                    Pass
                    # opts["dft__hessian"] = True
                elif run_type == "optimization":
                    opts["dft__gopt"] = True
            elif exec == "molresponse":
                opts["dft__save"] = True
                opts["response__archive"] = "../moldft/restartdata"
                opts["response__xc"] = method

            print(opts)

            logger.debug("JOB_OPTS")
            logger.debug(pp.pformat(opts))
            exec_commands = format_keywords(opts)
            exec_rec["infiles"] = {}
            if exec == "moldft":
                exec_rec["infiles"]["input"] = exec_commands + molcmd
            else:
                exec_rec["infiles"]["molresponse.in"] = exec_commands
            print(exec_rec["infiles"])

            madnessrec[exec] = exec_rec
        print(madnessrec)

        return madnessrec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:
        success, dexe = execute(
            command=inputs["command"],
            infiles=inputs["infiles"],
            outfiles=["mad.calc_info.json", "mad.scf_info.json", "response_base.json"],
            scratch_exist_ok=True,
            scratch_name=inputs.get("scratch_name", None),
            scratch_directory=inputs["scratch_directory"],
            scratch_messy=True,
        )
        print("success", success)
        print("dexe", dexe)

        return success, dexe

    def extract_properties(moldft_calcinfo: Dict[str, Any], scfinfo) -> Dict[str, Any]:
        """Translate MRChem output to QCSChema properties.

        Parameters
        ----------

        Returns
        -------
        """
        properties = {}

        properties["calcinfo_nmo"] = moldft_calcinfo["calcinfo_nmo"]
        properties["calcinfo_nalpha"] = moldft_calcinfo["calcinfo_nalpha"]
        properties["calcinfo_nbeta"] = moldft_calcinfo["calcinfo_nbeta"]
        properties["calcinfo_natom"] = moldft_calcinfo["calcinfo_natom"]
        properties["return_energy"] = moldft_calcinfo["return_energy"]

        properties["scf_dipole_moment"] = tensor_to_numpy(scfinfo["scf_dipole_moment"])

        properties["nuclear_repulsion_energy"] = moldft_calcinfo["e_nrep"][-1]
        properties["scf_xc_energy"] = moldft_calcinfo["e_xc"][-1]
        properties["scf_total_energy"] = moldft_calcinfo["e_tot"][-1]
        properties["scf_iterations"] = moldft_calcinfo["iterations"]
        properties["scf_dipole_moment"] = scfinfo = ["scf_dipole_moment"]

        properties["scf_gradient"] = tensor_to_numpy(scfinfo["gradient"])

        return properties

    def parse_output(self, outfiles, input_model: "AtomicInput") -> "AtomicResult":  # lgtm: [py/similar-function]

        output_data = {}
        for app, outfiles in outfiles.items():
            if app == "moldft":

                # Get the qcvars
                moldft_out = {}
                moldft_out["stdout"] = outfiles["stdout"]
                moldft_out["stderr"] = outfiles["stderr"]

                calcinfo = json.loads(outfiles["outfiles"]["calc_info"])
                scfinfo = json.loads(outfiles["outfiles"]["scf_info"])

                output_data["properties"] = extract_properties(calcinfo, scfinfo)

                return_result = scfinfo["scf_energy"]
                if isinstance(return_result, Decimal):
                    retres = float(return_result)

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
                    k.upper(): str(v) if isinstance(v, Decimal) else v
                    for k, v in qcel.util.unnp(qcvars, flat=True).items()
                }
                output_data["extras"]["outfiles"] = {
                    "input": native_files["input"],
                    "calc_info": json.loads(native_files["calc_info"]),
                    "scf_info": json.loads(native_files["scf_info"]),
                }

                output_data["success"] = True
        return AtomicResult(**{**input_model.dict(), **output_data})
