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
from qcengine.programs.mrchem import extract_properties
from .harvester import extract_formatted_properties, tensor_to_numpy
from .keywords import format_keywords
from ..model import ProgramHarness
from ...exceptions import InputError
from ...util import execute, create_mpi_invocation, popen
from ..util import error_stamp


pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
logger = logging.getLogger(__name__)

extra_outfiles = {"moldft": ["mad.calc_info.json", "mad.scf_info.json"], "molresponse": ["response_base.json"]}
calc_info_json = "mad.calc_info.json"
scfinfo_json = "mad.scf_info.json"
response_json = "response_base.json"


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

            outfiles = extra_outfiles[app]
            command = job_inputs["command"]
            infiles = job_inputs["infiles"]

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
            output = {"stdout": dexe["stdout"], "stderr": dexe["stderr"], "outfiles": dexe["outfiles"]}

            stdin = job_inputs["infiles"]
            if "There is an error in the input file" in output["stdout"]:
                raise InputError(error_stamp(stdin, output["stdout"], output["stderr"]))
            if "not compiled" in output["stdout"]:
                # recoverable with a different compilation with optional modules
                raise InputError(error_stamp(stdin, output["stdout"], output["stderr"]))

            if success:
                output["outfiles"]["stdout"] = output["stdout"]
                output["outfiles"]["stderr"] = output["stderr"]
                output["outfiles"]["input"] = stdin
                output["success"] = success
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
                if len(input_model.keywords["moldft"]) == 0:
                    opts["dft__save"] = True
                # if the number of execs is greater than 1, we need to save moldft archive
                if len(input_model.keywords) > 1:
                    opts["dft__save"] = True
            elif exec == "molresponse":
                opts["response__archive"] = "../mad.restartdata"
                opts["response__xc"] = method

            logger.debug("JOB_OPTS")
            logger.debug(pp.pformat(opts))
            exec_commands = format_keywords(opts)
            exec_rec["infiles"] = {}
            if exec == "moldft":
                exec_rec["infiles"]["input"] = exec_commands + molcmd
            else:
                exec_rec["infiles"]["response.in"] = exec_commands

            madnessrec[exec] = exec_rec

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

        return success, dexe

    def extract_moldft_properties(self, moldft_calcinfo: Dict[str, Any], scfinfo: Dict[str, Any]) -> Dict[str, Any]:
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

        scf_data = moldft_calcinfo["scf_e_data"]

        properties["nuclear_repulsion_energy"] = scf_data["e_nrep"][-1]
        properties["scf_xc_energy"] = scf_data["e_xc"][-1]
        properties["scf_total_energy"] = scf_data["e_tot"][-1]
        properties["scf_iterations"] = scf_data["iterations"]

        extra_properties = {}
        extra_properties["cpu_time"] = moldft_calcinfo["time_tag"]["cpu_time"]
        extra_properties["wall_time"] = moldft_calcinfo["time_tag"]["wall_time"]
        extra_properties["iterations"] = scf_data["iterations"]

        # properties["scf_gradient"] = tensor_to_numpy(scfinfo["gradient"])

        return properties, extra_properties

    def extract_response_properties(self, response_json):

        extra_properties = {}
        polarizability = tensor_to_numpy(response_json["response_data"]["data"]["alpha"][0])
        extra_properties["polarizability"] = polarizability[-1, :].reshape(3, 3)

        extra_properties["omega"] = response_json["parameters"]["omega"]
        extra_properties["cpu_time"] = response_json["time_data"]["cpu_time"]
        extra_properties["wall_time"] = response_json["time_data"]["wall_time"]
        extra_properties["iterations"] = response_json["response_data"]["iterations"]

        properties = {}
        properties["scf_iterations"] = response_json["response_data"]["iterations"]

        return properties, extra_properties
        # prepare a list of computed response properties

        # fill up return_result

    def parse_output(self, outfiles, input_model: "AtomicInput") -> "AtomicResult":  # lgtm: [py/similar-function]

        output_data = {}

        driver = input_model.driver

        # things to collect
        # 1. If the application succeded
        # 2. Standard output and error for each application
        # 4. The native files for each application
        # 3. The properties for each application from the native files

        # Then based on the driver we will also collect either moldft energy or all properties
        all_success = {}
        all_std_out = {}
        all_std_err = {}
        all_native_files = {}
        all_properties = {}
        all_extras = {}
        total_iterations = []

        for app, output in outfiles.items():
            all_success[app] = output["success"]
            all_std_out[app] = output["stdout"]
            all_std_err[app] = output["stderr"]
            all_native_files[app] = output["outfiles"]

            # collect native files from outfiles
            if app == "moldft":
                moldft_calcinfo = json.loads(output["outfiles"][calc_info_json])
                scfinfo = json.loads(output["outfiles"][scfinfo_json])
                all_properties[app], all_extras[app] = self.extract_moldft_properties(moldft_calcinfo, scfinfo)
            elif app == "molresponse":
                response_base_file = json.loads(output["outfiles"][response_json])
                (
                    all_properties[app],
                    extra_properties,
                ) = self.extract_response_properties(response_base_file)
                all_extras[app] = extra_properties
                all_success[app] = response_base_file["converged"]

            total_iterations.append(all_properties[app]["scf_iterations"])
            del all_properties[app]["scf_iterations"]
        output_data["success"] = True if np.array([v for v in all_success.values()]).all() else False
        output_data["stdout"] = " ".join(["\n"] + [v for v in all_std_out.values()])
        output_data["stderr"] = " ".join(["\n"] + [v for v in all_std_err.values()])
        output_data["native_files"] = all_native_files
        output_data["extras"] = all_extras

        #
        all_prop = {}
        for app, prop in all_properties.items():
            all_prop.update(prop)
        # first sum up total interations for each application
        total_scf_iter = sum(total_iterations)
        all_prop["scf_iterations"] = total_scf_iter
        output_data["properties"] = all_prop

        if input_model.driver == "energy":
            output_data["return_result"] = output_data["properties"]["return_energy"]
        elif input_model.driver == "properties":
            output_data["return_result"] = output_data["properties"]
        else:
            raise InputError(f"Driver {input_model.driver} not implemented for MRChem.")
        # Get the qcvars

        return AtomicResult(**{**input_model.dict(), **output_data})
