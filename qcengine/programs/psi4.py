"""
Calls the Psi4 executable.
"""
import json
import os
from typing import Dict

from qcelemental.models import AtomicResult
from qcelemental.util import deserialize, parse_version, safe_version, which

from ..exceptions import InputError, RandomError, ResourceError, UnknownError
from ..util import execute, popen, temporary_directory
from .model import ProgramHarness


class Psi4Harness(ProgramHarness):

    _defaults = {
        "name": "Psi4",
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
            "psi4",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install psi4 -c psi4`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("psi4")
        if which_prog not in self.version_cache:
            with popen([which_prog, "--version"]) as exc:
                exc["proc"].wait(timeout=30)
            self.version_cache[which_prog] = safe_version(exc["stdout"])

        candidate_version = self.version_cache[which_prog]

        if "undef" in candidate_version:
            raise TypeError(
                "Using custom build without tags. Please pull git tags with `git pull origin master --tags`."
            )

        return candidate_version

    def compute(self, input_model: "AtomicInput", config: "JobConfig") -> "AtomicResult":
        """
        Runs Psi4 in API mode
        """
        self.found(raise_error=True)
        pversion = parse_version(self.get_version())

        if pversion < parse_version("1.2"):
            raise ResourceError("Psi4 version '{}' not understood.".format(self.get_version()))

        # Location resolution order config.scratch_dir, $PSI_SCRATCH, /tmp
        parent = config.scratch_directory
        if parent is None:
            parent = os.environ.get("PSI_SCRATCH", None)

        error_type = None
        error_message = None
        compute_success = False

        with temporary_directory(parent=parent, suffix="_psi_scratch") as tmpdir:

            caseless_keywords = {k.lower(): v for k, v in input_model.keywords.items()}
            if (input_model.molecule.molecular_multiplicity != 1) and ("reference" not in caseless_keywords):
                input_model.keywords["reference"] = "uhf"

            # Old-style JSON-based command line
            if pversion < parse_version("1.4a2.dev160"):

                # Setup the job
                input_data = input_model.dict(encoding="json")
                input_data["nthreads"] = config.ncores
                input_data["memory"] = int(config.memory * 1024 * 1024 * 1024 * 0.95)  # Memory in bytes
                input_data["success"] = False
                input_data["return_output"] = True

                if input_data["schema_name"] == "qcschema_input":
                    input_data["schema_name"] = "qc_schema_input"

                # Execute the program
                success, output = execute(
                    [which("psi4"), "--scratch", tmpdir, "--json", "data.json"],
                    {"data.json": json.dumps(input_data)},
                    ["data.json"],
                    scratch_directory=tmpdir,
                )

                if success:
                    output_data = json.loads(output["outfiles"]["data.json"])
                    if "extras" not in output_data:
                        output_data["extras"] = {}

                    # Check QCVars
                    local_qcvars = output_data.pop("psi4:qcvars", None)
                    if local_qcvars:
                        # Edge case where we might already have qcvars, should not happen
                        if "qcvars" in output_data["extras"]:
                            output_data["extras"]["local_qcvars"] = local_qcvars
                        else:
                            output_data["extras"]["qcvars"] = local_qcvars

                    if output_data["success"] is False:
                        if "error_message" not in output_data["error"]:
                            error_message = output_data["error"]
                            error_type = "internal_error"
                        else:
                            error_message = output_data["error"]["error_message"]
                            error_type = output_data["error"]["error_type"]

                    else:
                        compute_success = True

                else:
                    error_message = output["stderr"]
                    error_type = "execution_error"

                # Reset the schema if required
                output_data["schema_name"] = "qcschema_output"
                output_data.pop("memory", None)
                output_data.pop("nthreads", None)
                output_data["stdout"] = output_data.pop("raw_output", None)

            else:

                if ("psiapi" in input_model.extras) and input_model.extras["psiapi"]:
                    import psi4

                    psi4.core.set_num_threads(config.ncores, quiet=True)
                    psi4.set_memory(f"{config.memory}GB", quiet=True)
                    output_data = psi4.schema_wrapper.run_qcschema(input_model).dict()
                    output_data["extras"]["psiapi_evaluated"] = True
                    success = True
                else:
                    run_cmd = [
                        which("psi4"),
                        "--scratch",
                        tmpdir,
                        "--nthread",
                        str(config.ncores),
                        "--memory",
                        f"{config.memory}GB",
                        "--qcschema",
                        "data.msgpack",
                    ]
                    input_files = {"data.msgpack": input_model.serialize("msgpack-ext")}
                    success, output = execute(
                        run_cmd, input_files, ["data.msgpack"], as_binary=["data.msgpack"], scratch_directory=tmpdir
                    )
                    if success:
                        output_data = deserialize(output["outfiles"]["data.msgpack"], "msgpack-ext")

                if success:
                    if output_data["success"] is False:
                        error_message = output_data["error"]["error_message"]
                        error_type = output_data["error"]["error_type"]
                    else:
                        compute_success = True
                else:
                    error_message = output["stderr"]
                    error_type = "execution_error"

        # Dispatch errors, PSIO Errors are not recoverable for future runs
        if compute_success is False:

            if "PSIO Error" in error_message:
                if "scratch directory" in error_message:
                    # Psi4 cannot access the folder or file
                    raise ResourceError(error_message)
                else:
                    # Likely a random error, worth retrying
                    raise RandomError(error_message)
            elif ("SIGSEV" in error_message) or ("SIGSEGV" in error_message) or ("segmentation fault" in error_message):
                raise RandomError(error_message)
            elif ("TypeError: set_global_option" in error_message) or (error_type == "ValidationError"):
                raise InputError(error_message)
            elif "RHF reference is only for singlets" in error_message:
                raise InputError(error_message)
            else:
                raise UnknownError(error_message)

        # Move several pieces up a level
        output_data["provenance"]["memory"] = round(config.memory, 3)
        output_data["provenance"]["nthreads"] = config.ncores

        # Delete keys
        output_data.pop("return_output", None)

        return AtomicResult(**output_data)
