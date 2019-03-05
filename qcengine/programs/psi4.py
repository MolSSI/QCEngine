"""
Calls the Psi4 executable.
"""

from qcelemental.models import FailedOperation, Result

from .executor import ProgramExecutor


class Psi4Executor(ProgramExecutor):

    _defaults = {
        "name": "Psi4",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }

    class Config(ProgramExecutor.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs Psi4 in API mode
        """

        try:
            import psi4
        except ModuleNotFoundError:
            raise ModuleNotFoundError("Could not find Psi4 in the Python path.")

        # Setup the job
        input_model = input_model.copy().dict()
        input_model["nthreads"] = config.ncores
        input_model["memory"] = int(config.memory * 1024 * 1024 * 1024 * 0.95)  # Memory in bytes
        input_model["success"] = False
        input_model["return_output"] = True

        if input_model["schema_name"] == "qcschema_input":
            input_model["schema_name"] = "qc_schema_input"

        scratch = config.scratch_directory
        if scratch is not None:
            input_model["scratch_location"] = scratch

        psi_version = self.parse_version(psi4.__version__)

        if psi_version > self.parse_version("1.2"):

            mol = psi4.core.Molecule.from_schema(input_model)
            if (mol.multiplicity() != 1) and ("reference" not in input_model["keywords"]):
                input_model["keywords"]["reference"] = "uhf"

            output_data = psi4.json_wrapper.run_json(input_model)
            if "extras" not in output_data:
                output_data["extras"] = {}

            output_data["extras"]["local_qcvars"] = output_data.pop("psi4:qcvars", None)
            if output_data["success"] is False:
                output_data["error"] = {"error_type": "internal_error", "error_message": output_data["error"]}

        else:
            raise TypeError("Psi4 version '{}' not understood.".format(psi_version))

        # Reset the schema if required
        output_data["schema_name"] = "qcschema_output"

        # Dispatch errors, PSIO Errors are not recoverable for future runs
        if output_data["success"] is False:

            if "PSIO Error" in output_data["error"]:
                raise ValueError(output_data["error"])

        # Move several pieces up a level
        if output_data["success"]:
            output_data["provenance"]["memory"] = round(output_data.pop("memory") / (1024**3), 3)  # Move back to GB
            output_data["provenance"]["nthreads"] = output_data.pop("nthreads")
            output_data["stdout"] = output_data.pop("raw_output", None)

            # Delete keys
            output_data.pop("return_ouput", None)

            return Result(**output_data)
        else:
            return FailedOperation(
                success=output_data.pop("success", False), error=output_data.pop("error"), input_data=output_data)

    def found(self) -> bool:
        try:
            import psi4
            return True
        except ModuleNotFoundError:
            return False
