import logging
from io import StringIO
from typing import TYPE_CHECKING, Any, ClassVar, Dict, Union

from qcelemental.util import parse_version, safe_version, which_import

from ..exceptions import InputError, RandomError, ResourceError, UnknownError
from .model import ProcedureHarness

if TYPE_CHECKING:
    from ..config import TaskConfig


def _handle_errors(output_data: Dict[str, Any], default_error: str) -> tuple[str, str]:
    error_block = output_data.get("error")
    if isinstance(error_block, dict):
        if "error_message" in error_block:
            return error_block["error_message"], error_block.get("error_type", "unknown_error")
        return str(error_block), "unknown_error"
    if isinstance(error_block, str):
        return error_block, "unknown_error"
    return default_error, "unknown_error"


class OptKingProcedure(ProcedureHarness):

    _defaults: ClassVar[Dict[str, Any]] = {"name": "OptKing", "procedure": "optimization"}

    version_cache: Dict[str, str] = {}

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "optking",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install optking -c conda-forge`.",
        )

    def build_input_model(
        self, data: Union[Dict[str, Any], "OptimizationInput"], *, return_input_schema_version: bool = False
    ) -> "OptimizationInput":
        return self._build_model(data, "OptimizationInput", return_input_schema_version=return_input_schema_version)

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("optking")
        if which_prog not in self.version_cache:
            import optking

            self.version_cache[which_prog] = safe_version(optking.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_model: "OptimizationInput", config: "TaskConfig") -> "OptimizationResult":
        if self.found(raise_error=True):
            import optking

        log_stream = StringIO()
        logname = f"{optking.log_name}{optking.__name__}"
        log = logging.getLogger(logname)
        log.addHandler(logging.StreamHandler(log_stream))
        log.setLevel("INFO")

        if parse_version(self.get_version()) < parse_version("0.5"):
            from qcelemental.models.v1 import OptimizationResult

            input_data_v1 = input_model.convert_v(1).dict()

            # Set retries to two if zero while respecting local_config
            local_config = config.dict()
            local_config["retries"] = local_config.get("retries", 2) or 2
            input_data_v1["input_specification"]["extras"]["_qcengine_local_config"] = local_config

            # Run the program
            output_v1 = optking.optwrapper.optimize_qcengine(input_data_v1)
            output_v1["stdout"] = log_stream.getvalue()

            output_v1["input_specification"]["extras"].pop("_qcengine_local_config", None)
            if output_v1["success"]:
                output_v1 = OptimizationResult(**output_v1)
                output = output_v1.convert_v(2, external_input_data=input_model)
            else:
                error_message, error_type = _handle_errors(output_v1, "OptKing failed without an error message.")
                error_cls = {
                    "random_error": RandomError,
                    "unknown_error": UnknownError,
                    "input_error": InputError,
                    "resource_error": ResourceError,
                }.get(error_type, UnknownError)
                raise error_cls(
                    error_message,
                    extras={
                        "failed_result": output_v1,
                    },
                )

        else:
            # optking v0.5 can run QCSchema 1<->1 and 2<->2, so update to 2
            from qcelemental.models.v2 import OptimizationResult

            input_data_v2 = input_model.model_dump()

            # Set retries to two if zero while respecting local_config
            local_config = config.model_dump()
            local_config["retries"] = local_config.get("retries", 2) or 2
            input_data_v2["specification"]["specification"]["extras"]["_qcengine_local_config"] = local_config

            # Run the program
            output_v2 = optking.optwrapper.optimize_qcengine(input_data_v2)
            output_v2["stdout"] = log_stream.getvalue()

            output_v2["input_data"]["specification"]["specification"]["extras"].pop("_qcengine_local_config", None)
            if output_v2["success"]:
                output = OptimizationResult(**output_v2)
            else:
                error_message, error_type = _handle_errors(output_v2, "OptKing failed without an error message.")
                error_cls = {
                    "random_error": RandomError,
                    "unknown_error": UnknownError,
                    "input_error": InputError,
                    "resource_error": ResourceError,
                }.get(error_type, UnknownError)
                raise error_cls(
                    error_message,
                    extras={
                        "failed_result": output_v2,
                    },
                )

        return output
