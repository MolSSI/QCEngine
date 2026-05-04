from typing import Any, ClassVar, Dict, Union

from qcelemental.models._v1v2 import OptimizationResult
from qcelemental.models.v2 import OptimizationInput
from qcelemental.util import safe_version, which_import

from ..exceptions import RandomError, UnknownError
from ..util import QCEL_V1V2_SHIM_CODE, environ_context
from .model import ProcedureHarness


def _handle_errors(output_data: Dict[str, Any], default_error: str) -> tuple[str, str]:
    error_block = output_data.get("error")
    if isinstance(error_block, dict):
        if "error_message" in error_block:
            return error_block["error_message"], error_block.get("error_type", "unknown_error")
        return str(error_block), "unknown_error"
    if isinstance(error_block, str):
        return error_block, "unknown_error"
    return default_error, "unknown_error"


def _attach_retry_provenance(failed_result: Dict[str, Any], retries: int) -> None:
    trajectory = failed_result.get("trajectory")
    if isinstance(trajectory, list):
        for point in trajectory:
            if isinstance(point, dict):
                point.setdefault("provenance", {}).setdefault("retries", retries)


class GeometricProcedure(ProcedureHarness):

    _defaults: ClassVar[Dict[str, Any]] = {"name": "geomeTRIC", "procedure": "optimization"}

    version_cache: Dict[str, str] = {}

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "geometric",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install geometric -c conda-forge`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("geometric")
        if which_prog not in self.version_cache:
            import geometric

            self.version_cache[which_prog] = safe_version(geometric.__version__)

        return self.version_cache[which_prog]

    def build_input_model(
        self, data: Union[Dict[str, Any], "OptimizationInput"], *, return_input_schema_version: bool = False
    ) -> "OptimizationInput":
        return self._build_model(data, "OptimizationInput", return_input_schema_version=return_input_schema_version)

    def compute(self, input_model: "OptimizationInput", config: "TaskConfig") -> "OptimizationResult":
        try:
            import geometric
        except ModuleNotFoundError:
            raise ModuleNotFoundError("Could not find geomeTRIC in the Python path.")

        input_data_v2 = input_model.model_dump()

        # Temporary patch for geomeTRIC
        input_data_v2["initial_molecule"]["symbols"] = list(input_data_v2["initial_molecule"]["symbols"])

        # Set retries to two if zero while respecting local_config
        local_config = config.model_dump()
        local_config["retries"] = local_config.get("retries", 2) or 2
        input_data_v2["specification"]["specification"]["extras"]["_qcengine_local_config"] = local_config

        # Run the program
        # * geomeTRIC only speaks v1 as of v1.1
        # * envvar allows geomeTRIC to call back QCEngine for gradients
        input_data__v1v2 = OptimizationInput(**input_data_v2).convert_v(QCEL_V1V2_SHIM_CODE).model_dump()
        with environ_context(env={"QCNG_USE_V1V2_SHIM": "1"}):
            output_v1 = geometric.run_json.geometric_run_json(input_data__v1v2)

        output_v1["provenance"] = {
            "creator": "geomeTRIC",
            "routine": "geometric.run_json.geometric_run_json",
            "version": geometric.__version__,
        }

        output_v1["schema_name"] = "qcschema_optimization_output"  # overwrites OptIn value
        output_v1["input_specification"]["extras"].pop("_qcengine_local_config", None)
        if output_v1["success"]:
            output_model__v1v2 = OptimizationResult(**output_v1)
            output = output_model__v1v2.convert_v(2, external_input_data=input_model)
        else:
            error_message, error_type = _handle_errors(output_v1, "geomeTRIC failed without an error message.")
            is_random_path = error_type == "random_error"
            if is_random_path:
                _attach_retry_provenance(output_v1, local_config["retries"])
            error_cls = RandomError if is_random_path else UnknownError
            raise error_cls(error_message, extras={"failed_result": output_v1})

        return output
