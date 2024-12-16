from typing import Any, ClassVar, Dict, Union

from qcelemental.models.v1 import OptimizationResult
from qcelemental.models.v2 import OptimizationInput
from qcelemental.util import safe_version, which_import

from .model import ProcedureHarness


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

        input_data_v1 = input_model.convert_v(1).dict()

        # Temporary patch for geomeTRIC
        input_data_v1["initial_molecule"]["symbols"] = list(input_data_v1["initial_molecule"]["symbols"])

        # Set retries to two if zero while respecting local_config
        local_config = config.dict()
        local_config["retries"] = local_config.get("retries", 2) or 2
        input_data_v1["input_specification"]["extras"]["_qcengine_local_config"] = local_config

        # Run the program
        output_v1 = geometric.run_json.geometric_run_json(input_data_v1)

        output_v1["provenance"] = {
            "creator": "geomeTRIC",
            "routine": "geometric.run_json.geometric_run_json",
            "version": geometric.__version__,
        }

        output_v1["schema_name"] = "qcschema_optimization_output"  # overwrites OptIn value
        output_v1["input_specification"]["extras"].pop("_qcengine_local_config", None)
        if output_v1["success"]:
            output_v1 = OptimizationResult(**output_v1)
            output = output_v1.convert_v(2, external_input_data=input_model)
        else:
            output = output_v1  # TODO almost certainly wrong, needs v2 conv

        return output
