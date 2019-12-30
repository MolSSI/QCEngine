from typing import Any, Dict, Union

from qcelemental.models import OptimizationResult, OptimizationInput
from qcelemental.util import which_import

from .model import ProcedureHarness


class OptKingProcedure(ProcedureHarness):

    _defaults = {"name": "OptKing", "procedure": "optimization"}

    class Config(ProcedureHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "optking",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install optking -c psi4`.",
        )

    def build_input_model(self, data: Union[Dict[str, Any], "OptimizationInput"]) -> "OptimizationInput":
        return self._build_model(data, OptimizationInput)

    def compute(self, input_model: "OptimizationInput", config: "JobConfig") -> "Optimization":
        try:
            import optking
        except ModuleNotFoundError:
            raise ModuleNotFoundError("Could not find optking in the Python path.")

        input_data = input_model.dict()

        # Set retries to two if zero while respecting local_config
        local_config = config.dict()
        local_config["retries"] = local_config.get("retries", 2) or 2
        input_data["input_specification"]["extras"]["_qcengine_local_config"] = local_config

        # Run the program
        output_data = optking.optwrapper.optimize_qcengine(input_data)

        output_data["schema_name"] = "qcschema_optimization_output"
        output_data["input_specification"]["extras"].pop("_qcengine_local_config", None)
        if output_data["success"]:
            output_data = OptimizationResult(**output_data)

        return output_data
