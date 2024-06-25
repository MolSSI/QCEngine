from typing import TYPE_CHECKING, Any, Dict, Union

from qcelemental.models import OptimizationInput, OptimizationResult
from qcelemental.util import safe_version, which_import

from .model import ProcedureHarness

if TYPE_CHECKING:
    from qcengine.config import TaskConfig
    from qcmanybody.models.generalized_optimization import GeneralizedOptimizationInput, GeneralizedOptimizationResult


class OptKingProcedure(ProcedureHarness):

    _defaults = {"name": "OptKing", "procedure": "optimization"}

    version_cache: Dict[str, str] = {}

    class Config(ProcedureHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "optking",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install optking -c conda-forge`.",
        )

    def build_input_model(self, data: Union[Dict[str, Any], "OptimizationInput"]) -> "OptimizationInput":
        return self._build_model(data, OptimizationInput)

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


class GenOptKingProcedure(OptKingProcedure):

    # note that "procedure" value below not used
    _defaults = {"name": "GenOptKing", "procedure": "genoptimization"}

    version_cache: Dict[str, str] = {}

    def found(self, raise_error: bool = False) -> bool:
        qc = which_import(
            "optking",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install optking -c conda-forge`.",
        )
        dep = which_import(
            "qcmanybody",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="For GenOptKing harness, please install via `conda install qcmanybody -c conda-forge`.",
        )

        return qc and dep

    def build_input_model(
        self, data: Union[Dict[str, Any], "GeneralizedOptimizationInput"]
    ) -> "GeneralizedOptimizationInput":
        from qcmanybody.models.generalized_optimization import GeneralizedOptimizationInput

        return self._build_model(data, GeneralizedOptimizationInput)

    def compute(
        self, input_model: "GeneralizedOptimizationInput", config: "TaskConfig"
    ) -> "GeneralizedOptimizationResult":
        self.found(raise_error=True)

        import optking
        from qcmanybody.models.generalized_optimization import GeneralizedOptimizationResult

        input_data = input_model.dict()

        # Set retries to two if zero while respecting local_config
        local_config = config.dict()
        local_config["retries"] = local_config.get("retries", 2) or 2
        input_data["input_specification"]["extras"]["_qcengine_local_config"] = local_config

        # Run the program
        output_data = optking.optimize_qcengine(input_data)

        output_data["schema_name"] = "qcschema_generalizedoptimizationresult"
        output_data["input_specification"]["extras"].pop("_qcengine_local_config", None)
        if output_data["success"]:
            output_data = GeneralizedOptimizationResult(**output_data)

        return output_data
