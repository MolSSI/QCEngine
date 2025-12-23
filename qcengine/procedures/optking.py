import logging
import sys
from io import StringIO
from typing import TYPE_CHECKING, Any, ClassVar, Dict, Union

from qcelemental.models.v1 import OptimizationResult
from qcelemental.models.v2 import OptimizationInput
from qcelemental.util import safe_version, which_import

from .model import ProcedureHarness

if TYPE_CHECKING:
    from ..config import TaskConfig


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

    def compute(self, input_model: "OptimizationInput", config: "TaskConfig") -> "Optimization":
        if self.found(raise_error=True):
            import optking

        log_stream = StringIO()
        logname = "psi4.optking" if "psi4" in sys.modules else "optking"
        log = logging.getLogger(logname)
        log.addHandler(logging.StreamHandler(log_stream))
        log.setLevel("INFO")

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
            output = output_v1  # TODO almost certainly wrong -- need v2 conversion?

        return output
