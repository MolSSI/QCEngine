from typing import TYPE_CHECKING, Any, Dict, Union

from qcelemental.util import safe_version, which_import

from .model import ProcedureHarness

if TYPE_CHECKING:
    from ..config import TaskConfig
    from qcmanybody.models import ManyBodyInput, ManyBodyResult


class QCManyBodyProcedure(ProcedureHarness):

    # v2: ClassVar[Dict[str, Any]]
    _defaults: Dict[str, Any] = {"name": "QCManyBody", "procedure": "manybody"}

    version_cache: Dict[str, str] = {}

    class Config(ProcedureHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "qcmanybody",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install qcmanybody -c conda-forge`.",
        )

    def build_input_model(self, data: Union[Dict[str, Any], "ManyBodyInput"]) -> "ManyBodyInput":
        from qcmanybody.models import ManyBodyInput

        return self._build_model(data, ManyBodyInput)

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("qcmanybody")
        if which_prog not in self.version_cache:
            import qcmanybody

            self.version_cache[which_prog] = safe_version(qcmanybody.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_model: "ManyBodyInput", config: "TaskConfig") -> "ManyBodyResult":
        from qcmanybody import ManyBodyComputer

        output_model = ManyBodyComputer.from_manybodyinput(input_model)

        return output_model
