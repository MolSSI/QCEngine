from typing import TYPE_CHECKING, Any, ClassVar, Dict, Tuple, Union

from qcelemental.util import safe_version, which_import

from ..util import model_wrapper
from .model import ProcedureHarness

if TYPE_CHECKING:
    from qcmanybody.models import ManyBodyInput, ManyBodyResult

    from ..config import TaskConfig


class QCManyBodyProcedure(ProcedureHarness):

    _defaults: ClassVar[Dict[str, Any]] = {"name": "QCManyBody", "procedure": "manybody"}

    version_cache: Dict[str, str] = {}

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "qcmanybody",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install qcmanybody -c conda-forge`.",
        )

    def build_input_model(
        self, data: Union[Dict[str, Any], "ManyBodyInput"], *, return_input_schema_version: bool = False
    ) -> "ManyBodyInput":
        from qcmanybody.models import ManyBodyInput

        return self._build_model(data, "ManyBodyInput", return_input_schema_version=return_input_schema_version)

    # temporary can't use procedure/model.py b/c different import location than qcel.models
    def _build_model(
        self, data: Dict[str, Any], model: "BaseModel", /, *, return_input_schema_version: bool = False
    ) -> Union["BaseModel", Tuple["BaseModel", int]]:
        """
        Quick wrapper around util.model_wrapper for inherited classes
        """
        import qcmanybody

        v1_model = getattr(qcmanybody.models, model)
        v2_model = None

        if isinstance(data, v1_model):
            mdl = model_wrapper(data, v1_model)
        elif isinstance(data, v2_model):
            mdl = model_wrapper(data, v2_model)
        elif isinstance(data, dict):
            # remember these are user-provided dictionaries, so they'll have the mandatory fields,
            #   like driver, not the helpful discriminator fields like schema_version.
            # so long as versions distinguishable by a *required* field, id by dict is reliable.

            # if data.get("specification", False) or data.get("schema_version") == 2:
            #     mdl = model_wrapper(data, v2_model)
            # else:
            mdl = model_wrapper(data, v1_model)

        input_schema_version = mdl.schema_version
        if input_schema_version != 1:
            raise InputError("Can't use v2 ManyBody")

        if return_input_schema_version:
            return mdl, input_schema_version
            # return mdl.convert_v(2), input_schema_version
        else:
            return mdl
            # return mdl.convert_v(2)

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
