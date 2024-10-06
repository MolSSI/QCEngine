import abc
from typing import Any, Dict, Tuple, Union

import qcelemental
from pydantic import BaseModel, ConfigDict

from ..util import model_wrapper


class ProcedureHarness(BaseModel, abc.ABC):

    name: str
    procedure: str

    model_config = ConfigDict(
        frozen=True,
        extra="forbid",
    )

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @abc.abstractmethod
    def build_input_model(self, data: Union[Dict[str, Any], "BaseModel"], raise_error: bool = True) -> "BaseModel":
        """
        Build and validate the input model, passes if the data was a normal BaseModel input.

        Parameters
        ----------
        data : Union[Dict[str, Any], 'BaseModel']
            A data blob to construct the model from or the input model itself
        raise_error : bool, optional
            Raise an error or not if the operation failed.

        Returns
        -------
        BaseModel
            The input model for the procedure.
        """

    @abc.abstractmethod
    def compute(self, input_data: "BaseModel", config: "TaskConfig") -> "BaseModel":
        pass

    @abc.abstractmethod
    def found(self, raise_error: bool = False) -> bool:
        """
        Checks if the program can be found.

        Returns
        -------
        bool
            If the proceudre was found or not.
        """

    def _build_model(
        self, data: Dict[str, Any], model: "BaseModel", /, *, return_input_schema_version: bool = False
    ) -> Union["BaseModel", Tuple["BaseModel", int]]:
        """
        Quick wrapper around util.model_wrapper for inherited classes
        """

        v1_model = getattr(qcelemental.models.v1, model)
        v2_model = getattr(qcelemental.models.v2, model)

        if isinstance(data, v1_model):
            mdl = model_wrapper(data, v1_model)
        elif isinstance(data, v2_model):
            mdl = model_wrapper(data, v2_model)
        elif isinstance(data, dict):
            # remember these are user-provided dictionaries, so they'll have the mandatory fields,
            #   like driver, not the helpful discriminator fields like schema_version.

            # for now, the two dictionaries look the same, so cast to the one we want
            # note that this prevents correctly identifying the user schema version when dict passed in, so either as_v1/None or as_v2 will fail
            mdl = model_wrapper(data, v1_model)  # TODO v2

        input_schema_version = mdl.schema_version
        if return_input_schema_version:
            return mdl.convert_v(1), input_schema_version  # non-psi4 return_dict=False fail w/o this
        else:
            return mdl.convert_v(1)

    def get_version(self) -> str:
        """Finds procedure, extracts version, returns normalized version string.
        Returns
        -------
        str
            Return a valid, safe python version string.
        """
