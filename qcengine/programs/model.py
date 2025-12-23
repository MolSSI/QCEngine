import abc
import logging
from typing import Any, ClassVar, Dict, List, Optional, Tuple, Union

from pydantic import BaseModel, ConfigDict

from qcengine.config import TaskConfig
from qcengine.exceptions import KnownErrorException

from ..util import model_wrapper

logger = logging.getLogger(__name__)


class ProgramHarness(BaseModel, abc.ABC):
    """Base class for analytic single-geometry capable harnesses."""

    _defaults: ClassVar[Dict[str, Any]] = {}
    name: str
    scratch: bool
    thread_safe: bool
    thread_parallel: bool
    node_parallel: bool
    managed_memory: bool
    extras: Optional[Dict[str, Any]] = None

    model_config = ConfigDict(
        frozen=True,
        extra="forbid",
    )

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @abc.abstractmethod
    def compute(self, input_data: "AtomicInput", config: TaskConfig) -> Union["AtomicResult", "FailedOperation"]:
        """Top-level compute method to be implemented for every ProgramHarness

        Note:
            This method behave in any of the following ways:
                1. Return AtomicResult upon successful completion of a calculation
                2. Return FailedOperation object if an operation was unsuccessful or raised an exception. This is most
                    likely to occur if the underlying QC package has a QCSchema API that catches exceptions and
                    returns them as FailedOperation objects to end users.
                3. Raise an exception if a computation failed. The raised exception will be handled by the
                    qcng.compute() method and either raised or packaged as a FailedOperation object.
        """
        pass

    @staticmethod
    @abc.abstractmethod
    def found(raise_error: bool = False) -> bool:
        """
        Checks if the program can be found.

        Parameters
        ----------
        raise_error : bool, optional
            If True, raises an error if the program cannot be found.

        Returns
        -------
        bool
            Returns True if the program was found, False otherwise.
        """

    ## Utility

    def build_input_model(
        self, data: Dict[str, Any], *, return_input_schema_version: bool = False
    ) -> Union["AtomicInput", Tuple["AtomicInput", int]]:
        """
        Quick wrapper around util.model_wrapper for inherited classes
        """
        # Note: Someday when the multiple QCSchema versions QCEngine supports are all within the
        #  Pydantic v2 API base class, this can use discriminated unions instead of logic.

        from qcelemental.models.v1 import AtomicInput as v1_model
        from qcelemental.models.v2 import AtomicInput as v2_model

        if isinstance(data, v1_model):
            mdl = model_wrapper(data, v1_model)
        elif isinstance(data, v2_model):
            mdl = model_wrapper(data, v2_model)
        elif isinstance(data, dict):
            # remember these are user-provided dictionaries, so they'll have the mandatory fields,
            #   like driver, not the helpful discriminator fields like schema_version.
            # when dictionaries look the same, we can't correctly identify the user schema version
            #   so have to default to one or the other. Note that can force paths for testing by
            #   -1 -> 2 in schema_versions.
            # so long as versions distinguishable by a *required* field, id by dict is reliable.

            if data.get("specification", False) or data.get("schema_version") == 2:
                mdl = model_wrapper(data, v2_model)
            else:
                try:
                    mdl = model_wrapper(data, v1_model)
                except RuntimeError:
                    # form a QCSchema v1 layout model with Pydantic v2 API
                    # this is for py314 and only safe b/c immediately converted to v2 next
                    from qcelemental.models._v1v2 import AtomicInput as v1v2_model

                    mdl = model_wrapper(data, v1v2_model)

        input_schema_version = mdl.schema_version
        if return_input_schema_version:
            return mdl.convert_v(2), input_schema_version
        else:
            return mdl.convert_v(2)

    def get_version(self) -> str:
        """Finds program, extracts version, returns normalized version string.

        Returns
        -------
        str
            Return a valid, safe python version string.
        """

    ## Computers

    def build_input(
        self, input_model: "AtomicInput", config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        raise ValueError("build_input is not implemented for {}.", self.__class__)

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:
        raise ValueError("execute is not implemented for {}.", self.__class__)

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":
        raise ValueError("parse_output is not implemented for {}.", self.__class__)


class ErrorCorrectionProgramHarness(ProgramHarness, abc.ABC):
    """Base class for Harnesses that include logic to correct common errors

    Classes which implement this Harness must override the :meth:`_compute` method
    rather than :meth:`compute`. The ``compute`` method from this class will make
    calls to ``_compute`` with different actions as it attempts to correct errors in the
    input files.

    The error corrections are defined by first implementing a :class:`KnownErrorException`
    that contains logic to determine if a certain error occurs and returns an appropriate
    update to the keywords of an AtomicInput.

    Then, modify the ``_compute`` method of your Harness to run the ``detect_error`` method
    of each ``KnownErrorException`` for the Harness after it finishes performing the computation.
    The ``detect_error`` method will raise an exception that is caught by the
    ``ErrorCorrectionProgramHarness`` and used to determine if/how to re-run the computation.
    """

    def _compute(self, input_data: "AtomicInput", config: TaskConfig) -> "AtomicResult":
        raise NotImplementedError()

    def compute(self, input_data: "AtomicInput", config: TaskConfig) -> "AtomicResult":
        # Get the error correction configuration
        error_policy = input_data.specification.protocols.error_correction

        # Create a local copy of the input data
        local_input_data = input_data

        # Run the method and, if it fails, assess if the failure is restartable
        observed_errors = {}  # Errors that have been observed previously
        while True:
            try:
                result = self._compute(local_input_data, config)
                break
            except KnownErrorException as e:
                logger.info(f"Caught a {type(e)} error.")

                # Determine whether this specific type of error is allowed
                correction_allowed = error_policy.allows(e.error_name)
                if not correction_allowed:
                    logger.info(f'Error correction for "{e.error_name}" is not allowed')
                    raise e
                logger.info(f'Error correction for "{e.error_name}" is allowed')

                # Check if it has run before
                # TODO (wardlt): Should we allow errors to be run >1 time?
                previously_run = e.error_name in observed_errors
                if previously_run:
                    logger.info(
                        "Error has been observed before and mitigation did not fix the issue. Raising exception"
                    )
                    raise e

                # Generate and apply the updated keywords
                keyword_updates = e.create_keyword_update(local_input_data)
                new_keywords = local_input_data.specification.keywords.copy()
                new_keywords.update(keyword_updates)
                local_input_data = input_data.__class__(
                    **local_input_data.model_dump(exclude={"specification"}),
                    specification={**local_input_data.specification.model_dump(), "keywords": new_keywords},
                )

                # Store the error details and mitigations employed
                observed_errors[e.error_name] = {"details": e.details, "keyword_updates": keyword_updates}

        # Add the errors observed and corrected for, if any
        if len(observed_errors) > 0:
            result.extras["observed_errors"] = observed_errors
        return result
