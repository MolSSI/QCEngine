import abc
import logging
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel
from pydantic.fields import Field
from qcelemental.models import AtomicInput, AtomicResult

from qcengine.exceptions import KnownErrorException

logger = logging.getLogger(__name__)


class ProgramHarness(BaseModel, abc.ABC):
    """Base class for objects that wrap calling external codes.

    Each subclass should implement how to find whether the code can be found,
    how to execute a calculation according to an AtomicInput,
    and how to parse the results."""

    _defaults: Dict[str, Any] = {}
    name: str
    scratch: bool
    thread_safe: bool
    thread_parallel: bool
    node_parallel: bool
    managed_memory: bool
    extras: Optional[Dict[str, Any]]

    class Config:
        allow_mutation = False
        extra = "forbid"
        arbitrary_types_allowed = True

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @abc.abstractmethod
    def compute(self, input_data: AtomicInput, config: "TaskConfig") -> AtomicResult:
        """Perform a computation given specification

        Parameters
        ----------
        input_data: AtomicInput
            Input specification for a calculation
        config: TaskConfig
            Resource configuration details
        observed_errors: List[str]
            List of errors that have been observed previously

        Returns
        -------
        AtomicResult
            Result of the calculation
        """

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

    def get_version(self) -> str:
        """Finds program, extracts version, returns normalized version string.

        Returns
        -------
        str
            Return a valid, safe python version string.
        """

    ## Computers

    def build_input(
        self, input_model: AtomicInput, config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Assemble all of the input settings needed to execute a computation given a specification

        Parameters
        ----------
        input_model: AtomicInput
            Input specification for a calculation
        config: TaskConfig
            Resource configuration details
        template: str
            Template for the calculation

        Returns
        -------
        Dict[str, Any]
            Any information needed by :meth:`execute` to perform the calculation
        """
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


class ErrorCorrectionProgramHarness(ProgramHarness):
    """Bases class for Harnesses that include logic to correct common errors

    **Implementation Guide**

    Each subclass should implement the ``_compute`` and ``build_input`` methods.

    These method should raise :class:`KnownErrorException`s when a task fails
    due to an error that your harness will be able to correct.
    The ``compute`` method for this class catches that exception,
    determines whether a correction is permitted based on the task configuration,
    and re-submits the computation with the details of the observed error.
    Your implementations of ``_build_input`` or ``_compute`` should alter
    their execution of the harnessed code based on the error details.

    The list of known errors for each class must be defined in the ``known_errors``
    class attribute of the harness.
    We recommend splitting error detection logic into many separate subclasses
    to make implementing the ``detect`` method for each error and any error
    mitigation logic simpler.
    """

    known_errors: List[KnownErrorException.__class__] = Field(
        None, description="List of classes that can recognize certain correctable errors from this Harness"
    )

    def compute(self, input_data: AtomicInput, config: "TaskConfig") -> AtomicResult:
        # Get the error correction configuration
        error_policy = input_data.protocols.error_correction

        # Run the method and, if it fails, assess if the failure is restartable
        observed_errors = {}  # Errors that have been observed previously
        while True:
            try:
                result = self._compute(input_data, config, observed_errors)
                break
            except KnownErrorException as e:
                logger.info(f"Caught a {type(e)} error.")

                # Determine whether this specific type of error is allowed
                correction_allowed = error_policy.allows(e.error_name)
                logger.info(f'Error correction for "{e.error_name}" is allowed')

                # Check if it has run before
                # TODO (wardlt): Should we allow errors to be run >1 time?
                previously_run = e.error_name in observed_errors

                # Run only if allowed and this error has not been seen before
                if correction_allowed and not previously_run:
                    logger.info("Adding error and restarting")
                    observed_errors[e.error_name] = e.details
                else:
                    logger.info(
                        "Error has been observed before and mitigation did not fix the issue. Raising exception"
                    )
                    raise e

        # Add the errors observed and corrected for, if any
        if len(observed_errors) > 0:
            result.extras["observed_errors"] = observed_errors
        return result

    @abc.abstractmethod
    def _compute(self, input_data: AtomicInput, config: "TaskConfig", observed_errors: Dict[str, dict]) -> AtomicResult:
        """Private method describing how to perform a computation given specification
        that includes a list of previously-observed errors

        Parameters
        ----------
        input_data: AtomicInput
            Input specification for a calculation
        config: TaskConfig
            Resource configuration details
        observed_errors: Dict[str, dict]
            Errors that have been observed and any details that go along with them

        Returns
        -------
        AtomicResult
            Result of the calculation
        """
        # TODO (wardlt): Should there be a default implementation of _compute: "build_inputs", "execute", parse_output"?

    @abc.abstractmethod
    def build_input(
        self,
        input_model: AtomicInput,
        config: "TaskConfig",
        template: Optional[str] = None,
        observed_errors: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Assemble all of the input settings needed to execute a computation given a specification
        that includes a list of previously-observed errors

        Parameters
        ----------
        input_model: AtomicInput
            Input specification for a calculation
        config: TaskConfig
            Resource configuration details
        observed_errors: List[str]
            List of errors that have been observed previously
        template: str
            Template for input file

        Returns
        -------
        Dict[str, Any]
            Any information needed by :meth:`execute` to perform the calculation
        """
