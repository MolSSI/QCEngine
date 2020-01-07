import abc
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel
from qcelemental.models import AtomicInput, AtomicResult

from qcengine.exceptions import KnownError


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
        allow_mutation: False
        extra: "forbid"

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

    Errors from the program being executed that are ``KnownError``s can be corrected.

    Each subclass should implement the ``_compute`` and ``build_input`` methods.
    """

    def compute(self, input_data: AtomicInput, config: "TaskConfig") -> AtomicResult:
        # TODO (wardlt): Add a protocol for controlling the error correction
        observed_errors = []  # List of the errors that have been caught
        while True:
            try:
                result = self._compute(input_data, config, observed_errors)
                break
            except KnownError as e:
                if len(observed_errors) > 2:
                    raise e
                else:
                    observed_errors.append(e.error_name)

        # Add the errors observed and corrected for, if any
        if len(observed_errors) > 0:
            result.extras["observed_errors"] = observed_errors
        return result

    @abc.abstractmethod
    def _compute(self, input_data: AtomicInput, config: "TaskConfig", observed_errors: List[str]) -> AtomicResult:
        """Private method describing how to perform a computation given specification
        that includes a list of previously-observed errors

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
