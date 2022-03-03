"""
Calls TeraChem in its "server mode" via a protobuf interface.
"""
import logging
from importlib import import_module
from os import getenv
from typing import TYPE_CHECKING, Any, Dict, Union

from qcelemental.models import AtomicResult, FailedOperation
from qcelemental.util import which_import

from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig

logger = logging.getLogger(__name__)

_pbs_defaults = {
    "name": "terachem_pbs",
    "scratch": False,
    "thread_safe": False,
    "thread_parallel": False,
    "node_parallel": False,
    "managed_memory": True,
}


class TeraChemPBSHarness(ProgramHarness):
    """QCEngine Harness for interfacing with the TeraChem running in Protocol Buffer Server Mode"""

    _defaults = _pbs_defaults
    _tcpb_package: str = "tcpb"
    _tcpb_min_version: str = "0.7.0"
    _tcpb_client: str = "TCProtobufClient"
    _env_vars: Dict[str, Any] = {
        "host": getenv("TERACHEM_PBS_HOST", "127.0.0.1"),
        "port": int(getenv("TERACHEM_PBS_PORT", 11111)),
    }
    _env_vars_external: str = "TERACHEM_PBS_HOST, TERACHEM_PBS_PORT"

    class Config(ProgramHarness.Config):
        pass

    @classmethod
    def found(cls, raise_error: bool = False) -> bool:
        """Whether TeraChemPBS harness is ready for operation.
        Parameters
        ----------
        raise_error: bool
            Passed on to control negative return between False and ModuleNotFoundError raised.
        Returns
        -------
        bool
            If tcpb package is found and server available, returns True.
            If raise_error is False and tcpb package missing and/or server us unavailable, returns False.
            If raise_error is True and tcpb package missing and/or server us unavailable, the error message is raised.
        """
        tcpb_pkg_available = which_import(
            cls._tcpb_package,
            return_bool=True,
            raise_error=raise_error,
            raise_msg=f"TeraChem protobuf client package (tcpb) not found. Please install tcpb>={cls._tcpb_min_version}.",
        )
        if not tcpb_pkg_available:
            return False

        tcpb = import_module(f"{cls._tcpb_package}")

        try:
            with getattr(tcpb, cls._tcpb_client)(**cls._env_vars) as client:
                return client.is_available()
        except tcpb.exceptions.ServerError as e:
            msg = (
                f"Unable to connect to TeraChem server at {cls._env_vars}. If you want to connect to a server elsewhere, consider "
                f"setting environment variables {cls._env_vars_external}."
            )
            logger.error(msg)
            if raise_error:
                raise OSError(msg) from e
        return False

    def get_version(self) -> str:
        """Returns version of TeraChem Protocol Buffer Server"""
        try:
            tcpb = import_module(self._tcpb_package)
        except ModuleNotFoundError:
            return None
        else:
            try:
                return tcpb.__version__
            except AttributeError:
                return None

    def compute(
        self, input_model: "AtomicInput", config: "TaskConfig" = None
    ) -> Union["AtomicResult", "FailedOperation"]:
        """
        Submit AtomicInput to TeraChem running in "server mode"
        """
        self.found(raise_error=True)

        tcpb = import_module(self._tcpb_package)

        with getattr(tcpb, self._tcpb_client)(**self._env_vars) as client:
            return client.compute(input_model)
