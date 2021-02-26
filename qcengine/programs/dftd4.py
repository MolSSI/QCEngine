"""
Harness for the DFT-D4 dispersion correction.
This implementation interfaces with the dftd4 Python-API, which provides
native support for QCSchema.

Therefore, this harness only has to provide a thin wrapper to integrate dftd4.
"""

from typing import Dict

from qcelemental.models import AtomicInput, AtomicResult
from qcelemental.util import safe_version, which_import

from ..config import TaskConfig
from .model import ProgramHarness


class DFTD4Harness(ProgramHarness):
    """Calculation harness for the DFT-D4 dispersion correction."""

    _defaults = {
        "name": "dftd4",
        "scratch": False,
        "thread_safe": True,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": False,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Check for the availability of the Python API of dftd4"""

        return which_import(
            "dftd4",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via a dftd4 version with enabled Python API",
        )

    def get_version(self) -> str:
        """Return the currently used version of dftd4"""
        self.found(raise_error=True)

        which_prog = which_import("dftd4")
        if which_prog not in self.version_cache:
            import dftd4

            self.version_cache[which_prog] = safe_version(dftd4.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_data: AtomicInput, config: TaskConfig) -> AtomicResult:
        """
        Actual interface to the dftd4 package. The compute function is just a thin
        wrapper around the native QCSchema interface of the dftd4 Python-API.
        """

        self.found(raise_error=True)

        import dftd4
        from dftd4.qcschema import run_qcschema

        # Run the Harness
        output = run_qcschema(input_data)

        # Make sure all keys from the initial input spec are sent along
        output.extras.update(input_data.extras)
        return output
