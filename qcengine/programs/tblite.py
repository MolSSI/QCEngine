"""
Harness for the tblite package.
This implementation interfaces with the tblite Python package, which provides
native QCSchema support for the GFN1-xTB, GFN2-xTB, and IPEA1-xTB methods.
"""

from typing import Any, ClassVar, Dict

from qcelemental.models.v2 import AtomicInput, AtomicResult
from qcelemental.util import safe_version, which_import

from ..config import TaskConfig
from ..exceptions import ResourceError, InputError
from ..util import environ_context
from .model import ProgramHarness


class TBLiteHarness(ProgramHarness):
    """Calculation harness for the tblite package."""

    _defaults: ClassVar[Dict[str, Any]] = {
        "name": "tblite",
        "scratch": False,
        "thread_safe": True,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Check for the availability of tblite"""

        return which_import(
            "tblite",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install tblite -c conda-forge`.",
        )

    def get_version(self) -> str:
        """Return the currently used version of tblite"""
        self.found(raise_error=True)

        which_prog = which_import("tblite")
        if which_prog not in self.version_cache:
            from tblite.qcschema import get_version

            self.version_cache[which_prog] = safe_version(".".join(map(str, get_version())))

        return self.version_cache[which_prog]

    def compute(self, input_data: AtomicInput, config: TaskConfig) -> AtomicResult:
        """
        Actual interface to the tblite package. The compute function is just a thin
        wrapper around the native QCSchema interface of tblite.
        """
        from ..testing import is_program_new_enough

        self.found(raise_error=True)

        if not is_program_new_enough("tblite", "0.6.0"):
            raise ResourceError(f"tblite version '{self.get_version()}' too old. Please update to at least '0.6.0'.")

        from tblite.qcschema import run_schema

        # tblite v0.6.0 natively speaks QCSchema v2 and returns AtomicResult or FailedOperation
        env = {
            "OMP_STACKSIZE": f"{config.memory} G",
            "OMP_NUM_THREADS": f"{config.ncores},1",
        }
        with environ_context(env=env):
            output = run_schema(input_data)

        if not output.success:
            raise InputError(output.error.error_message)
 
        return output
