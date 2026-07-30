"""Compute quantum chemistry using the Valeev group's MPQC4 executable."""

import pprint
from typing import Any, ClassVar, Dict

from qcelemental.util import safe_version, which

from ...config import TaskConfig
from ...exceptions import UnknownError
from ...util import execute
from ..model import ProgramHarness

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class MPQCHarness(ProgramHarness):
    """Interface for the MPQC4 project.

    Notes
    -----
    * Input is a JSON KeyVal document. Results are read from the
      ``Output KeyVal (format=JSON)`` block MPQC writes to stdout.
    * Only ``Energy`` and ``ExcitationEnergy`` properties are supported;
      ``gradient`` and ``hessian`` drivers raise ``InputError``.
    * Basis sets resolve from each directory in ``MPQC_BASIS_PATH`` first
      (colon-separated, in listed order), and only then from libint2's
      bundled library, whose location ``LIBINT_DATA_PATH`` sets. So
      ``MPQC_BASIS_PATH`` shadows ``LIBINT_DATA_PATH`` rather than replacing
      it -- the libint2 library stays the last-resort fallback. Both are
      inherited from the parent environment; see ``_build_environment``.
    """

    _defaults: ClassVar[Dict[str, Any]] = {
        "name": "MPQC",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": True,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "mpqc",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via https://github.com/ValeevGroup/mpqc4",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("mpqc")
        if which_prog not in self.version_cache:
            # `mpqc -v` prints e.g. "MPQC version 4.0.0-beta.1" and exits 0
            # without needing an input file.
            success, output = execute([which_prog, "-v"], {})
            if not success:
                raise UnknownError(output["stderr"])

            version = None
            for line in output["stdout"].splitlines():
                if line.startswith("MPQC version"):
                    version = line.split()[2]
                    break
            if version is None:
                raise UnknownError(f"Could not parse MPQC version from: {output['stdout']}")

            self.version_cache[which_prog] = safe_version(version)

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: TaskConfig) -> "AtomicResult":
        raise NotImplementedError("The MPQC harness cannot yet run calculations.")
