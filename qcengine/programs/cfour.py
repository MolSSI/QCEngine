"""Compute quantum chemistry using Mainz-Austin-Budapest-Gainesville's CFOUR executable."""

import pprint
from typing import Dict

from qcelemental.util import which, safe_version

from ..util import execute
from .executor import ProgramExecutor

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class CFOURExecutor(ProgramExecutor):

    _defaults = {
        "name": "CFOUR",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] ={}

    class Config(ProgramExecutor.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @staticmethod
    def found(raise_error: bool=False) -> bool:
        return which('xcfour', return_bool=True, raise_error=raise_error, raise_msg='Please install via http://cfour.de/')

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which('xcfour')
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "ZMAT"], {"ZMAT": "\nHe\n\n"})

            if success:
                for line in output["stdout"].splitlines():
                    if 'Version' in line:
                        branch = ' '.join(line.strip().split()[1:])
                self.version_cache[which_prog] = safe_version(branch)

        return self.version_cache[which_prog]

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        self.found(raise_error=True)
