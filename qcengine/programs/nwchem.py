"""
Calls the NWChem executable.
"""
from typing import Dict

from qcelemental.util import safe_version, which

from .model import ProgramHarness
from ..util import execute


class NWChemHarness(ProgramHarness):

    _defaults = {
        "name": "NWChem",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,  # T but common?
        "node_parallel": True,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] ={}

    class Config(ProgramHarness.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @staticmethod
    def found(raise_error: bool=False) -> bool:
        return which('nwchem', return_bool=True, raise_error=raise_error, raise_msg='Please install via http://www.nwchem-sw.org/index.php/Download')

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which('nwchem')
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "v.nw"], {"v.nw": ""})

            if success:
                for line in output["stdout"].splitlines():
                    if 'nwchem branch' in line:
                        branch = line.strip().split()[-1]
                    if 'nwchem revision' in line:
                        revision = line.strip().split()[-1]
                self.version_cache[which_prog] = safe_version(branch + '+' + revision)

        return self.version_cache[which_prog]

    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs NWChem in executable mode
        """
        self.found(raise_error=True)

