"""
Base file for the dqm_compute module.
"""

from . import dqm_config

from .dqm_config import get_config
from .psi_compute import run_psi4
from .compute import compute

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
