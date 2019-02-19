"""
Base file for the dqm_compute module.
"""

from . import config

from .compute import compute, compute_procedure
from .config import get_config
from .stock_mols import get_molecule

# Handle versioneer
from .extras import get_information
__version__ = get_information('version')
__git_revision__ = get_information('git_revision')
