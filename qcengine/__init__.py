"""
Base file for the dqm_compute module.
"""
from importlib.metadata import version

__version__ = version("qcengine")
__git_revision__ = "(no longer provided)"
del version

# isort: off
from . import config, exceptions
from .compute import compute, compute_procedure
from .config import get_config
from .extras import get_information
from .mdi_server import MDIServer
from .procedures import get_procedure, list_all_procedures, list_available_procedures
from .programs import get_program, list_all_programs, list_available_programs, register_program, unregister_program
from .stock_mols import get_molecule

# isort: on
