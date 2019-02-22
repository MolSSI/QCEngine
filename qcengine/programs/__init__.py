"""
Imports the various compute backends
"""

from . import psi4
from .rdkit import RDKitExecutor
from . import torchani

__all__ = ["register_program", "get_program", "get_programs"]

programs = {}


def register_program(entry_point):
    name = entry_point.name
    if name.lower() in programs.keys():
        raise ValueError('{} is already a registered program.'.format(name))

    programs[name.lower()] = {'entry_point': entry_point}


def get_program(name):
    return programs[name.lower()]['entry_point']


def get_programs():
    return programs.keys()


# register_program('psi4', psi4.psi4)
register_program(RDKitExecutor())
# register_program('torchani', torchani.torchani)
