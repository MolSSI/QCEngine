"""
Imports the various compute backends
"""

from .psi4 import Psi4Executor
from .rdkit import RDKitExecutor
from .torchani import TorchANIExecutor
from .molpro import MolproExecutor

__all__ = ["register_program", "get_program", "list_programs"]

programs = {}


def register_program(entry_point):
    name = entry_point.name
    if name.lower() in programs.keys():
        raise ValueError('{} is already a registered program.'.format(name))

    programs[name.lower()] = {'entry_point': entry_point}


def get_program(name):
    return programs[name.lower()]['entry_point']


def list_programs():
    return programs.keys()


register_program(Psi4Executor())
register_program(RDKitExecutor())
register_program(TorchANIExecutor())
register_program(MolproExecutor())
