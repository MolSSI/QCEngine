"""
Imports the various compute backends
"""

from . import psi4
from . import rdkit
from . import torchani

__all__ = ["register_program", "get_program", "get_programs"]

programs = {}


def register_program(name, entry_point):
    if name in programs.keys():
        raise ValueError('{} is already a registered program.'.format(name))

    programs[name] = {'entry_point': entry_point}


def get_program(name):
    return programs[name]['entry_point']


def get_programs():
    return programs.keys()


register_program('psi4', psi4.psi4)
register_program('rdkit', rdkit.rdkit)
register_program('torchani', torchani.torchani)
