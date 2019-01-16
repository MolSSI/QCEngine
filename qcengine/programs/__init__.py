"""
Imports the various compute backends
"""

programs = {}

def register_program(name, entry_point):
    if name in programs.keys():
        raise ValueError('{} is already a registered program.'.format(name))

    programs[name] = { 'entry_point': entry_point }


def get_program(name):
    return programs[name]['entry_point']


def get_programs():
    return programs.keys()


from . import psi4
register_program('psi4', psi4.psi4)

from . import rdkit
register_program('rdkit', rdkit.rdkit)

from . import torchani
register_program('torchani', torchani.torchani)
