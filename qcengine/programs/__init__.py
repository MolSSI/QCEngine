"""
Imports the various compute backends
"""

from typing import List, Set

from .psi4 import Psi4Executor
from .rdkit import RDKitExecutor
from .torchani import TorchANIExecutor
from .molpro import MolproExecutor
from .dftd3 import DFTD3Executor

__all__ = ["register_program", "get_program", "list_all_programs", "list_available_programs"]

programs = {}


def register_program(entry_point: 'ProgramExecutor') -> None:
    """
    Register a new ProgramExecutor with QCEngine
    """

    name = entry_point.name
    if name.lower() in programs.keys():
        raise ValueError('{} is already a registered program.'.format(name))

    programs[name.lower()] = entry_point


def get_program(name: str) -> 'ProgramExecutor':
    """
    Returns a programs executor class
    """
    return programs[name.lower()]


def list_all_programs() -> Set[str]:
    """
    List all programs registered by QCEngine.
    """
    return set(programs.keys())


def list_available_programs() -> Set[str]:
    """
    List all programs that can be exectued (found) by QCEngine.
    """

    ret = set()
    for k, p in programs.items():
        if p.found():
            ret.add(k)

    return ret


register_program(Psi4Executor())
register_program(RDKitExecutor())
register_program(TorchANIExecutor())
register_program(MolproExecutor())
register_program(DFTD3Executor())
