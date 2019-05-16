"""
Imports the various compute backends
"""

from typing import List, Set

from .cfour import CfourExecutor
from .dftd3 import DFTD3Executor
from .gamess import GAMESSExecutor
from .nwchem import NWChemExecutor
from .molpro import MolproExecutor
from .mp2d import MP2DExecutor
from .psi4 import Psi4Executor
from .rdkit import RDKitExecutor
from .terachem import TeraChemExecutor
from .torchani import TorchANIExecutor

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
register_program(TeraChemExecutor())
register_program(MP2DExecutor())
register_program(GAMESSExecutor())
register_program(NWChemExecutor())
register_program(CfourExecutor())
