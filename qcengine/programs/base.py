"""
Imports the various compute backends
"""

from typing import Set
from ..exceptions import InputError, ResourceError

from .cfour import CFOURHarness
from .dftd3 import DFTD3Harness
from .entos import EntosHarness
from .gamess import GAMESSHarness
from .nwchem import NWChemHarness
from .molpro import MolproHarness
from .mopac import MopacHarness
from .mp2d import MP2DHarness
from .psi4 import Psi4Harness
from .rdkit import RDKitHarness
from .terachem import TeraChemHarness
from .torchani import TorchANIHarness

__all__ = ["register_program", "get_program", "list_all_programs", "list_available_programs"]

programs = {}


def register_program(entry_point: 'ProgramHarness') -> None:
    """
    Register a new ProgramHarness with QCEngine.
    """

    name = entry_point.name
    if name.lower() in programs.keys():
        raise ValueError('{} is already a registered program.'.format(name))

    programs[name.lower()] = entry_point


def unregister_program(name: str) -> None:
    """
    Unregisters a given program.
    """

    ret = programs.pop(name.lower(), None)
    if ret is None:
        raise KeyError(f"Program {name} is not registered with QCEngine")


def get_program(name: str, check: bool = True) -> 'ProgramHarness':
    """
    Returns a program's executor class

    Parameters
    ----------
    check
        ``True`` Do raise error if program not found. ``False`` is handy for
        the specialized case of calling non-execution methods (like parsing for testing)
        on the returned ``Harness``.

    """
    name = name.lower()

    if name not in programs:
        raise InputError(f"Program {name} is not registered to QCEngine.")

    ret = programs[name]
    if check and not ret.found():
        raise ResourceError(f"Program {name} is registered with QCEngine, but cannot be found.")

    return ret


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


register_program(Psi4Harness())
register_program(RDKitHarness())
register_program(TorchANIHarness())
register_program(MolproHarness())
register_program(MopacHarness())
register_program(DFTD3Harness())
register_program(TeraChemHarness())
register_program(MP2DHarness())
#register_program(GAMESSHarness())
#register_program(NWChemHarness())
#register_program(CFOURHarness())
register_program(EntosHarness())
