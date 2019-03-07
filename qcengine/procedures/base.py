"""
Imports the various procedure backends
"""

from typing import List, Set

from .geometric import GeometricProcedure

__all__ = ["register_procedure", "get_procedure", "list_all_procedures", "list_available_procedures"]

procedures = {}


def register_procedure(entry_point: 'BaseProcedure') -> None:
    """
    Register a new BaseProcedure with QCEngine
    """

    name = entry_point.name
    if name.lower() in procedures.keys():
        raise ValueError('{} is already a registered procedure.'.format(name))

    procedures[name.lower()] = entry_point


def get_procedure(name: str) -> 'BaseProcedure':
    """
    Returns a procedures executor class
    """
    return procedures[name.lower()]


def list_all_procedures() -> Set[str]:
    """
    List all procedures registered by QCEngine.
    """
    return set(procedures.keys())


def list_available_procedures() -> Set[str]:
    """
    List all procedures that can be exectued (found) by QCEngine.
    """

    ret = set()
    for k, p in procedures.items():
        if p.found():
            ret.add(k)

    return ret


register_procedure(GeometricProcedure())
