"""
A small list of molecules used to validate and tests computation.
"""

import copy

_test_mols = {
    "hydrogen": {
        "symbols": ["H", "H"],
        "geometry": [0, 0, -0.65, 0.0, 0.0, 0.65],
        "molecular_multiplicity": 1,
        "connectivity": [[0, 1, 1]]
    },
    "lithium": {
        "symbols": ["Li"],
        "geometry": [0, 0, 0],
        "molecular_multiplicity": 2,
        "connectivity": []
    },
    "water": {
        "geometry": [
            0.0, 0.0, -0.1294769411935893, 0.0, -1.494187339479985, 1.0274465079245698, 0.0, 1.494187339479985,
            1.0274465079245698
        ],
        "symbols": ["O", "H", "H"],
        "connectivity": [[0, 1, 1], [0, 2, 1]]
    },
}


def get_molecule(name):
    """
    Returns a QC JSON representation of a test molecule.
    """
    if name not in _test_mols:
        raise KeyError("Molecule name '{}' not found".format(name))

    return copy.deepcopy(_test_mols[name])
