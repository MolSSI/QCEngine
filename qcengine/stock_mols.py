"""
A small list of molecules used to validate and tests computation.
"""

import copy

from qcelemental.models import Molecule

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
    "eneyne": {
        "symbols": ["C", "C", "H", "H", "H", "H", "C", "C", "H", "H"],
        "fragments": [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9]],
        "geometry": [0.000000,  -0.667578,  -2.124659,
                     0.000000,   0.667578,  -2.124659,
                     0.923621,  -1.232253,  -2.126185,
                    -0.923621,  -1.232253,  -2.126185,
                    -0.923621,   1.232253,  -2.126185,
                     0.923621,   1.232253,  -2.126185,
                     0.000000,   0.000000,   2.900503,
                     0.000000,   0.000000,   1.693240,
                     0.000000,   0.000000,   0.627352,
                     0.000000,   0.000000,   3.963929],
    },
} # yapf: disable


def get_molecule(name):
    """
    Returns a QC JSON representation of a test molecule.
    """
    if name not in _test_mols:
        raise KeyError("Molecule name '{}' not found".format(name))

    return copy.deepcopy(Molecule(**_test_mols[name]))
