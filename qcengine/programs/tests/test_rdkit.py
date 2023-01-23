"""Tests for RDKit functionality"""
import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using
from qcengine import stock_mols

@pytest.fixture
def ala():
    mol = stock_mols.get_molecule("L-alanine")
    return mol

@pytest.fixture
def gly():
    mol = stock_mols.get_molecule("glycine")
    return mol

@pytest.fixture
def his():
    mol = stock_mols.get_molecule("L-histidine")
    return mol

@pytest.fixture
def met():
    mol = stock_mols.get_molecule("L-methionine")
    return mol

@pytest.fixture
def phe():
    mol = stock_mols.get_molecule("L-phenylalanine")
    return mol

@pytest.fixture
def ser():
    mol = stock_mols.get_molecule("L-serine")
    return mol

@pytest.fixture
def trp():
    mol = stock_mols.get_molecule("L-tryptophan")
    return mol

@pytest.fixture
def tyr():
    mol = stock_mols.get_molecule("L-tyrosine")
    return mol

@pytest.mark.parametrize("molecule, symbols, connectivity", [
    ("L-tyrosine", [
        "C", "C", "C", "C", "C", "C", "O", "C", "C", "N", "C", "O", "O"
    ], [
        (1, 2, 2.0),
        (1, 6, 1.0),
        (2, 3, 1.0),
        (2, 8, 1.0),
        (3, 4, 2.0),
        (4, 5, 1.0),
        (5, 6, 2.0),
        (5, 7, 1.0),
        (8, 9, 1.0),
        (9, 10, 1.0),
        (9, 11, 1.0),
        (11, 12, 2.0),
        (11, 13, 1.0),
    ]),
    ("glycine", [
        "N", "C", "C", "O", "O"
    ], [
        (1, 2, 1.0),
        (2, 3, 1.0),
        (3, 4, 1.0),
        (3, 5, 2.0),
    ]),
    ("L-methionine", [
        "C", "S", "C", "C", "C", "C", "N", "O", "O"
    ], [
        (1, 2, 1.0),
        (2, 3, 1.0),
        (3, 4, 1.0),
        (4, 5, 1.0),
        (5, 6, 1.0),
        (5, 7, 1.0),
        (6, 8, 1.0),
        (6, 9, 2.0),
    ])
])
def test_rdkit_mol_object(molecule, symbols, connectivity):
    mol = stock_mols.get_molecule(molecule)
    sym = mol.symbols
    bonds = mol.connectivity
    assert np.array_equal(sym, symbols)
    assert np.array_equal(bonds, connectivity)

@pytest.mark.parametrize("molecule", [ala, gly, his, met, phe, ser, trp, tyr])
@using("rdkit")
def test_get_descriptors(molecule):
    resi = {
        "molecule": molecule,
        "driver": "properties",
        "model": {"method": "descriptors"},
    }

    result = qcng.compute(resi, "rdkit")
    print(result.return_result["canonical_smiles"])
