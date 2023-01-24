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


# @pytest.fixture
# def gly():
#     mol = stock_mols.get_molecule("glycine")
#     return mol
#
# @pytest.fixture
# def his():
#     mol = stock_mols.get_molecule("L-histidine")
#     return mol
#
# @pytest.fixture
# def met():
#     mol = stock_mols.get_molecule("L-methionine")
#     return mol
#
# @pytest.fixture
# def phe():
#     mol = stock_mols.get_molecule("L-phenylalanine")
#     return mol
#
# @pytest.fixture
# def ser():
#     mol = stock_mols.get_molecule("L-serine")
#     return mol
#
# @pytest.fixture
# def trp():
#     mol = stock_mols.get_molecule("L-tryptophan")
#     return mol
#
# @pytest.fixture
# def tyr():
#     mol = stock_mols.get_molecule("L-tyrosine")
#     return mol

# @pytest.mark.parametrize("molecule", [ala, gly, his, met, phe, ser, trp, tyr])
# Need to fix/add other amino acids to stock_mols
@using("rdkit")
def test_get_descriptors(ala):
    resi = {
        "molecule": ala,
        "driver": "properties",
        "model": {"method": "descriptors"},
    }

    result = qcng.compute(resi, "rdkit")
    expected = "C[C@H](N)C(=O)O"
    assert result.return_result["canonical_smiles"] == expected
