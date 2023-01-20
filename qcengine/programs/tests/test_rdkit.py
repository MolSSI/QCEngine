"""Tests for RDKit functionality"""
import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using
from qcengine import stock_mols

stock = stock_mols._test_mols

@pytest.fixture()
def test_molecule(name):
    mol = stock_mols.get_molecule(name)

    return mol

def test_mol_object(test_molecule("L-tyrosine"))
