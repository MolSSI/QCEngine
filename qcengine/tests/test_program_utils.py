"""
Tests the DQM compute dispatch module
"""

import pytest

import qcengine as qcng
from qcengine import testing

def test_list_programs():

    r = qcng.list_all_programs()
    assert r >= {"psi4", "rdkit", "molpro", "dftd3"}

@pytest.mark.parametrize("program", [
    pytest.param("psi4", marks=testing.using_psi4),
    pytest.param("torchani", marks=testing.using_torchani),
    pytest.param("rdkit", marks=testing.using_rdkit),
    ])
def test_check_program_avail(program):

    assert program in qcng.list_available_programs()

