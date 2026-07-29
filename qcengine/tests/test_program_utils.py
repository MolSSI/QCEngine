"""
Tests the DQM compute dispatch module
"""

import pytest

import qcengine as qcng
from qcengine.testing import using


def test_list_programs():

    r = qcng.list_all_programs()
    assert r >= {"psi4", "rdkit", "molpro", "dftd3"}


@pytest.mark.parametrize(
    "program",
    [
        pytest.param("psi4", marks=using("psi4")),
        pytest.param("torchani", marks=using("torchani")),
        pytest.param("rdkit", marks=using("rdkit")),
    ],
)
def test_check_program_avail(program):

    assert program in qcng.list_available_programs()


def test_program_version_check_is_targeted(monkeypatch):
    from qcengine.testing import is_program_new_enough

    class AvailableHarness:
        @staticmethod
        def found():
            return True

        @staticmethod
        def get_version():
            return "2.0"

    def unexpected_full_scan():
        pytest.fail("is_program_new_enough should not scan every registered program")

    monkeypatch.setattr(qcng, "list_all_procedures", lambda: set())
    monkeypatch.setattr(qcng, "list_available_programs", unexpected_full_scan)
    monkeypatch.setattr(qcng, "get_program", lambda name, check=False: AvailableHarness())

    assert is_program_new_enough("available", "1.0")


def test_program_avail_bounce():

    with pytest.raises(qcng.exceptions.InputError) as exc:
        qcng.compute({}, "bad_program", raise_error=True)

    assert "not registered" in str(exc.value)


def test_list_procedures():

    r = qcng.list_all_procedures()
    assert r >= {"geometric", "berny"}


@pytest.mark.parametrize(
    "procedure", [pytest.param("geometric", marks=using("geometric")), pytest.param("berny", marks=using("berny"))]
)
def test_check_procedure_avail(procedure):

    assert procedure in qcng.list_available_procedures()


def test_procedure_avail_bounce():

    with pytest.raises(qcng.exceptions.InputError) as exc:
        qcng.compute({}, "bad_program", raise_error=True)

    assert "not registered" in str(exc.value)
