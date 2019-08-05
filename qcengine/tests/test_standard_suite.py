"""
Tests the DQM compute dispatch module
"""

import copy

import pytest
from qcelemental.models import Molecule, ResultInput

import qcengine as qcng
from qcengine import testing

_canonical_methods = [
    ("dftd3", {"method": "b3lyp-d3"}),
    ("psi4", {"method": "hf", "basis": "6-31G"}),
    ("rdkit", {"method": "UFF"}),
    ("torchani", {"method": "ANI1x"}),
    ("mopac", {"method": "PM6"}),
] # yapf: disable

@pytest.mark.parametrize("program, model", _canonical_methods)
def test_compute_energy(program, model):
    if not testing.has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    inp = ResultInput(molecule=qcng.get_molecule("hydrogen"), driver="energy", model=model)
    ret = qcng.compute(inp, program, raise_error=True)

    assert ret.success is True
    assert isinstance(ret.return_result, float)


@pytest.mark.parametrize("program, model", _canonical_methods)
def test_compute_gradient(program, model):
    if not testing.has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    inp = ResultInput(molecule=qcng.get_molecule("hydrogen"),
                      driver="gradient",
                      model=model,
                      extras={"mytag": "something"})
    ret = qcng.compute(inp, program, raise_error=True)

    assert ret.success is True
    assert isinstance(ret.return_result, list)
    assert "mytag" in ret.extras, ret.extras


@pytest.mark.parametrize("program, model", [
    ("dftd3", {"method": "bad"}),
    ("psi4", {"method": "bad"}),
    ("rdkit", {"method": "bad"}),
    ("torchani", {"method": "bad"}),
    ("dftd3", {"method": "b3lyp-d3", "driver": "hessian"}),
    ("mopac", {"method": "bad"}),
]) # yapf: disable
def test_compute_bad_models(program, model):
    if not testing.has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    adriver = model.pop("driver", "energy")
    amodel = model
    inp = ResultInput(molecule=qcng.get_molecule("hydrogen"), driver=adriver, model=amodel)

    with pytest.raises(qcng.exceptions.InputError) as exc:
        ret = qcng.compute(inp, program, raise_error=True)
