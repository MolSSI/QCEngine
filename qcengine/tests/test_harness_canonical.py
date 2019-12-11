"""
Tests the DQM compute dispatch module
"""


import numpy as np
import pytest
from qcelemental.models import AtomicInput

import qcengine as qcng
from qcengine.testing import has_program

_canonical_methods = [
    ("dftd3", {"method": "b3lyp-d3"}),
    ("molpro", {"method": "hf", "basis": "6-31G"}),
    ("mopac", {"method": "PM6"}),
    ("mp2d", {"method": "MP2-DMP2"}),
    ("psi4", {"method": "hf", "basis": "6-31G"}),
    ("qchem", {"method": "hf", "basis": "6-31G"}),
    ("rdkit", {"method": "UFF"}),
    ("torchani", {"method": "ANI1x"}),
    ("entos", {"method": "pbe", "basis": "6-31G"}),
    ("turbomole", {"method": "pbe", "basis": "6-31G"}),
]


@pytest.mark.parametrize("program, model", _canonical_methods)
def test_compute_energy(program, model):
    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    inp = AtomicInput(molecule=qcng.get_molecule("hydrogen"), driver="energy", model=model)
    ret = qcng.compute(inp, program, raise_error=True)

    assert ret.success is True
    assert isinstance(ret.return_result, float)


@pytest.mark.parametrize("program, model", _canonical_methods)
def test_compute_gradient(program, model):
    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    inp = AtomicInput(
        molecule=qcng.get_molecule("hydrogen"), driver="gradient", model=model, extras={"mytag": "something"}
    )
    ret = qcng.compute(inp, program, raise_error=True)

    assert ret.success is True
    assert isinstance(ret.return_result, np.ndarray)
    assert len(ret.return_result.shape) == 2
    assert ret.return_result.shape[1] == 3
    assert "mytag" in ret.extras, ret.extras


@pytest.mark.parametrize(
    "program, model",
    [
        ("dftd3", {"method": "bad"}),
        ("dftd3", {"method": "b3lyp-d3", "driver": "hessian"}),
        ("mopac", {"method": "bad"}),
        ("mp2d", {"method": "bad"}),
        ("psi4", {"method": "bad"}),
        ("qchem", {"method": "bad"}),
        ("rdkit", {"method": "bad"}),
        ("torchani", {"method": "bad"}),
        ("entos", {"method": "bad"}),
        ("turbomole", {"method": "bad"}),
    ],
)
def test_compute_bad_models(program, model):
    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    adriver = model.pop("driver", "energy")
    amodel = model
    inp = AtomicInput(molecule=qcng.get_molecule("hydrogen"), driver=adriver, model=amodel)

    with pytest.raises(qcng.exceptions.InputError) as exc:
        ret = qcng.compute(inp, program, raise_error=True)
