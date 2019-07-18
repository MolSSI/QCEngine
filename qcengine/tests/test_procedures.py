"""
Tests the DQM compute dispatch module
"""

import copy

import numpy as np
import pytest
from qcelemental.models import OptimizationInput

import qcengine as qcng
from qcengine import testing

_base_json = {
    "schema_name": "qcschema_optimization_input",
    "schema_version": 1,
    "keywords": {
        "coordsys": "tric",
        "maxiter": 100,
        "program": None
    },
    "input_specification": {
        "schema_name": "qcschema_input",
        "schema_version": 1,
        "driver": "gradient",
        "model": None,
        "keywords": {},
    },
    "initial_molecule": None
}


@testing.using_psi4
@testing.using_geometric
def test_geometric_psi4():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = qcng.get_molecule("hydrogen")
    inp["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
    inp["input_specification"]["keywords"] = {"scf_properties": ["wiberg_lowdin_indices"]}
    inp["keywords"]["program"] = "psi4"

    inp = OptimizationInput(**inp)

    ret = qcng.compute_procedure(inp, "geometric", raise_error=True)
    assert 10 > len(ret.trajectory) > 1

    assert pytest.approx(ret.final_molecule.measure([0, 1]), 1.e-4) == 1.3459150737
    assert ret.provenance.creator.lower() == "geometric"
    assert ret.trajectory[0].provenance.creator.lower() == "psi4"

    # Check keywords passing
    for single in ret.trajectory:
        assert "scf_properties" in single.keywords
        assert "WIBERG_LOWDIN_INDICES" in single.extras["qcvars"]


@testing.using_psi4
@testing.using_geometric
def test_geometric_local_options():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = qcng.get_molecule("hydrogen")
    inp["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
    inp["keywords"]["program"] = "psi4"

    inp = OptimizationInput(**inp)

    # Set some extremely large number to test
    ret = qcng.compute_procedure(inp, "geometric", raise_error=True, local_options={"memory": "5000"})
    assert pytest.approx(ret.trajectory[0].provenance.memory, 1) == 4900

    # Make sure we cleaned up
    assert "_qcengine_local_config" not in ret.input_specification
    assert "_qcengine_local_config" not in ret.trajectory[0].extras


@testing.using_rdkit
@testing.using_geometric
def test_geometric_stdout():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = qcng.get_molecule("water")
    inp["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    inp["keywords"]["program"] = "rdkit"

    inp = OptimizationInput(**inp)

    ret = qcng.compute_procedure(inp, "geometric", raise_error=True)
    assert ret.success is True
    assert "Converged!" in ret.stdout


@testing.using_geometric
@testing.using_rdkit
def test_geometric_rdkit_error():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = qcng.get_molecule("water").copy(exclude={"connectivity"})
    inp["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    inp["keywords"]["program"] = "rdkit"

    inp = OptimizationInput(**inp)

    ret = qcng.compute_procedure(inp, "geometric")
    assert ret.success is False
    assert isinstance(ret.error.error_message, str)


@testing.using_geometric
@pytest.mark.parametrize("program, model, bench", [
    pytest.param("rdkit", {"method": "UFF"}, [1.87130923886072, 2.959448636243545, 104.50996425790237], marks=testing.using_rdkit),
    pytest.param("torchani", {"method": "ANI1x"}, [1.825818737501941, 2.8663765267932697, 103.4332610730292], marks=testing.using_torchani),
    pytest.param("mopac", {"method": "PM6"}, [1.793054066650722, 2.893333237502448, 107.57221117353501], marks=testing.using_mopac),
])
def test_geometric_generic(program, model, bench):
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = qcng.get_molecule("water")
    inp["input_specification"]["model"] = model
    inp["keywords"]["program"] = program

    ret = qcng.compute_procedure(inp, "geometric", raise_error=True)
    assert ret.success is True
    assert "Converged!" in ret.stdout

    r01, r02, r12, a102 = ret.final_molecule.measure([[0, 1], [0, 2], [1, 2], [1, 0, 2]])

    assert pytest.approx(r01, 1.e-4) == bench[0]
    assert pytest.approx(r02, 1.e-4) == bench[0]
    assert pytest.approx(r12, 1.e-4) == bench[1]
    assert pytest.approx(a102, 1.e-4) == bench[2]
