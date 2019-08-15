"""
Tests the DQM compute dispatch module
"""

import copy

import numpy as np
import pytest

import qcengine as qcng
from qcelemental.models import OptimizationInput
from qcengine import testing
from qcengine.testing import failure_engine

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
@testing.using_rdkit
def test_geometric_retries(failure_engine):
    inp = copy.deepcopy(_base_json)

    failure_engine.iter_modes = [
        "random_error", "pass", # Iter 1
        "random_error", "random_error", "pass", # Iter 2
    ] # yapf: disable
    failure_engine.iter_modes.extend(["pass"] * 20)

    inp["initial_molecule"] = {"symbols": ["He", "He"], "geometry": [0, 0, 0, 0, 0, failure_engine.start_distance]}
    inp["input_specification"]["model"] = {"method": "something"}
    inp["keywords"]["program"] = failure_engine.name

    inp = OptimizationInput(**inp)

    ret = qcng.compute_procedure(inp, "geometric", local_options={"ncores": 13}, raise_error=True)
    assert ret.success is True
    assert ret.trajectory[0].provenance.retries == 1
    assert ret.trajectory[0].provenance.ncores == 13
    assert ret.trajectory[1].provenance.retries == 2
    assert ret.trajectory[1].provenance.ncores == 13
    assert "retries" not in ret.trajectory[2].provenance.dict()

    # Ensure we still fail
    failure_engine.iter_modes = [
        "random_error", "pass", # Iter 1
        "random_error", "random_error", "pass", # Iter 2
    ] # yapf: disable
    ret = qcng.compute_procedure(inp, "geometric", local_options={"ncores": 13, "retries": 1})
    assert ret.success is False
    assert ret.input_data["trajectory"][0]["provenance"]["retries"] == 1
    assert len(ret.input_data["trajectory"]) == 2


@testing.using_geometric
@pytest.mark.parametrize("program, model, bench", [
    pytest.param("rdkit", {"method": "UFF"},
                 [1.87130923886072, 2.959448636243545, 104.5099642579023],
                 marks=testing.using_rdkit),
    pytest.param("torchani", {"method": "ANI1x"},
                 [1.82581873750194, 2.866376526793269, 103.4332610730292],
                 marks=testing.using_torchani),
    pytest.param("mopac", {"method": "PM6"},
                 [1.79305406665072, 2.893333237502448, 107.5722111735350],
                 marks=testing.using_mopac),
]) # yapf: disable
def test_geometric_generic(program, model, bench):
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = qcng.get_molecule("water")
    inp["input_specification"]["model"] = model
    inp["keywords"]["program"] = program
    inp["input_specification"]["extras"] = {"_secret_tags": {"mysecret_tag": "data1"}}

    ret = qcng.compute_procedure(inp, "geometric", raise_error=True)
    assert ret.success is True
    assert "Converged!" in ret.stdout

    r01, r02, r12, a102 = ret.final_molecule.measure([[0, 1], [0, 2], [1, 2], [1, 0, 2]])

    assert pytest.approx(r01, 1.e-4) == bench[0]
    assert pytest.approx(r02, 1.e-4) == bench[0]
    assert pytest.approx(r12, 1.e-4) == bench[1]
    assert pytest.approx(a102, 1.e-4) == bench[2]

    assert "_secret_tags" in ret.trajectory[0].extras
    assert "data1" == ret.trajectory[0].extras["_secret_tags"]["mysecret_tag"]