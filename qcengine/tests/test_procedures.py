"""
Tests the DQM compute dispatch module
"""

import pytest

import qcengine as qcng
from qcelemental.models import OptimizationInput
from qcengine import testing
from qcengine.testing import failure_engine


@pytest.fixture(scope="function")
def input_data():
    return {
        "keywords": {"coordsys": "tric", "maxiter": 100, "program": None},
        "input_specification": {"driver": "gradient", "model": None, "keywords": {}},
        "initial_molecule": None,
    }


@testing.using_psi4
@testing.using_geometric
def test_geometric_psi4(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("hydrogen")
    input_data["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
    input_data["input_specification"]["keywords"] = {"scf_properties": ["wiberg_lowdin_indices"]}
    input_data["keywords"]["program"] = "psi4"

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "geometric", raise_error=True)
    assert 10 > len(ret.trajectory) > 1

    assert pytest.approx(ret.final_molecule.measure([0, 1]), 1.0e-4) == 1.3459150737
    assert ret.provenance.creator.lower() == "geometric"
    assert ret.trajectory[0].provenance.creator.lower() == "psi4"

    # Check keywords passing
    for single in ret.trajectory:
        assert "scf_properties" in single.keywords
        assert "WIBERG_LOWDIN_INDICES" in single.extras["qcvars"]


@testing.using_psi4
@testing.using_geometric
def test_geometric_local_options(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("hydrogen")
    input_data["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
    input_data["keywords"]["program"] = "psi4"

    input_data = OptimizationInput(**input_data)

    # Set some extremely large number to test
    ret = qcng.compute_procedure(input_data, "geometric", raise_error=True, local_options={"memory": "5000"})
    assert pytest.approx(ret.trajectory[0].provenance.memory, 1) == 4900

    # Make sure we cleaned up
    assert "_qcengine_local_config" not in ret.input_specification
    assert "_qcengine_local_config" not in ret.trajectory[0].extras


@testing.using_rdkit
@testing.using_geometric
def test_geometric_stdout(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("water")
    input_data["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    input_data["keywords"]["program"] = "rdkit"

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "geometric", raise_error=True)
    assert ret.success is True
    assert "Converged!" in ret.stdout


@testing.using_geometric
@testing.using_rdkit
def test_geometric_rdkit_error(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("water").copy(exclude={"connectivity"})
    input_data["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    input_data["keywords"]["program"] = "rdkit"

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "geometric")
    assert ret.success is False
    assert isinstance(ret.error.error_message, str)


@testing.using_rdkit
@testing.using_geometric
def test_optimization_protocols(input_data):

    input_data["initial_molecule"] = qcng.get_molecule("water")
    input_data["input_specification"]["model"] = {"method": "UFF"}
    input_data["keywords"]["program"] = "rdkit"
    input_data["protocols"] = {"trajectory": "initial_and_final"}

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "geometric", raise_error=True)
    assert ret.success, ret.error.error_message

    assert len(ret.trajectory) == 2
    assert ret.initial_molecule.get_hash() == ret.trajectory[0].molecule.get_hash()
    assert ret.final_molecule.get_hash() == ret.trajectory[1].molecule.get_hash()


@testing.using_geometric
@testing.using_rdkit
def test_geometric_retries(failure_engine, input_data):

    failure_engine.iter_modes = ["random_error", "pass", "random_error", "random_error", "pass"]  # Iter 1  # Iter 2
    failure_engine.iter_modes.extend(["pass"] * 20)

    input_data["initial_molecule"] = {
        "symbols": ["He", "He"],
        "geometry": [0, 0, 0, 0, 0, failure_engine.start_distance],
    }
    input_data["input_specification"]["model"] = {"method": "something"}
    input_data["keywords"]["program"] = failure_engine.name

    input_data = OptimizationInput(**input_data)

    ret = qcng.compute_procedure(input_data, "geometric", local_options={"ncores": 13}, raise_error=True)
    assert ret.success is True
    assert ret.trajectory[0].provenance.retries == 1
    assert ret.trajectory[0].provenance.ncores == 13
    assert ret.trajectory[1].provenance.retries == 2
    assert ret.trajectory[1].provenance.ncores == 13
    assert "retries" not in ret.trajectory[2].provenance.dict()

    # Ensure we still fail
    failure_engine.iter_modes = ["random_error", "pass", "random_error", "random_error", "pass"]  # Iter 1  # Iter 2
    ret = qcng.compute_procedure(input_data, "geometric", local_options={"ncores": 13, "retries": 1})
    assert ret.success is False
    assert ret.input_data["trajectory"][0]["provenance"]["retries"] == 1
    assert len(ret.input_data["trajectory"]) == 2


@testing.using_geometric
@pytest.mark.parametrize(
    "program, model, bench",
    [
        pytest.param(
            "rdkit",
            {"method": "UFF"},
            [1.87130923886072, 2.959448636243545, 104.5099642579023],
            marks=testing.using_rdkit,
        ),
        pytest.param(
            "torchani",
            {"method": "ANI1x"},
            [1.82581873750194, 2.866376526793269, 103.4332610730292],
            marks=testing.using_torchani,
        ),
        pytest.param(
            "mopac",
            {"method": "PM6"},
            [1.7927843431811934, 2.893333237502448, 107.60441967992045],
            marks=testing.using_mopac,
        ),
    ],
)
def test_geometric_generic(input_data, program, model, bench):

    input_data["initial_molecule"] = qcng.get_molecule("water")
    input_data["input_specification"]["model"] = model
    input_data["keywords"]["program"] = program
    input_data["input_specification"]["extras"] = {"_secret_tags": {"mysecret_tag": "data1"}}

    ret = qcng.compute_procedure(input_data, "geometric", raise_error=True)
    assert ret.success is True
    assert "Converged!" in ret.stdout

    r01, r02, r12, a102 = ret.final_molecule.measure([[0, 1], [0, 2], [1, 2], [1, 0, 2]])

    assert pytest.approx(r01, 1.0e-4) == bench[0]
    assert pytest.approx(r02, 1.0e-4) == bench[0]
    assert pytest.approx(r12, 1.0e-4) == bench[1]
    assert pytest.approx(a102, 1.0e-4) == bench[2]

    assert "_secret_tags" in ret.trajectory[0].extras
    assert "data1" == ret.trajectory[0].extras["_secret_tags"]["mysecret_tag"]
