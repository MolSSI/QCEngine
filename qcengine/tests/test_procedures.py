"""
Tests the DQM compute dispatch module
"""

import copy

import pytest

import numpy as np

import qcengine as dc
from qcelemental.models import OptimizationInput
from . import addons

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


def _bond_dist(geom, a1, a2):
    """
    Computes a simple bond distance between two rows in a flat (n, 3) list of coordinates
    """
    if isinstance(geom, np.ndarray):
        geom = geom.flatten().tolist()
    a13 = a1 * 3
    a23 = a2 * 3

    xd = (geom[a13] - geom[a23])**2
    yd = (geom[a13 + 1] - geom[a23 + 1])**2
    zd = (geom[a13 + 2] - geom[a23 + 2])**2

    return (xd + yd + zd)**0.5


@addons.using_psi4
@addons.using_geometric
def test_geometric_psi4():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = dc.get_molecule("hydrogen")
    inp["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
    inp["keywords"]["program"] = "psi4"

    inp = OptimizationInput(**inp)

    ret = dc.compute_procedure(inp, "geometric", raise_error=True)
    assert 10 > len(ret["trajectory"]) > 1

    geom = ret["final_molecule"]["geometry"]
    assert pytest.approx(_bond_dist(geom, 0, 1), 1.e-4) == 1.3459150737


@addons.using_psi4
@addons.using_geometric
def test_geometric_local_options():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = dc.get_molecule("hydrogen")
    inp["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
    inp["keywords"]["program"] = "psi4"

    inp = OptimizationInput(**inp)

    # Set some extremely large number to test
    ret = dc.compute_procedure(inp, "geometric", raise_error=True, local_options={"memory": "5000"})
    assert pytest.approx(ret["trajectory"][0]["provenance"]["memory"], 1) == 4900

    # Make sure we cleaned up
    assert "_qcengine_local_config" not in ret["input_specification"]
    assert "_qcengine_local_config" not in ret["trajectory"][0]


@addons.using_rdkit
@addons.using_geometric
def test_geometric_stdout():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = dc.get_molecule("water")
    inp["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    inp["keywords"]["program"] = "rdkit"

    inp = OptimizationInput(**inp)

    ret = dc.compute_procedure(inp, "geometric", raise_error=True)
    assert ret["success"] is True
    assert "Converged!" in ret["stdout"]
    assert ret["stderr"] == "No stderr recieved."

    with pytest.raises(ValueError):
        _ = dc.compute_procedure(inp, "rdkit", raise_error=True)


@addons.using_rdkit
@addons.using_geometric
def test_geometric_rdkit_error():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = dc.get_molecule("water").copy(exclude="connectivity")
    inp["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    inp["keywords"]["program"] = "rdkit"

    inp = OptimizationInput(**inp)

    ret = dc.compute_procedure(inp, "geometric")
    assert ret["success"] is False
    assert isinstance(ret["error"]["error_message"], str)

    with pytest.raises(ValueError):
        _ = dc.compute_procedure(inp, "rdkit", raise_error=True)


@addons.using_torchani
@addons.using_geometric
def test_geometric_torchani():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = dc.get_molecule("water")
    inp["input_specification"]["model"] = {"method": "ANI1", "basis": None}
    inp["keywords"]["program"] = "torchani"

    ret = dc.compute_procedure(inp, "geometric", raise_error=True)
    assert ret["success"] is True
    assert "Converged!" in ret["stdout"]
    assert ret["stderr"] == "No stderr recieved."
