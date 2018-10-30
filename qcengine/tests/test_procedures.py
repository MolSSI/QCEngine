"""
Tests the DQM compute dispatch module
"""

import copy
import pytest

import qcengine as dc
from . import addons

_base_json = {
    "schema_name": "qc_schema_optimization_input",
    "schema_version": 1,
    "keywords": {
        "coordsys": "tric",
        "maxiter": 100,
        "program": None
    },
    "input_specification": {
        "schema_name": "qc_schema_input",
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

    ret = dc.compute_procedure(inp, "geometric")
    assert 10 > len(ret["trajectory"]) > 1

    geom = ret["final_molecule"]["geometry"]
    assert pytest.approx(_bond_dist(geom, 0, 1), 1.e-4) == 1.3459150737

@addons.using_rdkit
@addons.using_geometric
def test_geometric_stdout():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = dc.get_molecule("water")
    inp["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    inp["keywords"]["program"] = "rdkit"

    ret = dc.compute_procedure(inp, "geometric")
    assert ret["success"] is True
    assert "Converged!" in ret["stdout"]
    assert ret["stderr"] == "No stderr recieved."

    with pytest.raises(ValueError):
        ret = dc.compute_procedure(inp, "rdkit", raise_error=True)

@addons.using_rdkit
@addons.using_geometric
def test_geometric_rdkit_error():
    inp = copy.deepcopy(_base_json)

    inp["initial_molecule"] = dc.get_molecule("water")
    del inp["initial_molecule"]["connectivity"]
    inp["input_specification"]["model"] = {"method": "UFF", "basis": ""}
    inp["keywords"]["program"] = "rdkit"

    ret = dc.compute_procedure(inp, "geometric")
    assert ret["success"] is False
    assert isinstance(ret["error_message"], str)

    with pytest.raises(ValueError):
        ret = dc.compute_procedure(inp, "rdkit", raise_error=True)
