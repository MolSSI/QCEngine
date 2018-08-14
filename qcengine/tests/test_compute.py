"""
Tests the DQM compute dispatch module
"""

import copy
import pytest

import qcengine as dc
from . import addons

_base_json = {
    "schema_name": "qc_schema_input",
    "schema_version": 1
}

@addons.using_psi4
def test_psi4():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = dc.get_molecule("water")
    json_data["driver"] = "energy"
    json_data["model"] = {"method": "SCF", "basis": "sto-3g"}
    json_data["keywords"] = {"scf_type": "df"}
    json_data["return_output"] = False

    ret = dc.compute(json_data, "psi4")

    assert ret["driver"] == "energy"
    assert "provenance" in ret

    for key in ["cpu", "hostname", "username", "wall_time"]:
        assert key in ret["provenance"]

    assert ret["success"] is True

@addons.using_psi4
def test_psi4_ref_switch():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = dc.get_molecule("lithium")
    json_data["driver"] = "energy"
    json_data["model"] = {"method": "SCF", "basis": "sto-3g"}
    json_data["keywords"] = {"scf_type": "df"}
    json_data["return_output"] = False

    ret = dc.compute(json_data, "psi4")

    assert ret["success"] is True

@addons.using_rdkit
def test_rdkit():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = dc.get_molecule("water")
    json_data["molecule"]["connectivity"] = [[0, 1, 1], [0, 2, 1]]
    json_data["driver"] = "gradient"
    json_data["model"] = {"method": "UFF", "basis": ""}
    json_data["keywords"] = {}
    json_data["return_output"] = False

    ret = dc.compute(json_data, "rdkit")

    assert ret["success"] is True

@addons.using_psi4
@addons.using_geometric
def test_geometric():
    qc_schema_input = {
        "schema_name": "qc_schema_input",
        "schema_version": 1,
        "driver": "gradient",
        "model": {
            "method": "HF",
            "basis": "sto-3g"
        },
        "keywords": {},
    }

    json_data = {
        "schema_name": "qc_schema_optimization_input",
        "schema_version": 1,
        "keywords": {
            "coordsys": "tric",
            "maxiter": 100,
            "program": "psi4"
        },
        "input_specification": qc_schema_input,
        "initial_molecule": {
            "geometry": [
                0.0,  0.0, -0.6,
                0.0,  0.0,  0.6,
            ],
            "symbols": ["H", "H"],
            "connectivity": [[0, 1, 1]]
        },
    }

    ret = dc.compute_procedure(json_data, "geometric")
    assert 10 > len(ret["trajectory"]) > 1

    geom = ret["final_molecule"]["geometry"]
    assert pytest.approx(abs(geom[2] - geom[5]), 3) == -1.3459150737
