"""
Tests the DQM compute dispatch module
"""

import copy
import pytest

import qcengine as dc
from . import addons

_base_json = {"schema_name": "qc_schema_input", "schema_version": 1}


def test_missing_key():
    ret = dc.compute({"hello": "hi"}, "bleh")
    assert ret["success"] is False
    assert "hello" in ret


@addons.using_psi4
def test_psi4_task():
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
def test_rdkit_task():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = dc.get_molecule("water")
    json_data["driver"] = "gradient"
    json_data["model"] = {"method": "UFF", "basis": ""}
    json_data["keywords"] = {}
    json_data["return_output"] = False

    ret = dc.compute(json_data, "rdkit")

    assert ret["success"] is True


@addons.using_rdkit
def test_rdkit_connectivity_error():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = dc.get_molecule("water")
    json_data["driver"] = "gradient"
    json_data["model"] = {"method": "UFF", "basis": ""}
    json_data["keywords"] = {}
    json_data["return_output"] = False
    del json_data["molecule"]["connectivity"]

    ret = dc.compute(json_data, "rdkit")
    assert ret["success"] is False
    assert "connectivity" in ret["error_message"]

    with pytest.raises(ValueError):
        ret = dc.compute(json_data, "rdkit", raise_error=True)
