"""
Tests the DQM compute dispatch module
"""

# Important this is first
import addons

import os
import dqm_compute as dc

@addons.using_psi4
def test_psi4():
    json_data = {}
    json_data["molecule"] = """He 0 0 0\n--\nHe 0 0 1"""
    json_data["driver"] = "energy"
    json_data["method"] = "SCF"
    json_data["options"] = {"BASIS": "STO-3G"}
    json_data["return_output"] = False

    ret = dc.compute(json_data, "psi4")

    assert ret["driver"] == "energy"
    for key in ["provenance", "wall_time"]:
        assert key in ret

    for key in ["cpu", "hostname", "username"]:
        assert key in ret["provenance"]

    assert ret["success"] == True

@addons.using_psi4
def test_psi4_switch():
    json_data = {}
    json_data["molecule"] = """Li"""
    json_data["driver"] = "energy"
    json_data["method"] = "SCF"
    json_data["options"] = {"BASIS": "STO-3G"}
    json_data["return_output"] = False

    ret = dc.compute(json_data, "psi4")

    assert ret["success"] == True

def test_dispatch():

    data = {"something": "whatever"}
    results = dc.compute(data, "pass")

    assert results["success"] == True
    assert results["something"] == "whatever"


@addons.using_rdkit
def test_rdkit():
    json_data = {}
    json_data["molecule"] = addons.test_mols["water"]
    json_data["molecule"]["connectivity"] = [[0, 1, 1], [0, 2, 1]]
    json_data["driver"] = "gradient"
    json_data["model"] = {"method": "UFF", "basis": ""}
    json_data["options"] = {}
    json_data["return_output"] = False

    ret = dc.compute(json_data, "rdkit")

    assert ret["success"] == True
