"""
Tests the DQM compute dispatch module
"""

# Important this is first
import addons

import os
import dqm_compute as dc

@addons.using_psi4
def test_psi4():
    os.environ["PSI_SCRATCH"] = "/tmp"
    ret = dc.test_psi4()

    assert ret["driver"] == "energy"
    for key in ["provenance", "wall_time"]:
        assert key in ret

    for key in ["cpu", "hostname", "username"]:
        assert key in ret["provenance"]
    print(ret)

@addons.using_psi4
def test_psi4_switch():
    json_data = {}
    json_data["molecule"] = """Li"""
    json_data["driver"] = "energy"
    json_data["method"] = 'SCF'
    json_data["options"] = {"BASIS": "STO-3G"}
    json_data["return_output"] = False
    os.environ["PSI_SCRATCH"] = "/tmp"

    ret = dc.compute(json_data, "psi4")

    assert ret["success"] == True

def test_dispatch():

    data = {"something": "whatever"}
    results = dc.compute(data, "pass")

    assert results["success"] == True
    assert results["something"] == "whatever"

