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

def test_dispatch():

    data = {"something": "whatever"}
    results = dc.compute(data, "pass")
    
    assert results["success"] == True
    assert results["something"] == "whatever"
    
