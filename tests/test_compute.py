"""
Tests the DQM compute dispatch module
"""

# Important this is first
import addons

import os
import dqm_compute as dc

def test_dispatch():

    data = {"something": "whatever"}
    results = dc.compute(data, "pass")
    
    assert results["success"] == True
    assert results["something"] == "whatever"
    
