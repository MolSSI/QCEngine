"""
Tests the DQM compute module configuration
"""

import os
import dqm_compute as dc

os.environ["FW_CONFIG_FILE"] = os.path.dirname(os.path.abspath(__file__))


def test_get_default():
    assert dc.get_config()["memory_per_job"] == 4
    assert dc.get_config(key="memory_per_job") == 4


def test_hostname_matches():
    assert dc.get_config(hostname="dt5")["memory_per_job"] == 60
    assert dc.get_config(key="memory_per_job", hostname="dt5") == 60


def test_default_matches():
    """
    Tests that defaults properly come through on different compute
    """
    bench = "/home/user/psi4"
    assert dc.get_config(hostname="nr5")["psi_path"] == bench
    assert dc.get_config(key="psi_path", hostname="nr5") == bench
