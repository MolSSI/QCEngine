"""
Tests the DQM compute module configuration
"""

# Important this is first
import addons

import os
import qcengine as dc

def test_config_path():
    cpath = dc.config.get_global("config_path")
    test_path = os.path.join("qcengine", "test")
    assert test_path in cpath

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


def test_environmental_vars():

    assert dc.get_config("scratch_directory") == "/tmp/"
    assert dc.get_config("scratch_directory", hostname="dt5") is None
