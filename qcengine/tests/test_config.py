"""
Tests the DQM compute module configuration
"""


import pytest
import copy
import os

import qcengine as dc
from . import addons

@pytest.fixture
def opt_state():
    """
    Capture the options state and temporarily override.
    """
    tmp = copy.deepcopy(dc.config._globals)

    base_path = os.path.dirname(os.path.abspath(__file__))
    dc.load_options(os.path.join(base_path, "conf_test.yaml")) 

    yield

    dc.config._globals = tmp

def test_config_path(opt_state):
    cpath = dc.config.get_global("config_path")
    test_path = os.path.join("qcengine", "test")
    assert test_path in cpath

def test_get_default(opt_state):
    assert dc.get_config()["memory_per_job"] == 4
    assert dc.get_config(key="memory_per_job") == 4


def test_hostname_matches(opt_state):
    assert dc.get_config(hostname="dt5")["memory_per_job"] == 60
    assert dc.get_config(key="memory_per_job", hostname="dt5") == 60


def test_default_matches(opt_state):
    """
    Tests that defaults properly come through on different compute
    """
    bench = "/home/user/psi4"
    assert dc.get_config(hostname="nr5")["psi_path"] == bench
    assert dc.get_config(key="psi_path", hostname="nr5") == bench


def test_environmental_vars(opt_state):

    assert dc.get_config("scratch_directory") == "/tmp/"
    assert dc.get_config("scratch_directory", hostname="dt5") is None
