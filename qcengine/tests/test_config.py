"""
Tests the DQM compute module configuration
"""

import copy
import os

import pytest

import qcengine as dc


@pytest.fixture
def opt_state_basic():
    """
    Capture the options state and temporarily override.
    """
    tmp = copy.deepcopy(dc.config._globals)

    base_path = os.path.dirname(os.path.abspath(__file__))
    dc.load_options(os.path.join(base_path, "conf_basic.yaml"))

    yield

    dc.config._globals = tmp


def test_config_path(opt_state_basic):
    cpath = dc.config.get_global("config_path")
    test_path = os.path.join("qcengine", "test")
    assert test_path in cpath


def test_get_default(opt_state_basic):
    assert dc.get_config()["memory_per_job"] == 4
    assert dc.get_config(key="memory_per_job") == 4


def test_hostname_matches(opt_state_basic):
    assert dc.get_config(hostname="dt5")["memory_per_job"] == 60
    assert dc.get_config(key="memory_per_job", hostname="dt5") == 60


def test_default_matches(opt_state_basic):
    """
    Tests that defaults properly come through on different compute
    """
    bench = "/home/user/psi4"
    assert dc.get_config(hostname="nr5")["psi_path"] == bench
    assert dc.get_config(key="psi_path", hostname="nr5") == bench


def test_environmental_vars(opt_state_basic):

    assert dc.get_config("scratch_directory") == "/tmp/"
    assert dc.get_config("scratch_directory", hostname="dt5") is None


@pytest.fixture
def opt_state_auto():
    """
    Capture the options state and temporarily override.
    """
    tmp = copy.deepcopy(dc.config._globals)

    base_path = os.path.dirname(os.path.abspath(__file__))
    dc.load_options(os.path.join(base_path, "conf_auto.yaml"))

    yield

    dc.config._globals = tmp


def test_auto_threads(opt_state_auto):

    assert dc.get_config("jobs_per_node") == 1
    assert isinstance(dc.get_config("nthreads_per_job"), int)
    assert dc.get_config("nthreads_per_job") > 0
    assert dc.get_config("nthreads_per_job") < 100

    assert isinstance(dc.get_config("memory_per_job"), (int, float))
    assert dc.get_config("memory_per_job") > 0.01  # Always more than 1OMB free?


def test_global_repr(opt_state_auto):

    dc.config.global_repr()