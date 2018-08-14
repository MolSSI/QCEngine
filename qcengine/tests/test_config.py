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

    # Snapshot env
    old_globals = copy.deepcopy(dc.config._globals)
    old_environ = dict(os.environ)

    os.environ["TMPDIR"] = "/tmp/"

    config = {
      "default_compute": {
        "psi_path": "/home/user/psi4",
        "jobs_per_node": 1,
        "nthreads_per_job": 2,
        "memory_per_job": 4,
        "scratch_directory": "$TMPDIR"
      },
      "other_compute": {
        "dragonsooth": {
          "psi_path": "/home/user/dt/psi4",
          "hostname": "dt*",
          "jobs_per_node": 2,
          "nthreads_per_job": 6,
          "memory_per_job": 60,
          "scratch_directory": "$NOVAR_RANDOM_ABC123"
        },
        "new_river": {
          "hostname": "nr*",
          "jobs_per_node": 2,
          "nthreads_per_job": 12,
          "memory_per_job": 120
        }
      }
    }

    dc.load_options(config)

    yield

    # Reset env
    os.environ.update(old_environ)
    dc.config._globals = old_globals


def test_config_path(opt_state_basic):
    cpath = dc.config.get_global("config_path")
    assert cpath is None


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

    print(dc.get_config("scratch_directory"))
    assert dc.get_config("scratch_directory") == "/tmp/"
    assert dc.get_config("scratch_directory", hostname="dt5") is None


@pytest.fixture
def opt_state_auto():
    """
    Capture the options state and temporarily override.
    """

    # Snapshot env
    old_globals = copy.deepcopy(dc.config._globals)
    old_environ = dict(os.environ)

    config = {
        "default_compute": {
            "psi_path": "/home/user/psi4",
            "jobs_per_node": 1,
            "nthreads_per_job": "auto",
            "memory_per_job": "auto",
            "scratch_directory": "$TMPDIR"
        }
    }
    

    os.environ["TMPDIR"] = "/tmp/"
    dc.load_options(config)

    yield

    # Reset env
    os.environ.update(old_environ)
    dc.config._globals = old_globals


def test_auto_threads(opt_state_auto):

    assert dc.get_config("jobs_per_node") == 1
    assert isinstance(dc.get_config("nthreads_per_job"), int)
    assert dc.get_config("nthreads_per_job") > 0
    assert dc.get_config("nthreads_per_job") < 100

    assert isinstance(dc.get_config("memory_per_job"), (int, float))
    assert dc.get_config("memory_per_job") > 0.01  # Always more than 1OMB free?


def test_global_repr(opt_state_auto):

    dc.config.global_repr()
