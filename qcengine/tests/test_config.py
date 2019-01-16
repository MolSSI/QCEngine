"""
Tests the DQM compute module configuration
"""

import copy
import os

import pytest

import qcengine

from qcengine.testing import environ_context


def test_node_blank():
    node = qcengine.config.NodeDescriptor(name="something", hostname_pattern="*")

    assert node.nthreads_per_job > 1
    assert node.available_memory > 0.1
    assert node.memory_per_job > 0.1


def test_node_auto():

    description = {
        "name": "something",
        "hostname_pattern": "*",
        "jobs_per_node": 1,
        "nthreads_per_job": "auto",
        "memory_per_job": "auto",
        "total_cores": 4,
        "available_memory": 10,
        "memory_safety_factor": 0,
    }

    node1 = qcengine.config.NodeDescriptor(**description)
    assert node1.nthreads_per_job == 4
    assert pytest.approx(node1.memory_per_job) == 10.0

    description["jobs_per_node"] = 2
    node2 = qcengine.config.NodeDescriptor(**description)
    assert node2.nthreads_per_job == 2
    assert pytest.approx(node2.memory_per_job) == 5


def test_node_environ():

    scratch_name = "myscratch1234"
    with environ_context({"QCA_SCRATCH_DIR": scratch_name}):
        description = {
            "name": "something",
            "hostname_pattern": "*",
            "scratch_directory": "$QCA_SCRATCH_DIR",
        }

        node = qcengine.config.NodeDescriptor(**description)
        assert node.scratch_directory == scratch_name


def test_node_skip_environ():
    description = {
        "name": "something",
        "hostname_pattern": "*",
        "scratch_directory": "$RANDOM_ENVIRON",
    }

    node = qcengine.config.NodeDescriptor(**description)
    assert node.scratch_directory is None


@pytest.fixture
def opt_state_basic():
    """
    Capture the options state and temporarily override.
    """

    # Snapshot env
    old_node = copy.deepcopy(dc.config.NODE_DESCRIPTORS)

    scratch_name = "myscratch1234"
    with environ_context({"QCA_SCRATCH_DIR": scratch_name}):

        configs = [{
                "name": "default",
                "hostname_pattern": "*",
                "jobs_per_node": 1,
                "nthreads_per_job": 2,
                "memory_per_job": 4,
                "scratch_directory": "$QCA_SCRATCH_DIR"
            },{
                "name": "dragonsooth":,
                "hostname_pattern": "dt*",
                "jobs_per_node": 2,
                "nthreads_per_job": 6,
                "memory_per_job": 60,
                "scratch_directory": "$NOVAR_RANDOM_ABC123"
            },{
                "name": "newriver",
                "hostname_pattern": "nr*",
                "jobs_per_node": 2,
                "nthreads_per_job": 12,
                "memory_per_job": 120
            }
        ]
        for desc in configs:
            node = qcengine.config.NodeDescriptor(**descc)
            dc.config.NODE_DESCRIPTORS[desc["name"]] = node

        yield

        # Reset env
        dc.config.NODE_DESCRIPTORS = old_node


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


def test_global_repr(opt_state_auto):

    assert isinstance(dc.config.global_repr(), str)
