"""
Tests the DQM compute module configuration
"""

import copy
import os

import pytest
import pydantic

import qcengine

from qcengine.testing import environ_context


def test_node_blank():
    node = qcengine.config.NodeDescriptor(name="something", hostname_pattern="*")

    assert node.nthreads_per_job > 0
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
    old_node = copy.deepcopy(qcengine.config.NODE_DESCRIPTORS)

    scratch_name = "myscratch1234"
    with environ_context({"QCA_SCRATCH_DIR": scratch_name}):

        configs = [{
            "name": "default",
            "hostname_pattern": "*",
            "jobs_per_node": 1,
            "nthreads_per_job": 2,
            "memory_per_job": 4,
            "scratch_directory": "$QCA_SCRATCH_DIR"
        }, {
            "name": "dragonstooth",
            "hostname_pattern": "dt*",
            "jobs_per_node": 2,
            "nthreads_per_job": 6,
            "memory_per_job": 60,
            "scratch_directory": "$NOVAR_RANDOM_ABC123"
        }, {
            "name": "newriver",
            "hostname_pattern": "nr*",
            "jobs_per_node": 2,
            "nthreads_per_job": 12,
            "memory_per_job": 120
        }]
        for desc in configs:
            node = qcengine.config.NodeDescriptor(**desc)
            qcengine.config.NODE_DESCRIPTORS[desc["name"]] = node

        yield

        # Reset env
        qcengine.config.NODE_DESCRIPTORS = old_node


def test_node_matching(opt_state_basic):
    node = qcengine.config.get_node_descriptor("nomatching")
    assert node.name == "default"

    node = qcengine.config.get_node_descriptor("dt149")
    assert node.name == "dragonstooth"

    node = qcengine.config.get_node_descriptor("nr149")
    assert node.name == "newriver"


def test_node_env(opt_state_basic):
    node = qcengine.config.get_node_descriptor("dt")
    assert node.name == "dragonstooth"
    assert node.scratch_directory is None

    node = qcengine.config.get_node_descriptor("nomatching")
    assert node.name == "default"
    assert node.scratch_directory == "myscratch1234"


def test_config_default(opt_state_basic):
    config = qcengine.config.get_config("something")
    assert config.nthreads == 2
    assert config.memory == 4

    config = qcengine.config.get_config("dt149")
    assert config.nthreads == 6
    assert config.memory == 60


def test_config_local(opt_state_basic):
    config = qcengine.config.get_config("something", {"nthreads": 10})
    assert config.nthreads == 10
    assert config.memory == 4


def test_config_validation(opt_state_basic):
    with pytest.raises(pydantic.ValidationError):
        config = qcengine.config.get_config(hostname="something", local_options={"bad": 10})


def test_global_repr():

    assert isinstance(qcengine.config.global_repr(), str)
