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


def test_node_auto():

    desc = {
        "name": "something",
        "hostname_pattern": "*",
        "jobs_per_node": 1,
        "ncores": 4,
        "memory": 10,
        "memory_safety_factor": 0,
    }
    node1 = qcengine.config.NodeDescriptor(**desc)
    job1 = qcengine.get_config(hostname=node1)
    assert job1.ncores == 4
    assert pytest.approx(job1.memory) == 10.0

    desc["jobs_per_node"] = 2
    node2 = qcengine.config.NodeDescriptor(**desc)
    job2 = qcengine.get_config(hostname=node2)
    assert job2.ncores == 2
    assert pytest.approx(job2.memory) == 5.0


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
            "name": "dragonstooth",
            "hostname_pattern": "dt*",
            "jobs_per_node": 2,
            "ncores": 12,
            "memory": 120,
            "scratch_directory": "$NOVAR_RANDOM_ABC123"
        }, {
            "name": "newriver",
            "hostname_pattern": "nr*",
            "jobs_per_node": 2,
            "ncores": 24,
            "memory": 240
        }, {
            "name": "default",
            "hostname_pattern": "*",
            "jobs_per_node": 1,
            "memory": 4,
            "memory_safety_factor": 0,
            "ncores": 5,
            "scratch_directory": "$QCA_SCRATCH_DIR"
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
    config = qcengine.config.get_config(hostname="something")
    assert config.ncores == 5
    assert config.memory == 4

    config = qcengine.config.get_config(hostname="dt149")
    assert config.ncores == 6
    assert pytest.approx(config.memory, 0.1) == 54


def test_config_local_ncores(opt_state_basic):
    config = qcengine.config.get_config(hostname="something", local_options={"ncores": 10})
    assert config.ncores == 10
    assert config.memory == 4


def test_config_local_njobs(opt_state_basic):
    config = qcengine.config.get_config(hostname="something", local_options={"jobs_per_node": 5})
    assert config.ncores == 1
    assert pytest.approx(config.memory) == 0.8


def test_config_local_njob_ncore(opt_state_basic):
    config = qcengine.config.get_config(hostname="something", local_options={"jobs_per_node": 3, "ncores": 1})
    assert config.ncores == 1
    assert pytest.approx(config.memory, 0.1) == 1.33


def test_config_local_njob_ncore(opt_state_basic):
    config = qcengine.config.get_config(
        hostname="something", local_options={"jobs_per_node": 3,
                                             "ncores": 1,
                                             "memory": 6})
    assert config.ncores == 1
    assert pytest.approx(config.memory, 0.1) == 6


def test_config_validation(opt_state_basic):
    with pytest.raises(pydantic.ValidationError):
        config = qcengine.config.get_config(hostname="something", local_options={"bad": 10})


def test_global_repr():

    assert isinstance(qcengine.config.global_repr(), str)
