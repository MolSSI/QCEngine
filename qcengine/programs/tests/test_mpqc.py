"""Tests for the MPQC4 harness."""

import pytest

import qcengine as qcng
from qcengine.testing import uusing


def test_mpqc_registration():
    """Registered under the lowercase name users pass to compute(), with the
    resource declarations qcengine.config uses to allocate cores and memory."""
    assert "mpqc" in qcng.list_all_programs()

    harness = qcng.get_program("mpqc", check=False)
    assert harness.name == "MPQC"
    assert isinstance(harness.found(), bool)  # answers without raising when absent
    assert (harness.thread_parallel, harness.node_parallel, harness.managed_memory) == (True, True, True)


@uusing("mpqc")
def test_mpqc_version():
    """get_version parses `mpqc -v` output into a comparable version."""
    ver = qcng.get_program("mpqc").get_version()
    assert ver.startswith("4."), ver
