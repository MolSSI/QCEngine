"""Tests for the MPQC4 harness."""

import pytest

import qcengine as qcng
from qcengine.programs.mpqc.keywords import deep_merge, extract_reserved, format_keywords
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


@pytest.mark.parametrize(
    "flat, expected",
    [
        pytest.param({"scf__max_iter": 100}, {"scf": {"max_iter": 100}}, id="one-level"),
        pytest.param({"wfn__eom__manifold": "2h1p"}, {"wfn": {"eom": {"manifold": "2h1p"}}}, id="two-level"),
        pytest.param(
            {"scf__max_iter": 100, "scf__spin_restricted": False},
            {"scf": {"max_iter": 100, "spin_restricted": False}},
            id="siblings-merge",
        ),
        pytest.param({"units": "2018CODATA"}, {"units": "2018CODATA"}, id="bare-key"),
    ],
)
def test_format_keywords(flat, expected):
    assert format_keywords(flat) == expected


def test_deep_merge_and_extract_reserved():
    """deep_merge: overlay wins, neither argument mutated. extract_reserved:
    the three reserved channels leave `remaining`, and `property__*` in
    particular must not survive into the tree as a stray top-level `property`
    block, which MPQC reads only from the task subtree."""
    base, overlay = {"scf": {"type": "SD", "max_iter": 30}}, {"scf": {"max_iter": 100}}
    assert deep_merge(base, overlay) == {"scf": {"type": "SD", "max_iter": 100}}
    assert (base, overlay) == ({"scf": {"type": "SD", "max_iter": 30}}, {"scf": {"max_iter": 100}})

    assert extract_reserved(
        {
            "scf__max_iter": 100,
            "mpqc_input": {"atoms": {}},
            "mpqc_env": {"MAD_BUFFER_SIZE": "64MB"},
            "property__n_roots": 4,
        }
    ) == ({"scf__max_iter": 100}, {"atoms": {}}, {"MAD_BUFFER_SIZE": "64MB"}, {"n_roots": 4})

    assert extract_reserved({"scf__max_iter": 100}) == ({"scf__max_iter": 100}, {}, {}, {})
