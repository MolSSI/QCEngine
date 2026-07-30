"""Tests for the MPQC4 harness."""

import pytest

import qcengine as qcng
from qcengine.exceptions import InputError
from qcengine.programs.mpqc.germinate import muster_modelchem
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


@pytest.mark.parametrize(
    "method, expected",
    [
        # For HF the wfn is the SCF, so no separate ref block.
        pytest.param("hf", ("SD", False, "Energy", {}), id="hf"),
        pytest.param("mp2", ("MP2", True, "Energy", {"method": "standard"}), id="mp2"),
        pytest.param("CCSD(T)", ("CCSD(T)", True, "Energy", {}), id="case-insensitive"),
        pytest.param("mpqc-mp2", ("MP2", True, "Energy", {"method": "standard"}), id="mpqc-prefix-stripped"),
        pytest.param("cck", ("CCk", True, "Energy", {"k": 2}), id="cck-default-rank"),
        pytest.param("sci", ("sCI", True, "Energy", {}), id="sci-defaults-to-energy"),
        pytest.param("eom-ccsd", ("EOM-CCSD", True, "ExcitationEnergy", {}), id="eom-ccsd"),
        pytest.param("eom-ip-ccsd", ("EOM-IP-CCSD", True, "ExcitationEnergy", {}), id="eom-ip"),
        pytest.param("eom-ea-ccsd", ("EOM-EA-CCSD", True, "ExcitationEnergy", {}), id="eom-ea"),
        pytest.param(
            "eom-cck", ("CCk", True, "ExcitationEnergy", {"k": 2, "eom": {"manifold": "2h2p"}}), id="eom-cck"
        ),
        pytest.param("cis", ("CIS", True, "ExcitationEnergy", {}), id="native-cis"),
        # native passthrough keeps MPQC's exact casing; DF-RHF is its own SCF
        pytest.param("df-rhf", ("DF-RHF", False, "Energy", {}), id="native-scf"),
    ],
)
def test_muster_modelchem(method, expected):
    assert muster_modelchem(method, "energy") == expected
    # driver=properties resolves the same way
    assert muster_modelchem(method, "properties") == expected


@pytest.mark.parametrize(
    "method, requested, expected",
    [
        pytest.param("sci", "ExcitationEnergy", ("sCI", True, "ExcitationEnergy", {}), id="sci-excitation"),
        pytest.param("sci", "Energy", ("sCI", True, "Energy", {}), id="sci-energy-explicit"),
        pytest.param("cck", "ExcitationEnergy", ("CCk", True, "ExcitationEnergy", {"k": 2}), id="cck-excitation"),
        pytest.param(
            "eom-cck", "Energy", ("CCk", True, "Energy", {"k": 2, "eom": {"manifold": "2h2p"}}), id="cck-energy"
        ),
        pytest.param("sCI", "excitationenergy", ("sCI", True, "ExcitationEnergy", {}), id="case-insensitive"),
        pytest.param("hf", None, ("SD", False, "Energy", {}), id="none-takes-default"),
    ],
)
def test_muster_honors_requested_property(method, requested, expected):
    """property__type selects among the properties a wfn type actually provides,
    rather than the method name fixing one."""
    assert muster_modelchem(method, "energy", requested) == expected


@pytest.mark.parametrize(
    "method, requested, match",
    [
        pytest.param("hf", "ExcitationEnergy", "cannot provide ExcitationEnergy", id="scf-has-no-roots"),
        pytest.param("mp2", "ExcitationEnergy", "cannot provide ExcitationEnergy", id="mp2-has-no-roots"),
        pytest.param("ccsd", "ExcitationEnergy", "cannot provide ExcitationEnergy", id="ccsd-has-no-roots"),
        pytest.param("cis", "Energy", "cannot provide Energy", id="cis-has-no-total-energy"),
        pytest.param("eom-ccsd", "Energy", "cannot provide Energy", id="eom-has-no-total-energy"),
        pytest.param("sci", "RDM", "not supported", id="out-of-scope-property"),
    ],
)
def test_muster_rejects_unavailable_property(method, requested, match):
    with pytest.raises(InputError, match=match):
        muster_modelchem(method, "energy", requested)


@pytest.mark.parametrize(
    "method",
    [
        # unknown: test_compute_bad_models requires InputError before any subprocess
        "bad",
        # complex/periodic. A leading-z deny-pattern would have leaked DFJ-zRHF
        # and friends, hence the allow-list.
        "zrhf",
        "zmp2",
        "zcck",
        "zsci",
        "zsd",
        "dfj-zrhf",
        "maj-cadfk-zrhf",
        "gammapointmp2",
        # basis-free MRA
        "mra::spsolver",
        # spin-orbital hard-coded CC
        "ccsd-so",
        "ccsd_so-gpu",
    ],
)
def test_muster_rejects_out_of_scope(method):
    with pytest.raises(InputError):
        muster_modelchem(method, "energy")


@pytest.mark.parametrize("driver", ["gradient", "hessian"])
def test_muster_rejects_derivative_drivers(driver):
    """The message must contain `gradient not implemented`."""
    with pytest.raises(InputError, match=f"{driver} not implemented"):
        muster_modelchem("hf", driver)
