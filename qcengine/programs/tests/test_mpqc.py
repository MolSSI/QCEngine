"""Tests for the MPQC4 harness."""

import json

import pytest
from qcelemental.models.v2 import AtomicInput, AtomicSpecification, Molecule
from qcelemental.tests.test_model_results import center_data

import qcengine as qcng
from qcengine.config import get_config
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


def _atomic_input(method="hf", basis="6-31G", driver="energy", keywords=None, molecule=None):
    return AtomicInput(
        molecule=molecule or Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
        specification=AtomicSpecification(
            driver=driver, model={"method": method, "basis": basis}, keywords=keywords or {}
        ),
    )


def _built_tree(inp):
    harness = qcng.get_program("mpqc", check=False)
    job = harness.build_input(inp, get_config(task_config={"ncores": 4, "memory": 2.0}))
    assert job["command"][-2:] == ["-i", "mpqc.json"]
    assert "-o" not in job["command"]  # no -o, so results land on stdout
    return json.loads(job["infiles"]["mpqc.json"])


_ESCAPE_HATCH_TREE = {
    "units": "2018CODATA",
    "obs": {"name": "sto-3g", "atoms": "$:atoms"},
    "wfn_world": {"atoms": "$:atoms", "basis": "$:obs"},
    "wfn": {"type": "SD", "wfn_world": "$:wfn_world", "atoms": "$:atoms"},
}

_OH_RADICAL = Molecule(
    symbols=["O", "H"], geometry=[0.0, 0.0, 0.0, 0.0, 0.0, 1.8], molecular_charge=0.0, molecular_multiplicity=2
)


@pytest.mark.parametrize(
    "kwargs, expected",
    [
        # Exactly one property, singular `property` key, nested under `mpqc`.
        pytest.param({}, {"mpqc": {"property": {"type": "Energy"}}, "property": None}, id="single-property-block"),
        # geometry in bohr and unsorted, because the atoms block defaults to
        # angstrom and would otherwise reorder atoms
        pytest.param(
            {},
            {"atoms": {"units": "bohr", "sort_input": False, "atoms": [{"element": "H"}, {"element": "H"}]}},
            id="molecule-bohr-unsorted",
        ),
        # HF: the wfn is the SCF, no separate scf block and no ref
        pytest.param({"method": "hf"}, {"scf": None, "wfn": {"type": "SD", "ref": None}}, id="hf-self-referencing"),
        pytest.param(
            {"method": "mp2"},
            {"scf": {"type": "SD", "spin_restricted": True}, "wfn": {"type": "MP2", "ref": "$:scf"}},
            id="mp2-references-scf",
        ),
        pytest.param(
            {"basis": "cc-pVDZ"},
            {"obs": {"name": "cc-pVDZ", "atoms": "$:atoms"}, "wfn_world": {"atoms": "$:atoms", "basis": "$:obs"}},
            id="basis-and-world",
        ),
        # multiplet -> UHF by default
        pytest.param(
            {"method": "mp2", "molecule": _OH_RADICAL},
            {"scf": {"charge": 0, "multiplicity": 2, "spin_restricted": False}},
            id="charge-and-multiplicity",
        ),
        # user keyword wins, generated sibling keys survive the merge
        pytest.param(
            {"method": "mp2", "keywords": {"scf__max_iter": 250}},
            {"scf": {"max_iter": 250, "type": "SD"}},
            id="keywords-override",
        ),
        # sCI provides both properties; property__type picks one and
        # property__* keys land in the property block alongside it
        pytest.param(
            {"method": "sci", "keywords": {"property__type": "ExcitationEnergy", "property__n_roots": 4}},
            {"mpqc": {"property": {"type": "ExcitationEnergy", "n_roots": 4}}},
            id="property-keywords",
        ),
        # without property__type, sCI is a ground-state energy
        pytest.param(
            {"method": "sci"},
            {"mpqc": {"property": {"type": "Energy", "n_roots": None}}},
            id="sci-defaults-to-energy",
        ),
        # MPQC's own ExcitationEnergy default
        pytest.param({"method": "eom-ccsd"}, {"mpqc": {"property": {"n_roots": 3}}}, id="default-n-roots"),
        pytest.param({"method": "hf"}, {"mpqc": {"property": {"n_roots": None}}}, id="energy-has-no-n-roots"),
        # escape hatch: user tree preserved, molecule and property still injected
        pytest.param(
            {"keywords": {"mpqc_input": _ESCAPE_HATCH_TREE}},
            {"obs": {"name": "sto-3g"}, "atoms": {"units": "bohr"}, "mpqc": {"property": {"type": "Energy"}}},
            id="mpqc_input-escape-hatch",
        ),
        # a harness-owned property block must not erase user mpqc__* keywords
        pytest.param(
            {"keywords": {"mpqc__file_prefix": "job"}},
            {"mpqc": {"file_prefix": "job", "property": {"type": "Energy"}}},
            id="mpqc-subtree-preserved",
        ),
        pytest.param(
            {"method": "mp2", "keywords": {"wfn__method": "df", "dfbs__name": "cc-pVDZ-RI"}},
            {"dfbs": {"name": "cc-pVDZ-RI"}, "wfn_world": {"df_basis": "$:dfbs"}},
            id="df-with-df-basis",
        ),
        # the constants set is a default, not a fixture
        pytest.param({}, {"units": "2018CODATA"}, id="units-default"),
        pytest.param({"keywords": {"units": "2014CODATA"}}, {"units": "2014CODATA"}, id="units-override"),
        pytest.param(
            {"keywords": {"mpqc_input": {**_ESCAPE_HATCH_TREE, "units": "2014CODATA"}}},
            {"units": "2014CODATA"},
            id="units-override-escape-hatch",
        ),
    ],
)
def test_build_input_tree(kwargs, expected, request):
    _assert_subtree(_built_tree(_atomic_input(**kwargs)), expected, request.node.name)


def _assert_subtree(actual, expected, tnm):
    """Compare a partial tree: a None leaf asserts the key is absent, a dict
    recurses, a list recurses element-wise, anything else compares equal."""
    for key, want in expected.items():
        if want is None:
            assert key not in actual, f"{tnm}: {key} should be absent, found {actual[key]!r}"
            continue
        assert key in actual, f"{tnm}: {key} missing"
        got = actual[key]
        if isinstance(want, dict):
            _assert_subtree(got, want, tnm)
        elif isinstance(want, list):
            assert len(got) == len(want), f"{tnm}: {key} has {len(got)} items, want {len(want)}"
            for g, w in zip(got, want):
                _assert_subtree(g, w, tnm)
        else:
            assert got == want, f"{tnm}: {key} is {got!r}, want {want!r}"


@pytest.mark.parametrize(
    "kwargs, match",
    [
        # exact string test_compute_energy_qcsk_basis greps for
        pytest.param(
            {"basis": {"name": "custom", "center_data": center_data, "atom_map": ["bs_sto3g_h", "bs_sto3g_h"]}},
            "QCSchema BasisSet for model.basis not implemented",
            id="qcschema-basisset",
        ),
        pytest.param({"basis": None}, "basis", id="no-basis"),
        pytest.param(
            {"molecule": Molecule(symbols=["H", "H"], geometry=[0.0, 0.0, 0.0, 0.0, 0.0, 1.4], real=[True, False])},
            "[Gg]host",
            id="ghost-atoms",
        ),
        # a df method with no dfbs leaves a dangling $:dfbs
        pytest.param({"method": "mp2", "keywords": {"wfn__method": "df"}}, "dfbs__name", id="df-without-df-basis"),
    ],
)
def test_build_input_rejects(kwargs, match):
    with pytest.raises(InputError, match=match):
        _built_tree(_atomic_input(**kwargs))
