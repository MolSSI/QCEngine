"""Tests for the Gaussian harness: germinate, keywords, harvester, and runner."""

import importlib.util
import re
from decimal import Decimal
from typing import Any, Dict
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.models.v2 import AtomicInput, AtomicSpecification, BasisSet, Molecule

from qcengine.exceptions import InputError, ResourceError, UnknownError
from qcengine.programs.gaussian.germinate import muster_modelchem
from qcengine.programs.gaussian.harvester import harvest, is_normal_termination
from qcengine.programs.gaussian.keywords import build_com_file, build_route_line
from qcengine.programs.gaussian import runner as gaussian_runner
from qcengine.programs.gaussian.runner import GaussianHarness
from qcengine.testing import uusing as using

_cclib_available = importlib.util.find_spec("cclib") is not None
requires_cclib = pytest.mark.skipif(not _cclib_available, reason="cclib not installed")

# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

_EV_TO_HARTREE = qcel.constants.conversion_factor("eV", "hartree")


@pytest.fixture
def water():
    """Simple water molecule for use in tests."""
    return Molecule.from_data(
        """
        0 1
        O  0.000000000  0.000000000  0.117176320
        H  0.000000000  0.757335953 -0.468705281
        H  0.000000000 -0.757335953 -0.468705281
        units angstrom
        """
    )


@pytest.fixture
def water_energy_input(water):
    """Minimal AtomicInput for a closed-shell HF/STO-3G energy on water."""
    return AtomicInput(molecule=water, specification=AtomicSpecification(program="gaussian", driver="energy", model={"method": "HF", "basis": "STO-3G"}))


@pytest.fixture
def water_energy_input_all_files(water):
    """AtomicInput with native_files protocol set to 'all' for testing file output."""
    return AtomicInput(
        molecule=water,
        specification=AtomicSpecification(
            program="gaussian",
            driver="energy",
            model={"method": "HF", "basis": "STO-3G"},
            protocols={"native_files": "all"},
        ),
    )


@pytest.fixture
def harness():
    return GaussianHarness()


# ---------------------------------------------------------------------------
# Mock ccdata helper
# ---------------------------------------------------------------------------


def _make_ccdata(
    scfenergies=None,
    mpenergies=None,
    ccenergies=None,
    atomcoords=None,
    atomnos=None,
    grads=None,
    hessian=None,
    moments=None,
    atomcharges=None,
):
    """Build a minimal mock ccdata object."""
    # Water geometry in Angstrom
    if atomcoords is None:
        atomcoords = np.array(
            [[[0.0, 0.0, 0.117176], [0.0, 0.757336, -0.468705], [0.0, -0.757336, -0.468705]]]
        )
    if atomnos is None:
        atomnos = np.array([8, 1, 1])
    if scfenergies is None:
        scfenergies = np.array([-2041.3])

    obj = MagicMock()
    obj.scfenergies = scfenergies
    obj.mpenergies = mpenergies
    obj.ccenergies = ccenergies
    obj.atomcoords = atomcoords
    obj.atomnos = atomnos
    obj.grads = grads
    obj.hessian = hessian
    obj.moments = moments
    obj.atomcharges = atomcharges if atomcharges is not None else {}
    # Make hasattr work correctly
    if grads is None:
        del obj.grads
    if hessian is None:
        del obj.hessian
    if moments is None:
        del obj.moments
    return obj


def _make_ccdata_with_attrs(**kwargs):
    """Like _make_ccdata but use spec=None so we can control hasattr."""
    ccdata = type("CCData", (), {})()
    defaults = dict(
        scfenergies=np.array([-2041.3]),
        atomcoords=np.array(
            [[[0.0, 0.0, 0.117176], [0.0, 0.757336, -0.468705], [0.0, -0.757336, -0.468705]]]
        ),
        atomnos=np.array([8, 1, 1]),
        mpenergies=None,
        ccenergies=None,
        grads=None,
        hessian=None,
        moments=None,
        atomcharges={},
    )
    defaults.update(kwargs)
    for k, v in defaults.items():
        setattr(ccdata, k, v)
    # cclib's Gaussian parser sets `natom` as an int. Mirror that here so the
    # harvester's `ccdata.natom` reads work without each test having to set it.
    if not hasattr(ccdata, "natom"):
        setattr(ccdata, "natom", len(ccdata.atomnos))
    return ccdata


def _patch_cclib(ccdata):
    """Return a context manager that patches cclib to return ccdata."""
    mock_ccopen = MagicMock()
    mock_ccopen.return_value.parse.return_value = ccdata
    return patch("cclib.io.ccopen", mock_ccopen)


# ---------------------------------------------------------------------------
# Tests: germinate.muster_modelchem — method mapping
# ---------------------------------------------------------------------------


# rq-0199bb02
def test_muster_modelchem_hf():
    result = muster_modelchem("hf", "energy", 1)
    assert result["method_string"] == "HF"
    assert result["job_type"] == ""
    assert result["extra_keywords"] == {}


# rq-7049aaa5
def test_muster_modelchem_scf_alias():
    result = muster_modelchem("scf", "energy", 1)
    assert result["method_string"] == "HF"


# rq-546faafb
def test_muster_modelchem_mp2():
    result = muster_modelchem("mp2", "energy", 1)
    assert result["method_string"] == "MP2"


# rq-2052a2fa
def test_muster_modelchem_ccsd():
    result = muster_modelchem("ccsd", "energy", 1)
    assert result["method_string"] == "CCSD"


# rq-92492893
def test_muster_modelchem_ccsd_t():
    result = muster_modelchem("ccsd(t)", "energy", 1)
    assert result["method_string"] == "CCSD(T)"


# rq-851d782c
def test_muster_modelchem_passthrough():
    result = muster_modelchem("B3LYP", "energy", 1)
    assert result["method_string"] == "B3LYP"


# rq-079b1e22
def test_muster_modelchem_case_insensitive():
    result = muster_modelchem("CCSD(T)", "energy", 1)
    assert result["method_string"] == "CCSD(T)"


# ---------------------------------------------------------------------------
# Tests: germinate.muster_modelchem — driver mapping
# ---------------------------------------------------------------------------


# rq-c266b9bb
def test_muster_modelchem_energy_job_type():
    result = muster_modelchem("hf", "energy", 1)
    assert result["job_type"] == ""


# rq-f1c116db
def test_muster_modelchem_gradient_job_type():
    result = muster_modelchem("hf", "gradient", 1)
    assert result["job_type"] == "Force=NoStep"


# rq-82942b90
def test_muster_modelchem_hessian_job_type():
    result = muster_modelchem("hf", "hessian", 1)
    assert result["job_type"] == "Freq"


# rq-f512b60b
def test_muster_modelchem_properties_driver():
    result = muster_modelchem("hf", "properties", 1)
    assert result["job_type"] == ""
    assert result["extra_keywords"].get("Pop") == "Full"


# rq-e008c6fe
def test_muster_modelchem_unsupported_driver():
    with pytest.raises(InputError):
        muster_modelchem("hf", "optimization", 1)


# ---------------------------------------------------------------------------
# Tests: germinate.muster_modelchem — open-shell handling
# ---------------------------------------------------------------------------


# rq-972b6470
def test_muster_modelchem_singlet_no_u_prefix():
    result = muster_modelchem("hf", "energy", 1)
    assert result["method_string"] == "HF"


# rq-6d73bb9e
def test_muster_modelchem_doublet_auto_u_prefix():
    result = muster_modelchem("hf", "energy", 2)
    assert result["method_string"] == "UHF"


# rq-c946c285
def test_muster_modelchem_triplet_dft_auto_u_prefix():
    result = muster_modelchem("B3LYP", "energy", 3)
    assert result["method_string"] == "UB3LYP"


# rq-25c4ad49
def test_muster_modelchem_no_double_u_prefix():
    result = muster_modelchem("UHF", "energy", 2)
    assert result["method_string"] == "UHF"


# rq-fe13245b
def test_muster_modelchem_rohf_no_u_prefix():
    result = muster_modelchem("ROHF", "energy", 2)
    assert result["method_string"] == "ROHF"


# ---------------------------------------------------------------------------
# Tests: keywords.build_route_line
# ---------------------------------------------------------------------------


# rq-af032242
def test_build_route_line_hf_energy():
    line = build_route_line("HF", "STO-3G", "", {}, {})
    assert line == "#P HF/STO-3G"


# rq-9acd1ce5
def test_build_route_line_gradient():
    line = build_route_line("HF", "STO-3G", "Force=NoStep", {}, {})
    assert line.startswith("#P HF/STO-3G")
    assert "Force=NoStep" in line


# rq-cd0d6637
def test_build_route_line_properties_pop_full():
    line = build_route_line("HF", "STO-3G", "", {"Pop": "Full"}, {})
    assert "Pop=Full" in line


# rq-7152455e
def test_build_route_line_user_keyword_appended():
    line = build_route_line("HF", "STO-3G", "", {}, {"SCF": "Tight"})
    assert "SCF=Tight" in line


# rq-e7897f09
def test_build_route_line_reserved_keywords_excluded():
    line = build_route_line("HF", "STO-3G", "", {}, {"memory": "8GB", "nprocs": "4"})
    assert "memory" not in line
    assert "nprocs" not in line


# ---------------------------------------------------------------------------
# Tests: keywords.build_com_file
# ---------------------------------------------------------------------------


# rq-d7e8161d
def test_build_com_file_link0_section():
    com = build_com_file(
        {"NProcShared": 4, "Mem": "8GB"},
        "#P HF/STO-3G",
        "title",
        0,
        1,
        "H  0.0  0.0  0.0\nH  0.0  0.0  0.74",
    )
    assert com.startswith("%NProcShared=4\n%Mem=8GB")


# rq-342a2ca3
def test_build_com_file_charge_multiplicity_line():
    com = build_com_file({"NProcShared": 1, "Mem": "1GB"}, "#P HF/STO-3G", "t", 1, 2, "H 0 0 0")
    assert "1 2" in com


# rq-b819d8e0
def test_build_com_file_coordinate_precision():
    atom_block = "O   0.0000000000   0.1171763200   0.0000000000"
    com = build_com_file({"NProcShared": 1, "Mem": "1GB"}, "#P HF/STO-3G", "t", 0, 1, atom_block)
    # Check there are at least 10 decimal digits in the atom line
    assert re.search(r"\d+\.\d{10}", com)


# rq-23d8b2c6
def test_build_com_file_trailing_blank_line():
    com = build_com_file({"NProcShared": 1, "Mem": "1GB"}, "#P HF/STO-3G", "t", 0, 1, "H 0 0 0")
    assert com.endswith("\n\n")


# rq-db92093a
def test_build_input_memory_rounds_down_minimum_1(water_energy_input, harness):
    from qcengine.config import TaskConfig

    cfg = TaskConfig(scratch_directory=None, scratch_messy=False, ncores=1, nnodes=1, memory=0.5, retries=0, mpiexec_command=None)
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/usr/local/g16/g16"):
        job = harness.build_input(water_energy_input, cfg)
    assert "%Mem=1GB" in job["infiles"]["gaussian.com"]


# rq-4ad16d95
def test_build_input_memory_whole_gb(water_energy_input, harness):
    from qcengine.config import TaskConfig

    cfg = TaskConfig(scratch_directory=None, scratch_messy=False, ncores=4, nnodes=1, memory=8.0, retries=0, mpiexec_command=None)
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/usr/local/g16/g16"):
        job = harness.build_input(water_energy_input, cfg)
    assert "%Mem=8GB" in job["infiles"]["gaussian.com"]


# rq-0aa99076
def test_build_input_basisset_object_raises(water, harness):
    from qcengine.config import TaskConfig

    bs = BasisSet(name="sto-3g", center_data={}, atom_map=[])
    inp = AtomicInput(molecule=water, specification=AtomicSpecification(program="gaussian", driver="energy", model={"method": "HF", "basis": bs}))
    cfg = TaskConfig(scratch_directory=None, scratch_messy=False, ncores=1, nnodes=1, memory=1.0, retries=0, mpiexec_command=None)
    with pytest.raises(InputError):
        harness.build_input(inp, cfg)


# ---------------------------------------------------------------------------
# Tests: runner.GaussianHarness.found()
# ---------------------------------------------------------------------------


# rq-4cf7c460
def test_found_true_when_g16_and_cclib_available(harness):
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/usr/local/g16/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            assert harness.found() is True


# rq-d3b8d8c3
def test_found_falls_back_to_g09(harness):
    def fake_find():
        # Simulates g16 absent, g09 present
        from unittest.mock import call
        return "/usr/local/g09/g09"

    with patch("qcengine.programs.gaussian.runner._find_gaussian", side_effect=fake_find):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            assert harness.found() is True


# rq-8e21b7d7
def test_found_false_when_no_gaussian(harness):
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value=None):
        assert harness.found() is False


# rq-b5e72c6d
def test_found_false_when_cclib_missing(harness):
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/usr/local/g16/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=False):
            assert harness.found() is False


# rq-bd98d417
def test_found_raises_when_gaussian_missing_raise_error(harness):
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value=None):
        with pytest.raises(ModuleNotFoundError):
            harness.found(raise_error=True)


# rq-242ef1d8
def test_found_raises_when_cclib_missing_raise_error(harness):
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/usr/local/g16/g16"):
        with patch(
            "qcengine.programs.gaussian.runner.which_import",
            side_effect=ModuleNotFoundError("cclib not found"),
        ):
            with pytest.raises(ModuleNotFoundError):
                harness.found(raise_error=True)


# ---------------------------------------------------------------------------
# Tests: runner.GaussianHarness.get_version()
# ---------------------------------------------------------------------------

_VERSION_BANNER_16 = (
    " Gaussian 16, Revision C.02,\n"
    " M. J. Frisch, et al.\n"
    " Copyright (c) 1988-2017, Gaussian, Inc.\n"
)

_VERSION_BANNER_09 = (
    " Gaussian 09, Revision D.01,\n"
    " M. J. Frisch, et al.\n"
)


# rq-186c5060
def test_get_version_g16(harness):
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch(
                "qcengine.programs.gaussian.runner.execute",
                return_value=(True, {"stdout": _VERSION_BANNER_16, "stderr": "", "outfiles": {}}),
            ):
                harness.version_cache.clear()
                version = harness.get_version()
    assert version == "16.C.02"


# rq-850ede4f
def test_get_version_g09(harness):
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g09"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch(
                "qcengine.programs.gaussian.runner.execute",
                return_value=(True, {"stdout": _VERSION_BANNER_09, "stderr": "", "outfiles": {}}),
            ):
                harness.version_cache.clear()
                version = harness.get_version()
    assert version == "09.D.01"


# rq-3e134642
def test_get_version_cached(harness):
    harness.version_cache["/g16"] = "16.C.02"
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch("qcengine.programs.gaussian.runner.execute") as mock_exe:
                version = harness.get_version()
                mock_exe.assert_not_called()
    assert version == "16.C.02"
    harness.version_cache.clear()


# rq-20d9a784
def test_get_version_unparseable_raises(harness):
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch(
                "qcengine.programs.gaussian.runner.execute",
                return_value=(False, {"stdout": "no version here", "stderr": "", "outfiles": {}}),
            ):
                harness.version_cache.clear()
                with pytest.raises(UnknownError):
                    harness.get_version()


# ---------------------------------------------------------------------------
# Tests: runner.GaussianHarness.compute() — error dispatch
# ---------------------------------------------------------------------------

_NORMAL_LOG = "Normal termination of Gaussian 16\n"


def _mock_compute(harness, water_energy_input, log_text):
    """Helper: patch found, build_input, execute, and get_version for compute tests."""
    job_record = {
        "infiles": {"gaussian.com": "%NProcShared=1\n#P HF/STO-3G\n"},
        "command": ["/g16", "gaussian.com"],
        "scratch_directory": None,
        "scratch_messy": False,
    }
    dexe = {
        "outfiles": {"gaussian.log": log_text},
        "stdout": "",
        "stderr": "",
    }
    return job_record, dexe


# rq-c77a2c13
def test_compute_raises_for_basisset_object(harness, water):
    bs = BasisSet(name="sto-3g", center_data={}, atom_map=[])
    inp = AtomicInput(molecule=water, specification=AtomicSpecification(program="gaussian", driver="energy", model={"method": "HF", "basis": bs}))
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with pytest.raises(InputError):
                harness.compute(inp, MagicMock())


# rq-a4c449da
def test_compute_raises_resource_error_on_memory_failure(harness, water_energy_input):
    log = "galloc: could not allocate memory\nError termination\n"
    job_record, dexe = _mock_compute(harness, water_energy_input, log)
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch.object(GaussianHarness, "build_input", return_value=job_record):
                with patch.object(GaussianHarness, "execute", return_value=(False, dexe)):
                    with pytest.raises(ResourceError):
                        harness.compute(water_energy_input, MagicMock())


# rq-165224fc
def test_compute_raises_input_error_on_ernie(harness, water_energy_input):
    log = "Ernie: is not defined as a Gaussian input file\nError termination\n"
    job_record, dexe = _mock_compute(harness, water_energy_input, log)
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch.object(GaussianHarness, "build_input", return_value=job_record):
                with patch.object(GaussianHarness, "execute", return_value=(False, dexe)):
                    with pytest.raises(InputError):
                        harness.compute(water_energy_input, MagicMock())


# rq-edb73b83
def test_compute_raises_input_error_on_internal_coord(harness, water_energy_input):
    log = "Error in internal coordinate system\nError termination\n"
    job_record, dexe = _mock_compute(harness, water_energy_input, log)
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch.object(GaussianHarness, "build_input", return_value=job_record):
                with patch.object(GaussianHarness, "execute", return_value=(False, dexe)):
                    with pytest.raises(InputError):
                        harness.compute(water_energy_input, MagicMock())


# rq-71aceb93
def test_compute_raises_unknown_error_on_abnormal_termination(harness, water_energy_input):
    log = "Something went wrong.\nError termination via Lnk1e\n"
    job_record, dexe = _mock_compute(harness, water_energy_input, log)
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch.object(GaussianHarness, "build_input", return_value=job_record):
                with patch.object(GaussianHarness, "execute", return_value=(False, dexe)):
                    with pytest.raises(UnknownError):
                        harness.compute(water_energy_input, MagicMock())


# ---------------------------------------------------------------------------
# Tests: runner.GaussianHarness.compute() — success paths
# ---------------------------------------------------------------------------

def _full_compute_patches(harness, water_energy_input, log_text, ccdata, return_result):
    """Return context-manager stack for a successful compute() call."""
    from contextlib import ExitStack

    job_record = {
        "infiles": {"gaussian.com": "%NProcShared=1\n#P HF/STO-3G\n"},
        "command": ["/g16", "gaussian.com"],
        "scratch_directory": None,
        "scratch_messy": False,
    }
    dexe = {
        "outfiles": {"gaussian.log": log_text},
        "stdout": "",
        "stderr": "",
    }

    stack = ExitStack()
    stack.enter_context(patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"))
    stack.enter_context(patch("qcengine.programs.gaussian.runner.which_import", return_value=True))
    stack.enter_context(patch.object(GaussianHarness, "build_input", return_value=job_record))
    stack.enter_context(patch.object(GaussianHarness, "execute", return_value=(True, dexe)))
    stack.enter_context(patch.object(GaussianHarness, "get_version", return_value="16.C.02"))
    stack.enter_context(_patch_cclib(ccdata))
    return stack


# rq-99807237
@requires_cclib
def test_compute_returns_atomic_result_energy(harness, water_energy_input_all_files, water):
    log_text = _NORMAL_LOG
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))

    with _full_compute_patches(harness, water_energy_input_all_files, log_text, ccdata, None):
        result = harness.compute(water_energy_input_all_files, MagicMock())

    assert result.success is True
    assert isinstance(result.return_result, float)
    assert "input" in result.native_files  # rq-f13e0d24
    assert "gaussian.log" in result.native_files
    assert result.provenance.creator == "Gaussian"


# rq-e62578e2
@requires_cclib
def test_compute_returns_atomic_result_gradient(harness, water):
    forces = np.array([[0.01, 0.0, -0.02], [0.0, 0.005, 0.01], [0.0, -0.005, 0.01]])
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        grads=forces[np.newaxis, :, :],  # shape (1, 3, 3)
    )
    inp = AtomicInput(molecule=water, specification=AtomicSpecification(program="gaussian", driver="gradient", model={"method": "HF", "basis": "STO-3G"}))
    log_text = _NORMAL_LOG

    job_record = {
        "infiles": {"gaussian.com": ""},
        "command": ["/g16", "gaussian.com"],
        "scratch_directory": None,
        "scratch_messy": False,
    }
    dexe = {"outfiles": {"gaussian.log": log_text}, "stdout": "", "stderr": ""}

    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch.object(GaussianHarness, "build_input", return_value=job_record):
                with patch.object(GaussianHarness, "execute", return_value=(True, dexe)):
                    with patch.object(GaussianHarness, "get_version", return_value="16.C.02"):
                        with _patch_cclib(ccdata):
                            result = harness.compute(inp, MagicMock())

    assert result.success is True
    # AtomicResult shapes gradient as (natoms, 3); 3 atoms → shape (3, 3)
    assert np.array(result.return_result).shape == (3, 3)


# rq-4d683dac
@requires_cclib
def test_compute_returns_atomic_result_hessian(harness, water):
    hess_matrix = np.eye(9) * 0.5
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        hessian=hess_matrix,
    )
    inp = AtomicInput(molecule=water, specification=AtomicSpecification(program="gaussian", driver="hessian", model={"method": "HF", "basis": "STO-3G"}))
    log_text = _NORMAL_LOG

    job_record = {
        "infiles": {"gaussian.com": ""},
        "command": ["/g16", "gaussian.com"],
        "scratch_directory": None,
        "scratch_messy": False,
    }
    dexe = {"outfiles": {"gaussian.log": log_text}, "stdout": "", "stderr": ""}

    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch.object(GaussianHarness, "build_input", return_value=job_record):
                with patch.object(GaussianHarness, "execute", return_value=(True, dexe)):
                    with patch.object(GaussianHarness, "get_version", return_value="16.C.02"):
                        with _patch_cclib(ccdata):
                            result = harness.compute(inp, MagicMock())

    assert result.success is True
    # AtomicResult shapes hessian as (3*natoms, 3*natoms); 3 atoms → (9, 9)
    assert np.array(result.return_result).shape == (9, 9)


# rq-1a051bc2
@requires_cclib
def test_compute_returns_atomic_result_properties(harness, water):
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        moments=[None, np.array([0.0, 1.85, 0.0])],
        atomcharges={"mulliken": [-0.5, 0.25, 0.25]},
    )
    inp = AtomicInput(molecule=water, specification=AtomicSpecification(program="gaussian", driver="properties", model={"method": "HF", "basis": "STO-3G"}))
    log_text = _NORMAL_LOG

    job_record = {
        "infiles": {"gaussian.com": ""},
        "command": ["/g16", "gaussian.com"],
        "scratch_directory": None,
        "scratch_messy": False,
    }
    dexe = {"outfiles": {"gaussian.log": log_text}, "stdout": "", "stderr": ""}

    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch.object(GaussianHarness, "build_input", return_value=job_record):
                with patch.object(GaussianHarness, "execute", return_value=(True, dexe)):
                    with patch.object(GaussianHarness, "get_version", return_value="16.C.02"):
                        with _patch_cclib(ccdata):
                            result = harness.compute(inp, MagicMock())

    assert result.success is True
    assert isinstance(result.return_result, float)
    assert "CURRENT ENERGY" in result.extras["qcvars"]


# ---------------------------------------------------------------------------
# Tests: harvester.harvest() — energy parsing
# ---------------------------------------------------------------------------


def _make_water_mol():
    return Molecule.from_data(
        """
        0 1
        O  0.000000000  0.000000000  0.117176320
        H  0.000000000  0.757335953 -0.468705281
        H  0.000000000 -0.757335953 -0.468705281
        units angstrom
        """
    )


# rq-71b27e68
@requires_cclib
def test_harvest_hf_energy():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    expected = -2041.3 * _EV_TO_HARTREE
    assert abs(float(qcvars["HF TOTAL ENERGY"]) - expected) < 1e-8
    assert abs(float(qcvars["CURRENT ENERGY"]) - float(qcvars["HF TOTAL ENERGY"])) < 1e-12


# rq-6bff3228
@requires_cclib
def test_harvest_dft_energy():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2082.1]))
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "b3lyp", "dummy log")
    expected = -2082.1 * _EV_TO_HARTREE
    assert abs(float(qcvars["HF TOTAL ENERGY"]) - expected) < 1e-8
    assert abs(float(qcvars["CURRENT ENERGY"]) - float(qcvars["HF TOTAL ENERGY"])) < 1e-12


# rq-578bb785
@requires_cclib
def test_harvest_mp2_energies():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        mpenergies=np.array([[-2044.8]]),
    )
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "mp2", "dummy log")
    hf_h = -2041.3 * _EV_TO_HARTREE
    mp2_h = -2044.8 * _EV_TO_HARTREE
    assert abs(float(qcvars["HF TOTAL ENERGY"]) - hf_h) < 1e-8
    assert abs(float(qcvars["MP2 TOTAL ENERGY"]) - mp2_h) < 1e-8
    assert abs(float(qcvars["MP2 CORRELATION ENERGY"]) - (mp2_h - hf_h)) < 1e-8
    assert abs(float(qcvars["CURRENT ENERGY"]) - mp2_h) < 1e-8


# rq-1b7ae62e
@requires_cclib
def test_harvest_ccsd_energies():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        ccenergies=np.array([-2045.2]),
    )
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "ccsd", "dummy log")
    hf_h = -2041.3 * _EV_TO_HARTREE
    ccsd_h = -2045.2 * _EV_TO_HARTREE
    assert abs(float(qcvars["CCSD TOTAL ENERGY"]) - ccsd_h) < 1e-8
    assert abs(float(qcvars["CCSD CORRELATION ENERGY"]) - (ccsd_h - hf_h)) < 1e-8
    assert abs(float(qcvars["CURRENT ENERGY"]) - ccsd_h) < 1e-8


# rq-eddef743
@requires_cclib
def test_harvest_ccsd_t_energies():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        ccenergies=np.array([-2045.5]),
    )
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "ccsd(t)", "dummy log")
    ccsdt_h = -2045.5 * _EV_TO_HARTREE
    assert abs(float(qcvars["CCSD(T) TOTAL ENERGY"]) - ccsdt_h) < 1e-8
    assert abs(float(qcvars["CURRENT ENERGY"]) - ccsdt_h) < 1e-8


# ---------------------------------------------------------------------------
# Tests: harvester.harvest() — gradient parsing
# ---------------------------------------------------------------------------


# rq-36736115
@requires_cclib
def test_harvest_gradient_negated():
    in_mol = _make_water_mol()
    forces = np.array([[0.01, 0.0, -0.02], [0.0, 0.005, 0.01], [0.0, -0.005, 0.01]])
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        grads=forces[np.newaxis, :, :],
    )
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    expected = -forces.flatten()
    assert grad is not None
    np.testing.assert_allclose(grad, expected, atol=1e-12)


# rq-3b0adf96
@requires_cclib
def test_harvest_no_gradient_when_absent():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))
    # grads is None (no attribute)
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    assert grad is None


# ---------------------------------------------------------------------------
# Tests: harvester.harvest() — Hessian parsing
# ---------------------------------------------------------------------------


# rq-f2e32162
@requires_cclib
def test_harvest_hessian():
    in_mol = _make_water_mol()
    hess_matrix = np.eye(9) * 0.5
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        hessian=hess_matrix,
    )
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    assert hess is not None
    assert hess.shape == (9, 9)
    np.testing.assert_allclose(hess, hess_matrix)


# rq-344f3432
@requires_cclib
def test_harvest_no_hessian_when_absent():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    assert hess is None


# ---------------------------------------------------------------------------
# Tests: harvester.harvest() — properties parsing
# ---------------------------------------------------------------------------


# rq-4c4a1542
@requires_cclib
def test_harvest_dipole_moment():
    in_mol = _make_water_mol()
    dipole_vec = np.array([0.0, 1.85, 0.0])
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        moments=[None, dipole_vec],
    )
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    assert abs(float(qcvars["DIPOLE MOMENT"]) - 1.85) < 1e-6


# rq-ffae3fa6
@requires_cclib
def test_harvest_mulliken_charges():
    in_mol = _make_water_mol()
    charges = [-0.5, 0.25, 0.25]
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        atomcharges={"mulliken": charges},
    )
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    assert "MULLIKEN CHARGES" in qcvars
    assert len(qcvars["MULLIKEN CHARGES"]) == 3


# rq-53f8ed12
@requires_cclib
def test_harvest_missing_mulliken_no_error():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        atomcharges={},
    )
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    assert "MULLIKEN CHARGES" not in qcvars


# ---------------------------------------------------------------------------
# Tests: harvester.harvest() — geometry and molecule output
# ---------------------------------------------------------------------------


# rq-a76c3d7b
@requires_cclib
def test_harvest_output_molecule_symbols():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    assert list(out_mol.symbols) == ["O", "H", "H"]


# rq-a9194afb
@requires_cclib
def test_harvest_nre_cross_check_passes():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))
    # Should not raise
    with _patch_cclib(ccdata):
        harvest(in_mol, "hf", "dummy log")


# rq-0b4cf843
@requires_cclib
def test_harvest_nre_mismatch_raises():
    in_mol = _make_water_mol()
    # Use very different geometry to force NRE mismatch
    far_coords = np.array([[[0.0, 0.0, 100.0], [0.0, 0.0, 200.0], [0.0, 0.0, 300.0]]])
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        atomcoords=far_coords,
    )
    with _patch_cclib(ccdata):
        with pytest.raises(ValueError, match="inconsistent"):
            harvest(in_mol, "hf", "dummy log")


# ---------------------------------------------------------------------------
# Tests: harvester.harvest() — error conditions
# ---------------------------------------------------------------------------


# rq-b7c7b8cb
@requires_cclib
def test_harvest_no_coords_raises():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        atomcoords=np.array([]),  # empty
    )
    with _patch_cclib(ccdata):
        with pytest.raises(ValueError):
            harvest(in_mol, "hf", "dummy log")


# rq-bb775969
@requires_cclib
def test_harvest_cclib_parse_failure_propagates():
    in_mol = _make_water_mol()
    with patch(
        "cclib.io.ccopen",
        side_effect=Exception("cclib parse error"),
    ):
        with pytest.raises(Exception, match="cclib parse error"):
            harvest(in_mol, "hf", "unparseable log")


# ---------------------------------------------------------------------------
# Tests: harvester.is_normal_termination()
# ---------------------------------------------------------------------------


# rq-90141f40
def test_is_normal_termination_true():
    log = "Some output\nNormal termination of Gaussian 16\n"
    assert is_normal_termination(log) is True


# rq-12c41d5e
def test_is_normal_termination_false():
    log = "Some output\nError termination via Lnk1e\n"
    assert is_normal_termination(log) is False


# ---------------------------------------------------------------------------
# Tests: ghost-atom support
# ---------------------------------------------------------------------------


def _make_hene_ghost_mol():
    """He real at origin, Ne ghost at (2.5, 0, 0) Å, with fix_com/fix_orientation."""
    return Molecule.from_data(
        """
        0 1
        He 0.0 0.0 0.0
        @Ne 2.5 0.0 0.0
        nocom
        noreorient
        """
    )


def _default_taskconfig(**overrides):
    from qcengine.config import TaskConfig

    cfg = dict(
        scratch_directory=None,
        scratch_messy=False,
        ncores=1,
        nnodes=1,
        memory=1.0,
        retries=0,
        mpiexec_command=None,
    )
    cfg.update(overrides)
    return TaskConfig(**cfg)


def _spec(driver="energy", method="HF", basis="STO-3G"):
    return AtomicSpecification(
        program="gaussian",
        driver=driver,
        model={"method": method, "basis": basis},
    )


# rq-2229f087
def test_ghost_real_atom_emits_plain_symbol(harness):
    mol = Molecule.from_data(
        {
            "symbols": ["H"],
            "geometry": [0.0, 0.0, 0.0],
            "real": [True],
            "molecular_charge": 0,
            "molecular_multiplicity": 2,
        }
    )
    inp = AtomicInput(molecule=mol, specification=_spec())
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    com = job["infiles"]["gaussian.com"]
    atom_lines = [ln for ln in com.splitlines() if ln.lstrip().startswith(("H ", "H-Bq"))]
    assert any(ln.lstrip().startswith("H ") for ln in atom_lines)
    assert not any(ln.lstrip().startswith("H-Bq ") for ln in atom_lines)


# rq-916e4919
def test_ghost_atom_emits_bq_suffix(harness):
    mol = Molecule.from_data(
        {
            "symbols": ["H"],
            "geometry": [0.0, 0.0, 0.0],
            "real": [False],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        }
    )
    inp = AtomicInput(molecule=mol, specification=_spec())
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    com = job["infiles"]["gaussian.com"]
    assert any(ln.lstrip().startswith("H-Bq ") for ln in com.splitlines())


# rq-2a015edf
def test_ghost_mixed_real_and_ghost(harness):
    mol = _make_hene_ghost_mol()
    inp = AtomicInput(molecule=mol, specification=_spec())
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    lines = [ln for ln in job["infiles"]["gaussian.com"].splitlines() if ln.strip()]
    he_idx = next(i for i, ln in enumerate(lines) if ln.lstrip().startswith("He "))
    ne_idx = next(i for i, ln in enumerate(lines) if ln.lstrip().startswith("Ne-Bq "))
    assert he_idx < ne_idx  # input order preserved


# rq-a100502c
def test_ghost_route_line_unchanged(harness):
    mol = _make_hene_ghost_mol()
    inp = AtomicInput(molecule=mol, specification=_spec(method="HF", basis="6-31G*"))
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    lines = job["infiles"]["gaussian.com"].splitlines()
    route_lines = [ln for ln in lines if ln.startswith("#P")]
    assert route_lines == ["#P HF/6-31G*"]


# rq-ce76f384
def test_ghost_all_atoms_ghost(harness):
    mol = Molecule.from_data(
        {
            "symbols": ["H", "H"],
            "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 1.4],
            "real": [False, False],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        }
    )
    inp = AtomicInput(molecule=mol, specification=_spec())
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    com = job["infiles"]["gaussian.com"]
    atom_lines = [
        ln for ln in com.splitlines()
        if ln.lstrip().startswith(("H ", "H-Bq"))
    ]
    assert len(atom_lines) == 2
    assert all("-Bq" in ln for ln in atom_lines)


# rq-ceef73bd
def test_ghost_coordinate_precision_preserved(harness):
    mol = Molecule.from_data(
        {
            "symbols": ["H"],
            "geometry": [1.23456789012 / qcel.constants.conversion_factor("bohr", "angstrom"), 0.0, 0.0],
            "real": [False],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        }
    )
    inp = AtomicInput(molecule=mol, specification=_spec())
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    com = job["infiles"]["gaussian.com"]
    ghost_line = next(ln for ln in com.splitlines() if "H-Bq" in ln)
    assert re.search(r"\d+\.\d{10}", ghost_line)


# rq-b27d4aff
@requires_cclib
def test_ghost_out_mol_real_preserved():
    in_mol = _make_hene_ghost_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-77.7]),
        atomcoords=np.array([[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]]]),
        atomnos=np.array([2, 10]),
    )
    with _patch_cclib(ccdata):
        _, _, _, out_mol = harvest(in_mol, "hf", "dummy log")
    assert list(out_mol.real) == [True, False]


# rq-731f01ba
@requires_cclib
def test_ghost_out_mol_symbols_from_in_mol():
    in_mol = _make_hene_ghost_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-77.7]),
        atomcoords=np.array([[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]]]),
        atomnos=np.array([2, 10]),
    )
    with _patch_cclib(ccdata):
        _, _, _, out_mol = harvest(in_mol, "hf", "dummy log")
    assert list(out_mol.symbols) == ["He", "Ne"]


# rq-7411efa0
@requires_cclib
def test_ghost_atom_count_sanity_check_passes():
    in_mol = _make_hene_ghost_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-77.7]),
        atomcoords=np.array([[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]]]),
        atomnos=np.array([2, 10]),
    )
    with _patch_cclib(ccdata):
        # Should not raise: 2 atoms in cclib output == 2 atoms in in_mol
        harvest(in_mol, "hf", "dummy log")


# rq-66067019
@requires_cclib
def test_ghost_atom_count_mismatch_raises():
    in_mol = _make_hene_ghost_mol()  # 2 atoms
    # cclib reports 3 atoms but in_mol has 2 → mismatch
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-77.7]),
        atomcoords=np.array([[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [5.0, 0.0, 0.0]]]),
        atomnos=np.array([2, 10, 10]),
    )
    with _patch_cclib(ccdata):
        with pytest.raises(ValueError, match="atom count"):
            harvest(in_mol, "hf", "dummy log")


# rq-a69d0898
@requires_cclib
def test_ghost_nre_check_excludes_ghosts():
    in_mol = _make_hene_ghost_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-77.7]),
        atomcoords=np.array([[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]]]),
        atomnos=np.array([2, 10]),
    )
    with _patch_cclib(ccdata):
        # NRE = 0 for both (single real atom), no ValueError
        harvest(in_mol, "hf", "dummy log")


# rq-b57c882b
@requires_cclib
def test_ghost_nre_mismatch_still_raises():
    # Use a 2 real + 1 ghost system so the real-atom-only NRE is nonzero
    # and can mismatch when the parsed geometry differs.
    in_mol = Molecule.from_data(
        {
            "symbols": ["H", "H", "He"],
            "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 1.4, 5.0, 0.0, 0.0],
            "real": [True, True, False],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
            "fix_com": True,
            "fix_orientation": True,
        }
    )
    # Parsed geometry has the two real H atoms far apart, so the
    # real-atom-only NRE differs from in_mol's by far more than 1e-3 Ha.
    far_coords = np.array([[[0.0, 0.0, 0.0], [0.0, 0.0, 50.0], [5.0, 0.0, 0.0]]])
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-1.0]),
        atomcoords=far_coords,
        atomnos=np.array([1, 1, 2]),
    )
    with _patch_cclib(ccdata):
        with pytest.raises(ValueError, match="inconsistent"):
            harvest(in_mol, "hf", "dummy log")


# rq-e5a98119
@requires_cclib
def test_ghost_nre_check_unchanged_for_all_real():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))
    with _patch_cclib(ccdata):
        # Existing behaviour: no ghosts, full-mol NRE comparison, no raise
        harvest(in_mol, "hf", "dummy log")


# rq-d97cc860
@requires_cclib
def test_ghost_gradient_spans_real_and_ghost():
    in_mol = _make_hene_ghost_mol()
    # Forces: nonzero on He, zero on Ne (ghost)
    forces = np.array([[0.001, 0.0, 0.0], [0.0, 0.0, 0.0]])
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-77.7]),
        atomcoords=np.array([[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]]]),
        atomnos=np.array([2, 10]),
        grads=forces[np.newaxis, :, :],
    )
    with _patch_cclib(ccdata):
        _, grad, _, _ = harvest(in_mol, "hf", "dummy log")
    assert grad is not None
    assert grad.shape == (6,)
    # Ghost atom (index 1) gradient is zero
    np.testing.assert_allclose(grad[3:6], np.zeros(3), atol=1e-12)


# rq-f8ce079b
@requires_cclib
def test_ghost_align_uses_generic_ghosts():
    # The harvester constructs calc_mol via qcelemental.models.Molecule (v1)
    # and dispatches align() on it, so spy on the same class.
    from qcelemental.models import Molecule as HarvesterMolecule

    in_mol = _make_hene_ghost_mol()
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-77.7]),
        atomcoords=np.array([[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]]]),
        atomnos=np.array([2, 10]),
    )

    captured = {"calls": []}
    real_align = HarvesterMolecule.align

    def spy_align(self, other, **kwargs):
        captured["calls"].append(kwargs)
        return real_align(self, other, **kwargs)

    with _patch_cclib(ccdata):
        with patch.object(HarvesterMolecule, "align", spy_align):
            harvest(in_mol, "hf", "dummy log")

    assert captured["calls"], "expected at least one Molecule.align call"
    assert all(call.get("generic_ghosts") is True for call in captured["calls"])


# rq-04028610
@requires_cclib
def test_ghost_compute_preserves_ghost_designation(harness):
    in_mol = _make_hene_ghost_mol()
    inp = AtomicInput(molecule=in_mol, specification=_spec(method="HF", basis="STO-3G"))

    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-77.7]),
        atomcoords=np.array([[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]]]),
        atomnos=np.array([2, 10]),
    )
    log_text = _NORMAL_LOG
    job_record = {
        "infiles": {"gaussian.com": "%NProcShared=1\n#P HF/STO-3G\n"},
        "command": ["/g16", "gaussian.com"],
        "scratch_directory": None,
        "scratch_messy": False,
    }
    dexe = {"outfiles": {"gaussian.log": log_text}, "stdout": "", "stderr": ""}

    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with patch("qcengine.programs.gaussian.runner.which_import", return_value=True):
            with patch.object(GaussianHarness, "build_input", return_value=job_record):
                with patch.object(GaussianHarness, "execute", return_value=(True, dexe)):
                    with patch.object(GaussianHarness, "get_version", return_value="16.C.02"):
                        with _patch_cclib(ccdata):
                            result = harness.compute(inp, MagicMock())

    assert list(result.molecule.real) == [True, False]
    assert list(result.molecule.symbols) == ["He", "Ne"]


# ---------------------------------------------------------------------------
# Tests: tricky-ghost (mp2(full) method + NoSymm)
# ---------------------------------------------------------------------------


# rq-2cb2c5c0
def test_tricky_ghost_mp2_full_method_string():
    """muster_modelchem("mp2(full)") emits Gaussian's "MP2(Full)" route token."""
    result = muster_modelchem("mp2(full)", "energy", 1)
    assert result["method_string"] == "MP2(Full)"


# rq-2cb2c5c0
def test_tricky_ghost_mp2_full_method_string_gradient():
    """muster_modelchem("mp2(full)", "gradient") uses Force=NoStep job type."""
    result = muster_modelchem("mp2(full)", "gradient", 1)
    assert result["method_string"] == "MP2(Full)"
    assert result["job_type"] == "Force=NoStep"


# rq-2cb2c5c0
def test_tricky_ghost_mp2_full_open_shell():
    """Open-shell mp2(full) gets a U prefix like other MP2 variants."""
    result = muster_modelchem("mp2(full)", "energy", 2)
    assert result["method_string"] == "UMP2(Full)"


# rq-2cb2c5c0
def test_tricky_ghost_nosymm_emitted_bare():
    """NoSymm in user_keywords with empty value emits as a bare token."""
    line = build_route_line("MP2(Full)", "6-31g*", "Force=NoStep", {}, {"NoSymm": ""})
    assert line.startswith("#P MP2(Full)/6-31g* Force=NoStep")
    assert " NoSymm" in line
    assert "NoSymm=" not in line  # bare, not "NoSymm="


# rq-2cb2c5c0
@requires_cclib
def test_tricky_ghost_mp2_full_harvested_as_mp2():
    """harvest() normalises method 'mp2(full)' to 'mp2' so MP2 qcvars are populated."""
    in_mol = _make_water_mol()
    log_with_mp2 = "EUMP2 =    -0.76276030623D+02\n"
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        mpenergies=np.array([[-2044.8]]),
    )
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "mp2(full)", log_with_mp2)
    assert "MP2 TOTAL ENERGY" in qcvars
    assert "MP2 CORRELATION ENERGY" in qcvars
    # CURRENT ENERGY should be the MP2 total (highest level)
    assert abs(float(qcvars["CURRENT ENERGY"]) - float(qcvars["MP2 TOTAL ENERGY"])) < 1e-12


# rq-dc9a34d9
def test_pgline_regex_matches_c2v():
    """The pgline regex extracts a non-* Gaussian point group."""
    line = " Full point group                 C2V     NOp   4"
    m = re.search(r"Full point group\s+(?P<pg>\S+)", line)
    assert m is not None
    assert m.group("pg") == "C2V"


# rq-dc9a34d9
def test_pgline_regex_matches_linear():
    """The pgline regex extracts a linear point group containing '*'."""
    line = " Full point group                 C*V     NOp   4"
    m = re.search(r"Full point group\s+(?P<pg>\S+)", line)
    assert m is not None
    assert m.group("pg") == "C*V"


# ---------------------------------------------------------------------------
# Tests: atom-label tolerance
# ---------------------------------------------------------------------------


def _build_labeled_h_mol(real_flags=None):
    """4 H atoms with the labels from test_atom_labels."""
    geom_bohr = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 5.0, 5.0, 0.0]
    return Molecule.from_data(
        {
            "symbols": ["H", "H", "H", "H"],
            "atom_labels": ["", "5", "_other", "_4sq"],
            "geometry": geom_bohr,
            "real": real_flags if real_flags is not None else [True, True, True, True],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        }
    )


# rq-20063525
def test_atom_labels_single_labeled_h_emits_bare_symbol(harness):
    mol = Molecule.from_data(
        {
            "symbols": ["H"],
            "atom_labels": ["_other"],
            "geometry": [0.0, 0.0, 0.0],
            "real": [True],
            "molecular_charge": 0,
            "molecular_multiplicity": 2,
        }
    )
    inp = AtomicInput(molecule=mol, specification=_spec())
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    com = job["infiles"]["gaussian.com"]
    atom_lines = [
        ln for ln in com.splitlines()
        if ln.lstrip().startswith("H ") or "_other" in ln
    ]
    assert any(ln.lstrip().startswith("H ") for ln in atom_lines)
    assert not any("_other" in ln for ln in atom_lines)


# rq-da7df6e9
def test_atom_labels_four_labeled_h_emit_bare_symbols(harness):
    mol = _build_labeled_h_mol()
    inp = AtomicInput(molecule=mol, specification=_spec(method="MP2", basis="aug-cc-pvdz"))
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    com = job["infiles"]["gaussian.com"]
    h_lines = [ln for ln in com.splitlines() if ln.lstrip().startswith("H ")]
    assert len(h_lines) == 4
    for forbidden in ("H5", "H_other", "H_4sq"):
        assert forbidden not in com


# rq-7a592160
def test_atom_labels_empty_labels_same_output_as_unlabeled(harness):
    geom = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 5.0, 5.0, 0.0]
    mol_with_empty = Molecule.from_data(
        {
            "symbols": ["H", "H", "H", "H"],
            "atom_labels": ["", "", "", ""],
            "geometry": geom,
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        }
    )
    mol_without = Molecule.from_data(
        {
            "symbols": ["H", "H", "H", "H"],
            "geometry": geom,
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        }
    )
    spec = _spec(method="MP2", basis="aug-cc-pvdz")
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        com_with = harness.build_input(AtomicInput(molecule=mol_with_empty, specification=spec), _default_taskconfig())["infiles"]["gaussian.com"]
        com_without = harness.build_input(AtomicInput(molecule=mol_without, specification=spec), _default_taskconfig())["infiles"]["gaussian.com"]
    assert com_with == com_without


# rq-090d144e
def test_atom_labels_route_line_clean(harness):
    mol = _build_labeled_h_mol()
    inp = AtomicInput(molecule=mol, specification=_spec(method="MP2", basis="aug-cc-pvdz"))
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    route = next(ln for ln in job["infiles"]["gaussian.com"].splitlines() if ln.startswith("#P"))
    for forbidden in ("5", "_other", "_4sq"):
        assert forbidden not in route or forbidden in "aug-cc-pvdz"


# rq-c0a7e84e
def test_atom_labels_build_input_does_not_raise(harness):
    mol = _build_labeled_h_mol()
    inp = AtomicInput(molecule=mol, specification=_spec(method="MP2", basis="aug-cc-pvdz"))
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        # Should not raise
        harness.build_input(inp, _default_taskconfig())


# rq-0892bdc0 — ghost atom with non-empty label still emits "<sym>-Bq", label dropped
def test_atom_labels_ghost_with_label_emits_bq(harness):
    mol = Molecule.from_data(
        {
            "symbols": ["He", "H"],
            "atom_labels": ["", "ghost_a"],
            "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 2.0],
            "real": [True, False],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        }
    )
    inp = AtomicInput(molecule=mol, specification=_spec())
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    com = job["infiles"]["gaussian.com"]
    ghost_lines = [ln for ln in com.splitlines() if "H-Bq" in ln]
    assert len(ghost_lines) == 1
    assert "ghost_a" not in com


# rq-d4eaf3c7
@requires_cclib
def test_atom_labels_out_mol_preserves_labels_when_fixed():
    """When fix_com=True and fix_orientation=True, out_mol is in_mol → labels preserved."""
    mol = Molecule.from_data(
        {
            "symbols": ["H", "H", "H", "H"],
            "atom_labels": ["", "5", "_other", "_4sq"],
            "geometry": [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 5.0, 5.0, 0.0],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
            "fix_com": True,
            "fix_orientation": True,
        }
    )
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-45.0]),
        atomcoords=np.array(
            [[[0.0, 0.0, 0.0], [2.645886, 0.0, 0.0], [0.0, 2.645886, 0.0], [2.645886, 2.645886, 0.0]]]
        ),
        atomnos=np.array([1, 1, 1, 1]),
    )
    with _patch_cclib(ccdata):
        _, _, _, out_mol = harvest(mol, "hf", "dummy log")
    assert list(out_mol.atom_labels) == ["", "5", "_other", "_4sq"]


# ---------------------------------------------------------------------------
# Tests: empirical dispersion (germinate)
# ---------------------------------------------------------------------------


# rq-0c1b5066
def test_dispersion_germinate_d3():
    """muster_modelchem strips '-d3' and emits EmpiricalDispersion=GD3."""
    result = muster_modelchem("b3lyp-d3", "energy", 1)
    assert result["method_string"] == "B3LYP"
    assert result["extra_keywords"]["EmpiricalDispersion"] == "GD3"
    assert result["dispersion_level"] == "D3"
    assert result["functional_name"] == "B3LYP"


# rq-0c1b5066
def test_dispersion_germinate_d3bj():
    result = muster_modelchem("b3lyp-d3bj", "energy", 1)
    assert result["method_string"] == "B3LYP"
    assert result["extra_keywords"]["EmpiricalDispersion"] == "GD3BJ"
    assert result["dispersion_level"] == "D3(BJ)"


# rq-0c1b5066
def test_dispersion_germinate_d2():
    result = muster_modelchem("b3lyp-d2", "energy", 1)
    assert result["method_string"] == "B3LYP"
    assert result["extra_keywords"]["EmpiricalDispersion"] == "GD2"
    assert result["dispersion_level"] == "D2"


# rq-0c1b5066
def test_dispersion_germinate_d3zero2b_alias():
    result = muster_modelchem("b3lyp-d3zero2b", "energy", 1)
    assert result["dispersion_level"] == "D3"
    assert result["functional_name"] == "B3LYP"


# rq-0c1b5066
def test_dispersion_germinate_bare_d_alias():
    """The bare '-d' alias maps to D2 (per QCEngine's get_dispersion_aliases)."""
    result = muster_modelchem("b3lyp-d", "energy", 1)
    assert result["dispersion_level"] == "D2"


# rq-0c1b5066
def test_dispersion_germinate_no_suffix_unchanged():
    result = muster_modelchem("b3lyp", "energy", 1)
    assert result["method_string"] == "B3LYP"
    assert "EmpiricalDispersion" not in result["extra_keywords"]
    assert "dispersion_level" not in result


# rq-0c1b5066
def test_dispersion_germinate_hyphenated_non_dispersion_preserved():
    """cam-B3LYP must not get falsely stripped (cam isn't an alias)."""
    result = muster_modelchem("cam-b3lyp", "energy", 1)
    assert result["method_string"] == "CAM-B3LYP"
    assert "EmpiricalDispersion" not in result["extra_keywords"]
    assert "dispersion_level" not in result


# rq-0c1b5066
def test_dispersion_germinate_hf_d3():
    result = muster_modelchem("hf-d3", "energy", 1)
    assert result["method_string"] == "HF"
    assert result["extra_keywords"]["EmpiricalDispersion"] == "GD3"
    assert result["functional_name"] == "HF"


# rq-0c1b5066
def test_dispersion_germinate_open_shell():
    result = muster_modelchem("b3lyp-d3", "energy", 2)
    assert result["method_string"] == "UB3LYP"
    assert result["extra_keywords"]["EmpiricalDispersion"] == "GD3"
    assert result["dispersion_level"] == "D3"
    assert result["functional_name"] == "B3LYP"


# rq-0c1b5066
def test_dispersion_germinate_rohf_no_double_prefix():
    result = muster_modelchem("rohf-d3", "energy", 2)
    assert result["method_string"] == "ROHF"
    assert result["extra_keywords"]["EmpiricalDispersion"] == "GD3"


# rq-0c1b5066
@pytest.mark.parametrize("method", ["b3lyp-d3m", "b3lyp-d3mbj", "b3lyp-d4", "b3lyp-d3op", "b3lyp-nl"])
def test_dispersion_germinate_unsupported_levels_raise(method):
    with pytest.raises(InputError, match=r"not natively supported"):
        muster_modelchem(method, "energy", 1)


# rq-0c1b5066 — explicit "longer suffix wins" test
def test_dispersion_germinate_longer_suffix_wins():
    """d3bj must be matched before d3 + stray 'bj'."""
    result = muster_modelchem("b3lyp-d3bj", "energy", 1)
    assert result["dispersion_level"] == "D3(BJ)"
    assert result["functional_name"] == "B3LYP"


# ---------------------------------------------------------------------------
# Tests: empirical dispersion (runner conflict detection)
# ---------------------------------------------------------------------------


# rq-59b816b7
def test_dispersion_conflict_raises(harness, water):
    inp = AtomicInput(
        molecule=water,
        specification=AtomicSpecification(
            program="gaussian",
            driver="energy",
            model={"method": "b3lyp-d3", "basis": "STO-3G"},
            keywords={"EmpiricalDispersion": "GD3BJ"},
        ),
    )
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with pytest.raises(InputError, match="overspecified"):
            harness.build_input(inp, _default_taskconfig())


# rq-59b816b7
def test_dispersion_conflict_case_insensitive(harness, water):
    inp = AtomicInput(
        molecule=water,
        specification=AtomicSpecification(
            program="gaussian",
            driver="energy",
            model={"method": "b3lyp-d3", "basis": "STO-3G"},
            keywords={"empiricaldispersion": "GD3"},
        ),
    )
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        with pytest.raises(InputError, match="overspecified"):
            harness.build_input(inp, _default_taskconfig())


# rq-59b816b7
def test_dispersion_user_keyword_only_passes_through(harness, water):
    """If the method has no dispersion suffix, a user EmpiricalDispersion kw passes through."""
    inp = AtomicInput(
        molecule=water,
        specification=AtomicSpecification(
            program="gaussian",
            driver="energy",
            model={"method": "B3LYP", "basis": "STO-3G"},
            keywords={"EmpiricalDispersion": "GD3"},
        ),
    )
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    assert "EmpiricalDispersion=GD3" in job["infiles"]["gaussian.com"]


# rq-59b816b7
def test_dispersion_method_suffix_only(harness, water):
    """method='b3lyp-d3' produces exactly one EmpiricalDispersion token and bare B3LYP/basis."""
    inp = AtomicInput(
        molecule=water,
        specification=AtomicSpecification(
            program="gaussian",
            driver="energy",
            model={"method": "b3lyp-d3", "basis": "STO-3G"},
        ),
    )
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/g16"):
        job = harness.build_input(inp, _default_taskconfig())
    com = job["infiles"]["gaussian.com"]
    assert com.count("EmpiricalDispersion=GD3") == 1
    # Route line has B3LYP/STO-3G, NOT B3LYP-D3/STO-3G
    route = next(ln for ln in com.splitlines() if ln.startswith("#P"))
    assert "B3LYP/STO-3G" in route
    assert "B3LYP-D3/" not in route


# ---------------------------------------------------------------------------
# Tests: empirical dispersion (harvester)
# ---------------------------------------------------------------------------


_DISP_LOG_SNIPPET = (
    " SCF Done:  E(RB3LYP) =  -76.4014123456     A.U. after   12 cycles\n"
    " Dispersion energy=       -0.0046617305 Hartrees.\n"
)


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_regex_path():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2079.3]))
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "b3lyp-d3", _DISP_LOG_SNIPPET)
    assert abs(float(qcvars["DISPERSION CORRECTION ENERGY"]) - (-0.0046617305)) < 1e-12


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_cclib_fallback():
    """When the log has no 'Dispersion energy=' line, fall back to ccdata.dispersionenergies."""
    in_mol = _make_water_mol()
    # eV values that round-trip to a known Hartree quantity
    disp_ev = -0.1269192
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2079.3]))
    ccdata.dispersionenergies = np.array([disp_ev])
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "b3lyp-d3", "SCF Done: E(RB3LYP) = -76.4 A.U. after 10 cycles\n")
    expected = disp_ev * _EV_TO_HARTREE
    assert abs(float(qcvars["DISPERSION CORRECTION ENERGY"]) - expected) < 1e-8


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_functional_specific_qcvars():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2079.3]))
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "b3lyp-d3", _DISP_LOG_SNIPPET)
    assert "B3LYP-D3 DISPERSION CORRECTION ENERGY" in qcvars
    assert abs(float(qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"]) - (-0.0046617305)) < 1e-12


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_functional_total_excludes_dispersion():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2079.3]))
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "b3lyp-d3", _DISP_LOG_SNIPPET)
    scf = -76.4014123456
    disp = -0.0046617305
    assert abs(float(qcvars["B3LYP FUNCTIONAL TOTAL ENERGY"]) - scf) < 1e-10
    assert abs(float(qcvars["B3LYP-D3 TOTAL ENERGY"]) - (scf + disp)) < 1e-10
    assert abs(float(qcvars["CURRENT ENERGY"]) - (scf + disp)) < 1e-10


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_d3bj_uses_parenthesised_label():
    in_mol = _make_water_mol()
    log = (
        " SCF Done:  E(RB3LYP) =  -76.4014123456     A.U. after   12 cycles\n"
        " Dispersion energy=       -0.0070000000 Hartrees.\n"
    )
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2079.3]))
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "b3lyp-d3bj", log)
    assert "B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY" in qcvars
    assert "B3LYP-D3(BJ) TOTAL ENERGY" in qcvars


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_scf_hf_unchanged():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2079.3]))
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "b3lyp-d3", _DISP_LOG_SNIPPET)
    scf = -76.4014123456
    # SCF TOTAL ENERGY / HF TOTAL ENERGY excludes the dispersion contribution
    assert abs(float(qcvars["HF TOTAL ENERGY"]) - scf) < 1e-10
    assert abs(float(qcvars["SCF TOTAL ENERGY"]) - scf) < 1e-10


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_method_without_suffix_no_functional_qcvars():
    """If user passes EmpiricalDispersion via keyword (no method suffix), generic qcvars only."""
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2079.3]))
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "b3lyp", _DISP_LOG_SNIPPET)
    assert "DISPERSION CORRECTION ENERGY" in qcvars
    # No formal dash label known → no functional-specific qcvars
    assert "B3LYP-D3 DISPERSION CORRECTION ENERGY" not in qcvars
    assert "B3LYP-D3 TOTAL ENERGY" not in qcvars
    # CURRENT ENERGY still reflects scf + disp for DFT
    scf = -76.4014123456
    disp = -0.0046617305
    assert abs(float(qcvars["CURRENT ENERGY"]) - (scf + disp)) < 1e-10


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_absent_no_qcvars():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2079.3]))
    # No dispersion line in log, no dispersionenergies attribute
    log = " SCF Done:  E(RB3LYP) =  -76.4014123456     A.U. after   12 cycles\n"
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "b3lyp", log)
    assert "DISPERSION CORRECTION ENERGY" not in qcvars
    assert abs(float(qcvars["CURRENT ENERGY"]) - float(qcvars["SCF TOTAL ENERGY"])) < 1e-12


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_nlc_method_no_dispersion_qcvars():
    """wB97M-V is self-consistent NLC; no 'Dispersion energy=' line in log."""
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2079.3]))
    log = " SCF Done:  E(RwB97M-V) =  -76.5     A.U. after   12 cycles\n"
    with _patch_cclib(ccdata):
        qcvars, _, _, _ = harvest(in_mol, "wb97m-v", log)
    assert "DISPERSION CORRECTION ENERGY" not in qcvars


# rq-1d5da737
@requires_cclib
def test_dispersion_harvest_gradient_unchanged():
    """Dispersion in log doesn't disturb the gradient pipeline; cclib's grads is the total."""
    in_mol = _make_water_mol()
    forces = np.array([[0.01, 0.0, -0.02], [0.0, 0.005, 0.01], [0.0, -0.005, 0.01]])
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2079.3]),
        grads=forces[np.newaxis, :, :],
    )
    with _patch_cclib(ccdata):
        qcvars, grad, _, _ = harvest(in_mol, "b3lyp-d3", _DISP_LOG_SNIPPET)
    assert grad is not None
    np.testing.assert_allclose(grad, -forces.flatten(), atol=1e-12)


# ---------------------------------------------------------------------------
# Integration tests (require Gaussian + cclib installed)
# ---------------------------------------------------------------------------


@using("gaussian")
def test_gaussian_energy_water(water):
    """End-to-end HF/STO-3G energy on water via the actual Gaussian executable."""
    import qcengine as qcng

    result = qcng.compute(
        AtomicInput(molecule=water, specification=AtomicSpecification(program="gaussian", driver="energy", model={"method": "HF", "basis": "STO-3G"})),
        "gaussian",
        raise_error=True,
    )
    assert result.success is True
    assert result.return_result < 0.0
    assert "HF TOTAL ENERGY" in result.extras["qcvars"]


@using("gaussian")
def test_gaussian_gradient_water(water):
    """End-to-end HF/STO-3G gradient on water via the actual Gaussian executable."""
    import qcengine as qcng

    result = qcng.compute(
        AtomicInput(molecule=water, specification=AtomicSpecification(program="gaussian", driver="gradient", model={"method": "HF", "basis": "STO-3G"})),
        "gaussian",
        raise_error=True,
    )
    assert result.success is True
    assert np.array(result.return_result).shape == (3, 3)


# rq-17959824
@using("gaussian")
def test_gaussian_b3lyp_d3_water(water):
    """B3LYP-D3/6-31G* energy on water; verify dispersion qcvars are populated and consistent."""
    import qcengine as qcng

    result = qcng.compute(
        AtomicInput(
            molecule=water,
            specification=AtomicSpecification(
                program="gaussian",
                driver="energy",
                model={"method": "b3lyp-d3", "basis": "6-31g*"},
            ),
        ),
        "gaussian",
        raise_error=True,
    )
    assert result.success is True

    qcvars = result.extras["qcvars"]
    # All three Psi4-convention qcvars must be present.
    assert "B3LYP-D3 DISPERSION CORRECTION ENERGY" in qcvars
    assert "B3LYP-D3 TOTAL ENERGY" in qcvars
    assert "B3LYP FUNCTIONAL TOTAL ENERGY" in qcvars
    # Empirical reference values (Gaussian 16, Revision C.02). Tolerance is loose
    # enough to absorb minor patch-level variation in the D3 parameter set.
    atol = 1.0e-6
    ref_func = -76.4087263383
    ref_disp = -0.0000078919
    ref_total = ref_func + ref_disp
    assert abs(float(qcvars["B3LYP FUNCTIONAL TOTAL ENERGY"]) - ref_func) < atol
    assert abs(float(qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"]) - ref_disp) < atol
    assert abs(float(qcvars["B3LYP-D3 TOTAL ENERGY"]) - ref_total) < atol
    # SCF + dispersion == total (numerical consistency).
    scf = float(qcvars["B3LYP FUNCTIONAL TOTAL ENERGY"])
    disp = float(qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"])
    tot = float(qcvars["B3LYP-D3 TOTAL ENERGY"])
    assert abs((scf + disp) - tot) < 1.0e-10
    # CURRENT ENERGY is the dispersion-corrected total.
    assert abs(float(qcvars["CURRENT ENERGY"]) - tot) < 1.0e-10
    # SCF TOTAL ENERGY / HF TOTAL ENERGY exclude dispersion.
    assert abs(float(qcvars["SCF TOTAL ENERGY"]) - scf) < 1.0e-10
