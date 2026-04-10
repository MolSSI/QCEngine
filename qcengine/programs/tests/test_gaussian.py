"""Tests for the Gaussian harness: germinate, keywords, harvester, and runner."""

import re
from decimal import Decimal
from typing import Any, Dict
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.models import AtomicInput, BasisSet, Molecule

from qcengine.exceptions import InputError, ResourceError, UnknownError
from qcengine.programs.gaussian.germinate import muster_modelchem
from qcengine.programs.gaussian.harvester import harvest, is_normal_termination
from qcengine.programs.gaussian.keywords import build_com_file, build_route_line
from qcengine.programs.gaussian import runner as gaussian_runner
from qcengine.programs.gaussian.runner import GaussianHarness
from qcengine.testing import using

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
    return AtomicInput(molecule=water, driver="energy", model={"method": "HF", "basis": "STO-3G"})


@pytest.fixture
def water_energy_input_all_files(water):
    """AtomicInput with native_files protocol set to 'all' for testing file output."""
    return AtomicInput(
        molecule=water,
        driver="energy",
        model={"method": "HF", "basis": "STO-3G"},
        protocols={"native_files": "all"},
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

    cfg = TaskConfig(scratch_directory=None, scratch_messy=False, ncores=1, memory=0.5, retries=0)
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/usr/local/g16/g16"):
        job = harness.build_input(water_energy_input, cfg)
    assert "%Mem=1GB" in job["infiles"]["gaussian.com"]


# rq-4ad16d95
def test_build_input_memory_whole_gb(water_energy_input, harness):
    from qcengine.config import TaskConfig

    cfg = TaskConfig(scratch_directory=None, scratch_messy=False, ncores=4, memory=8.0, retries=0)
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/usr/local/g16/g16"):
        job = harness.build_input(water_energy_input, cfg)
    assert "%Mem=8GB" in job["infiles"]["gaussian.com"]


# rq-5053204b
def test_build_input_ghost_atoms_raise(harness):
    from qcengine.config import TaskConfig

    ghost_mol = Molecule.from_data(
        {"symbols": ["O", "H", "H"], "geometry": [0, 0, 0.2, 0, 1.4, -0.9, 0, -1.4, -0.9],
         "real": [True, False, True]}
    )
    inp = AtomicInput(molecule=ghost_mol, driver="energy", model={"method": "HF", "basis": "STO-3G"})
    cfg = TaskConfig(scratch_directory=None, scratch_messy=False, ncores=1, memory=1.0, retries=0)
    with patch("qcengine.programs.gaussian.runner._find_gaussian", return_value="/usr/local/g16/g16"):
        with pytest.raises(InputError):
            harness.build_input(inp, cfg)


# rq-0aa99076
def test_build_input_basisset_object_raises(water, harness):
    from qcengine.config import TaskConfig

    bs = BasisSet(name="sto-3g", center_data={}, atom_map=[])
    inp = AtomicInput(molecule=water, driver="energy", model={"method": "HF", "basis": bs})
    cfg = TaskConfig(scratch_directory=None, scratch_messy=False, ncores=1, memory=1.0, retries=0)
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
    inp = AtomicInput(molecule=water, driver="energy", model={"method": "HF", "basis": bs})
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
def test_compute_returns_atomic_result_energy(harness, water_energy_input_all_files, water):
    log_text = _NORMAL_LOG
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))

    with _full_compute_patches(harness, water_energy_input_all_files, log_text, ccdata, None):
        result = harness.compute(water_energy_input_all_files, MagicMock())

    assert result.success is True
    assert isinstance(result.return_result, float)
    assert "gaussian.com" in result.native_files
    assert "gaussian.log" in result.native_files
    assert result.provenance.creator == "Gaussian"


# rq-e62578e2
def test_compute_returns_atomic_result_gradient(harness, water):
    forces = np.array([[0.01, 0.0, -0.02], [0.0, 0.005, 0.01], [0.0, -0.005, 0.01]])
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        grads=forces[np.newaxis, :, :],  # shape (1, 3, 3)
    )
    inp = AtomicInput(molecule=water, driver="gradient", model={"method": "HF", "basis": "STO-3G"})
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
def test_compute_returns_atomic_result_hessian(harness, water):
    hess_matrix = np.eye(9) * 0.5
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        hessian=hess_matrix,
    )
    inp = AtomicInput(molecule=water, driver="hessian", model={"method": "HF", "basis": "STO-3G"})
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
def test_compute_returns_atomic_result_properties(harness, water):
    ccdata = _make_ccdata_with_attrs(
        scfenergies=np.array([-2041.3]),
        moments=[None, np.array([0.0, 1.85, 0.0])],
        atomcharges={"mulliken": [-0.5, 0.25, 0.25]},
    )
    inp = AtomicInput(molecule=water, driver="properties", model={"method": "HF", "basis": "STO-3G"})
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
def test_harvest_hf_energy():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    expected = -2041.3 * _EV_TO_HARTREE
    assert abs(float(qcvars["HF TOTAL ENERGY"]) - expected) < 1e-8
    assert abs(float(qcvars["CURRENT ENERGY"]) - float(qcvars["HF TOTAL ENERGY"])) < 1e-12


# rq-6bff3228
def test_harvest_dft_energy():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2082.1]))
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "b3lyp", "dummy log")
    expected = -2082.1 * _EV_TO_HARTREE
    assert abs(float(qcvars["HF TOTAL ENERGY"]) - expected) < 1e-8
    assert abs(float(qcvars["CURRENT ENERGY"]) - float(qcvars["HF TOTAL ENERGY"])) < 1e-12


# rq-578bb785
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
def test_harvest_output_molecule_symbols():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))
    with _patch_cclib(ccdata):
        qcvars, grad, hess, out_mol = harvest(in_mol, "hf", "dummy log")
    assert list(out_mol.symbols) == ["O", "H", "H"]


# rq-a9194afb
def test_harvest_nre_cross_check_passes():
    in_mol = _make_water_mol()
    ccdata = _make_ccdata_with_attrs(scfenergies=np.array([-2041.3]))
    # Should not raise
    with _patch_cclib(ccdata):
        harvest(in_mol, "hf", "dummy log")


# rq-0b4cf843
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
# Integration tests (require Gaussian + cclib installed)
# ---------------------------------------------------------------------------


@using("gaussian")
def test_gaussian_energy_water(water):
    """End-to-end HF/STO-3G energy on water via the actual Gaussian executable."""
    import qcengine as qcng

    result = qcng.compute(
        {
            "molecule": water.dict(),
            "driver": "energy",
            "model": {"method": "HF", "basis": "STO-3G"},
        },
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
        {
            "molecule": water.dict(),
            "driver": "gradient",
            "model": {"method": "HF", "basis": "STO-3G"},
        },
        "gaussian",
        raise_error=True,
    )
    assert result.success is True
    assert np.array(result.return_result).shape == (3, 3)
