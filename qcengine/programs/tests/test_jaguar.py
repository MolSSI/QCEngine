from types import SimpleNamespace

import numpy as np
import pytest
from qcelemental import constants
from qcelemental.models.v2 import AtomicInput
from qcelemental.util import safe_version

import qcengine as qcng
from qcengine.exceptions import InputError
from qcengine.programs.jaguar import JaguarHarness
from qcengine.testing import uusing


def _input(driver="energy", method="hf", basis="sto-3g", keywords=None, extras=None, charge=0, multiplicity=1):
    return AtomicInput(
        molecule={
            "symbols": ["H", "H"],
            "geometry": [0.0, 0.0, -0.7, 0.0, 0.0, 0.7],
            "molecular_charge": charge,
            "molecular_multiplicity": multiplicity,
        },
        specification={
            "driver": driver,
            "model": {"method": method, "basis": basis},
            "keywords": keywords or {},
            "extras": extras or {},
            "protocols": {"native_files": "all"},
        },
    )


def test_registered():
    assert "jaguar" in qcng.list_all_programs()
    assert isinstance(qcng.get_program("jaguar", check=False), JaguarHarness)


@uusing("jaguar")
def test_versions():
    from schrodinger.application.jaguar.constants import VERSION
    from schrodinger.infra import mm

    encoded_version = int(mm.mmfile_get_product_version("jaguar"))
    major_minor, release = divmod(encoded_version, 1000)
    major, minor = divmod(major_minor, 10)
    harness = qcng.get_program("jaguar")

    assert harness.get_version() == safe_version(f"{major}.{minor}.{release:03d}")
    assert harness.get_suite_version() == safe_version(VERSION)


def test_property_mapping():
    dipole = SimpleNamespace(x=1.0, y=2.0, z=3.0)
    results = SimpleNamespace(
        energy=-1.1,
        nuclear_repulsion=0.7,
        scf_energy=-1.0,
        energy_one_electron=-2.0,
        energy_two_electron=0.3,
        rimp2_ss_energy=-0.01,
        rimp2_os_energy=-0.04,
        rimp2_corr_energy=-0.05,
        rimp2_energy=-1.05,
        dipole_qm=dipole,
    )
    output = SimpleNamespace(
        last_results=results,
        nbasis=5,
        num_occ_orbs_alpha=None,
        num_occ_orbs_beta=None,
        num_occ_orbs=1,
    )

    properties = JaguarHarness._optional_properties(output, _input())

    assert properties["calcinfo_natom"] == 2
    assert properties["calcinfo_nbasis"] == 5
    assert properties["calcinfo_nmo"] == 5
    assert properties["calcinfo_nalpha"] == properties["calcinfo_nbeta"] == 1
    assert properties["return_energy"] == -1.1
    assert properties["scf_total_energy"] == -1.0
    assert properties["mp2_total_energy"] == -1.05
    assert np.allclose(
        properties["scf_dipole_moment"],
        np.array([1.0, 2.0, 3.0]) / constants.dipmom_au2debye,
    )


def test_rohf_electron_counts():
    output = SimpleNamespace(
        nelectron=7,
        nbasis=5,
        num_occ_orbs_alpha=4,
        num_occ_orbs_beta=4,
        last_results=SimpleNamespace(energy=-1.1),
    )
    input_model = _input(charge=1, multiplicity=2)

    properties = JaguarHarness._optional_properties(output, input_model)

    assert properties["calcinfo_nalpha"] == 4
    assert properties["calcinfo_nbeta"] == 3


@pytest.mark.parametrize(
    "driver, derivative_label",
    [("energy", None), ("gradient", "GRADIENT"), ("hessian", "HESSIAN")],
)
def test_hf_qcvars(driver, derivative_label):
    properties = {"return_energy": -1.0, "scf_total_energy": -1.0}
    if derivative_label:
        properties[f"return_{derivative_label.lower()}"] = np.zeros((2, 3))

    qcvars = JaguarHarness._qcvars(properties, driver, "hf")

    assert qcvars["HF TOTAL ENERGY"] == -1.0
    assert qcvars["SCF TOTAL ENERGY"] == -1.0
    assert qcvars["CURRENT REFERENCE ENERGY"] == -1.0
    assert qcvars["CURRENT ENERGY"] == -1.0
    assert "CURRENT CORRELATION ENERGY" not in qcvars
    if derivative_label:
        assert f"HF TOTAL {derivative_label}" in qcvars
        assert f"CURRENT {derivative_label}" in qcvars


def test_rejects_unsupported_driver():
    with pytest.raises(InputError, match="not implemented"):
        JaguarHarness._validate_input(_input(driver="properties"))


def test_rejects_invalid_guess_input():
    input_model = _input(extras={"jaguar": {"guess_input": 42}})
    with pytest.raises(InputError, match="guess_input must be a non-empty string"):
        JaguarHarness._validate_input(input_model)


@uusing("jaguar")
def test_input_mapping(tmp_path):
    input_model = _input(
        driver="gradient",
        method="b3lyp-d3",
        basis="6-31g**",
        keywords={"maxit": 99, "dftname": "bad", "basis": "bad", "igeopt": 2, "isymm": 8},
    )

    jaguar_input = JaguarHarness._build_jaguar_input(input_model, str(tmp_path / "dispatch"))

    assert jaguar_input.getValue("dftname").lower() == "b3lyp-d3"
    assert jaguar_input.getValue("basis").lower() == "6-31g**"
    assert jaguar_input.getValue("maxit") == 99
    assert jaguar_input.getValue("igeopt") == -1
    assert jaguar_input.getValue("ifreq") == 0
    assert jaguar_input.getValue("isymm") == 0
    assert jaguar_input.getValue("molchg") == 0
    assert jaguar_input.getValue("multip") == 1
    assert np.allclose(
        [[atom.x, atom.y, atom.z] for atom in jaguar_input.getStructure().atom],
        np.asarray(input_model.molecule.geometry) * constants.bohr2angstroms,
    )


@uusing("jaguar")
def test_guess_input_mapping(tmp_path):
    guess_text = """&gen
maxit=7
basis=6-31g
&
&zmat
H1 0.0 0.0 -1.0
H2 0.0 0.0  1.0
&
&guess basgss=sto-3g numd=1
    1 Orbital Energy -0.500000 Occupation 1.000000 Symmetry A
  1.000000 0.000000
&
"""
    input_model = _input(
        driver="gradient",
        keywords={"maxit": 99},
        extras={"jaguar": {"guess_input": guess_text}},
    )
    job_base = tmp_path / "guess"

    JaguarHarness._validate_input(input_model)
    jaguar_input = JaguarHarness._build_jaguar_input(input_model, str(job_base))

    assert jaguar_input.getValue("maxit") == 99
    assert jaguar_input.getValue("igeopt") == -1
    assert jaguar_input.sectionDefined("guess")
    assert "Orbital Energy -0.500000" in jaguar_input.getSectionText("guess")

    jaguar_input.save()
    output_text = job_base.with_suffix(".in").read_text()
    assert "&guess" in output_text
    assert "maxit=99" in output_text.replace(" ", "")
    assert "maxit=7" not in output_text.replace(" ", "")


@uusing("jaguar")
def test_ghost_atom_mapping(tmp_path):
    input_model = AtomicInput(
        molecule={
            "symbols": ["H", "He"],
            "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 3.0],
            "real": [True, False],
            "molecular_charge": 0,
            "molecular_multiplicity": 2,
        },
        specification={
            "driver": "energy",
            "model": {"method": "hf", "basis": "sto-3g"},
        },
    )
    job_base = tmp_path / "ghost"

    JaguarHarness._validate_input(input_model)
    jaguar_input = JaguarHarness._build_jaguar_input(input_model, str(job_base))
    jaguar_atoms = list(jaguar_input.getStructure().atom)

    assert [atom.atomic_number for atom in jaguar_atoms] == [1, 2]

    jaguar_input.save()
    input_text = job_base.with_suffix(".in").read_text()
    assert "H1" in input_text
    assert "He2@" in input_text


@uusing("jaguar")
@pytest.mark.parametrize(
    "driver, shape",
    [
        ("energy", ()),
        ("gradient", (2, 3)),
        ("hessian", (6, 6)),
    ],
)
def test_compute(driver, shape, tmp_path):
    result = qcng.compute(
        _input(driver=driver),
        "jaguar",
        raise_error=True,
        task_config={"ncores": 1, "scratch_directory": str(tmp_path)},
        return_version=2,
    )

    harness = qcng.get_program("jaguar")
    assert result.success
    assert np.asarray(result.return_result).shape == shape
    assert result.properties.return_energy is not None
    assert result.provenance.version == harness.get_version()
    assert set(result.extras["jaguar"]) == {"suite_version", "point_group"}
    assert result.extras["jaguar"]["suite_version"] == harness.get_suite_version()
    assert "Jaguar version" in result.stdout
    assert "&gen" in result.native_files["input"]
    generated_input = result.native_files["dispatch.01.in"]
    assert "&gen" in generated_input
    if driver == "hessian":
        assert result.properties.return_hessian.shape == shape
        assert "&hess" in generated_input
