"""Tests for NWChem functionality"""
import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using


@pytest.fixture
def nh2():
    smol = """
 # R=1.008 #A=105.0
 0 2
 N   0.000000000000000   0.000000000000000  -0.145912918634892
 H   0.000000000000000  -1.511214298139000   1.013682596946108
 H   0.000000000000000   1.511214298139000   1.013682596946108
 units au
 symmetry c1
"""
    return qcel.models.Molecule.from_data(smol)


@using("nwchem")
def test_b3lyp(nh2):
    # Run NH2
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "b3lyp", "basis": "3-21g"}}
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)

    # Make sure the calculation completed successfully
    assert compare_values(-55.554037, res["return_result"], atol=1e-3)
    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # Check the other status information
    assert res["extras"]["qcvars"]["N ALPHA ELECTRONS"] == "5"
    assert res["extras"]["qcvars"]["N ATOMS"] == "3"
    assert res["extras"]["qcvars"]["N BASIS FUNCTIONS"] == "13"

    # Make sure the properties parsed correctly
    assert compare_values(-55.554037, res["properties"]["return_energy"], atol=1e-3)
    assert res["properties"]["calcinfo_natom"] == 3
    assert res["properties"]["calcinfo_nalpha"] == 5
    assert res["properties"]["calcinfo_nbasis"] == 13


@using("nwchem")
def test_hess(nh2):
    resi = {"molecule": nh2, "driver": "hessian", "model": {"method": "b3lyp", "basis": "3-21g"}}
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=False)
    assert compare_values(-3.5980754370e-02, res.return_result[0, 0], atol=1e-3)
    assert compare_values(0, res.return_result[1, 0], atol=1e-3)
    assert compare_values(0.018208307756, res.return_result[3, 0], atol=1e-3)
    assert np.allclose(res.return_result, res.return_result.T, atol=1e-8)  # Should be symmetric about diagonal

    # Test that the Hessian changes with rotation, but that its determinants remain the same
    shifted_nh2, _ = nh2.scramble(do_shift=False, do_mirror=False, do_rotate=True, do_resort=False)

    resi["molecule"] = shifted_nh2
    res_shifted = qcng.compute(resi, "nwchem", raise_error=True, return_dict=False)
    assert not np.allclose(res.return_result, res_shifted.return_result, atol=1e-8)
    assert np.isclose(np.linalg.det(res.return_result), np.linalg.det(res_shifted.return_result))


@using("nwchem")
def test_gradient(nh2):
    resi = {
        "molecule": nh2,
        "driver": "gradient",
        "model": {"method": "b3lyp", "basis": "3-21g"},
        "keywords": {"dft__convergence__gradient": "1e-6"},
    }
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)
    assert compare_values(4.22418267e-2, res["return_result"][2], atol=1e-7)  # Beyond accuracy of NWChem stdout

    # Rotate the molecule and verify that the gradient changes
    shifted_nh2, _ = nh2.scramble(do_shift=False, do_mirror=False, do_rotate=True, do_resort=False)

    resi["molecule"] = shifted_nh2
    res_shifted = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)

    assert not compare_values(4.22418267e-2, res_shifted["return_result"][2], atol=1e-7)

    # Make sure the two matrices still have the same determinant and norms, as they are just rotations of each other
    #  I am leveraging the fact that the gradients are square, 3x3 matrices just by happenstance of the
    #  test molecule having 3 atoms
    orig_grads = np.reshape(res["return_result"], (-1, 3))
    shif_grads = np.reshape(res_shifted["return_result"], (-1, 3))

    # Test that the magnitude of forces are the same
    assert np.allclose(np.linalg.norm(orig_grads, ord=2, axis=1), np.linalg.norm(shif_grads, ord=2, axis=1))

    # Test that the determinants are the same
    orig_det = np.linalg.det(orig_grads)
    shif_det = np.linalg.det(shif_grads)

    assert np.allclose(orig_det, shif_det)


@pytest.fixture
def h20():
    water = """
-1 2
O 0 0 0
H 0 0 1
H 0 1 0
    """
    return qcel.models.Molecule.from_data(water)


@using("nwchem")
def test_dipole(h20):
    # Run NH2
    resi = {
        "molecule": h20,
        "driver": "properties",
        "model": {"method": "dft", "basis": "3-21g"},
        "keywords": {"dft__xc": "b3lyp", "property__dipole": True},
    }
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)

    # Make sure the calculation completed successfully
    assert compare_values(-75.764944, res["return_result"], atol=1e-3)
    assert res["driver"] == "properties"
    assert "provenance" in res
    assert res["success"] is True

    # Check the other status information
    assert res["extras"]["qcvars"]["N ALPHA ELECTRONS"] == "6"
    assert res["extras"]["qcvars"]["N ATOMS"] == "3"
    assert res["extras"]["qcvars"]["N BASIS FUNCTIONS"] == "13"

    # Make sure the properties parsed correctly
    assert compare_values(-75.764944, res["properties"]["return_energy"], atol=1e-3)
    assert res["properties"]["calcinfo_natom"] == 3
    assert res["properties"]["calcinfo_nalpha"] == 6
    assert res["properties"]["calcinfo_nbasis"] == 13
    # Make sure Dipole Moment and center of charge parsed correctly
    assert compare_values(0.272949872, float(res["extras"]["qcvars"]["TOTAL DIPOLE MOMENT"]), atol=1e-5)
    assert compare_values(-0.00, float(res["extras"]["qcvars"]["DIPOLE MOMENT"][0]), atol=1e-3)
    assert compare_values(-0.00, float(res["extras"]["qcvars"]["DIPOLE MOMENT"][1]), atol=1e-3)
    assert compare_values(-0.272949872, float(res["extras"]["qcvars"]["DIPOLE MOMENT"][2]), atol=1e-5)


@pytest.fixture
def h20v2():
    water = """
O 0 0 0
H 0 0 1
H 0 1 0
    """
    return qcel.models.Molecule.from_data(water)


@using("nwchem")
def test_homo_lumo(h20v2):
    # Run NH2
    resi = {
        "molecule": h20v2,
        "driver": "energy",
        "model": {"method": "dft", "basis": "3-21g"},
        "keywords": {"dft__xc": "b3lyp"},
    }
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)

    # Make sure the calculation completed successfully
    assert compare_values(-75.968095, res["return_result"], atol=1e-3)
    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # Check the other status information
    assert res["extras"]["qcvars"]["N ALPHA ELECTRONS"] == "5"
    assert res["extras"]["qcvars"]["N ATOMS"] == "3"
    assert res["extras"]["qcvars"]["N BASIS FUNCTIONS"] == "13"

    # Make sure the properties parsed correctly
    assert compare_values(-75.968095, res["properties"]["return_energy"], atol=1e-3)
    assert res["properties"]["calcinfo_natom"] == 3
    assert res["properties"]["calcinfo_nalpha"] == 5
    assert res["properties"]["calcinfo_nbasis"] == 13
    # Make sure Dipole Moment and center of charge parsed correctly
    assert compare_values(-0.2636515, float(res["extras"]["qcvars"]["HOMO"][0]), atol=1e-5)
    assert compare_values(0.08207131, float(res["extras"]["qcvars"]["LUMO"][0]), atol=1e-5)
