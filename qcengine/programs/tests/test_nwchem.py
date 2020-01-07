"""Tests for NWChem functionality"""
from typing import Tuple

import pytest
import logging
import qcelemental as qcel
from qcelemental.testing import compare_values
from qcelemental.molparse import from_string, to_schema
from qcelemental.models import Molecule

import numpy as np
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


@pytest.fixture(
    params=[
        (
            "sym_geom_project",
            """5
problem molecule
C                    -0.023882416957     2.051854680480     0.015170845154
H                     0.101848034591    -0.056139192140    -0.168622795063
H                     1.991847353162     2.688758815555     0.170465320656
H                    -1.207660826274     2.655247931936    -1.636557172888
H                    -0.982047704522     2.919570674169     1.695146072141""",
        ),
        (
            "sym_map",
            """4
problem_molecule
N                     0.175125259758     1.766334351756    -0.364093068829
H                    -0.073384441932    -0.101051416417     0.167104340622
H                     1.805728057807     2.718491719551     0.155935622654
H                    -1.203844475633     2.681779035110    -1.412361704447""",
        ),
    ]
)
def problematic_mols(request) -> Tuple[Molecule, str]:
    """Convert a molecule to a Molecule object, assuming its units are Bohr"""
    mol_dict = from_string(request.param[1], dtype="xyz")["qm"]
    mol_dict["units"] = "Bohr"
    return to_schema(mol_dict, dtype=2), request.param[0]


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
    assert res["extras"]["qcvars"]["N BASIS"] == "13"

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


@using("nwchem")
def test_autosym_recovery(problematic_mols, caplog):
    molecule, error = problematic_mols
    resi = {
        "molecule": Molecule.from_data(molecule),
        "driver": "gradient",
        "model": {"method": "b3lyp", "basis": "3-21g"},
    }

    with caplog.at_level(logging.DEBUG):
        result = qcng.compute(resi, "nwchem", raise_error=True, return_dict=False)

    assert error in result.extras["observed_errors"]
