"""Tests for NWChem functionality"""
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
    assert res["extras"]["qcvars"]["N BASIS"] == "13"

    # Make sure the properties parsed correctly
    assert compare_values(-55.554037, res["properties"]["return_energy"], atol=1e-3)
    assert res["properties"]["calcinfo_natom"] == 3
    assert res["properties"]["calcinfo_nalpha"] == 5
    assert res["properties"]["calcinfo_nbasis"] == 13


@using("nwchem")
def test_hess(nh2):
    resi = {"molecule": nh2, "driver": "hessian", "model": {"method": "b3lyp", "basis": "3-21g"}}
    res = qcng.compute(resi, "nwchem", raise_error=True, return_dict=True)
    assert compare_values(-3.5980754370e-02, res["return_result"][0], atol=1e-3)


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
