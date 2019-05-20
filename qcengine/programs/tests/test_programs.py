"""
Tests the DQM compute dispatch module
"""

import copy

import pytest

import qcengine as qcng
from qcelemental.models import Molecule, ResultInput
from qcengine import testing

_base_json = {"schema_name": "qcschema_input", "schema_version": 1}


def test_missing_key():
    ret = qcng.compute({"hello": "hi"}, "bleh")
    assert ret.success is False
    assert "hello" in ret.input_data


def test_missing_key_raises():
    with pytest.raises(qcng.exceptions.InputError):
        ret = qcng.compute({"hello": "hi"}, "bleh", raise_error=True)


@testing.using_psi4
def test_psi4_task():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = qcng.get_molecule("water")
    json_data["driver"] = "energy"
    json_data["model"] = {"method": "SCF", "basis": "sto-3g"}
    json_data["keywords"] = {"scf_type": "df"}

    ret = qcng.compute(json_data, "psi4", raise_error=True, return_dict=True)

    assert ret["driver"] == "energy"
    assert "provenance" in ret
    assert "Final Energy" in ret["stdout"]

    for key in ["cpu", "hostname", "username", "wall_time"]:
        assert key in ret["provenance"]

    assert ret["success"] is True


@testing.using_psi4
def test_psi4_internal_failure():

    mol = Molecule.from_data("""0 3
     O    0.000000000000     0.000000000000    -0.068516245955
    """)
    psi4_task = {
        "molecule": mol,
        "driver": "energy",
        "model": {
            "method": "ccsd",
            "basis": "6-31g"
        },
        "keywords": {
            "reference": "rhf"
        }
    }
    with pytest.raises(qcng.exceptions.InputError) as exc:
        ret = qcng.compute(psi4_task, "psi4", raise_error=True)

    assert "reference is only" in str(exc.value)


@testing.using_psi4
def test_psi4_ref_switch():
    inp = ResultInput(
        **{
            "molecule": {
                "symbols": ["Li"],
                "geometry": [0, 0, 0],
                "molecular_multiplicity": 2
            },
            "driver": "energy",
            "model": {
                "method": "B3LYP",
                "basis": "sto-3g"
            },
            "keywords": {
                "scf_type": "df"
            }
        })

    ret = qcng.compute(inp, "psi4", raise_error=True, return_dict=False)

    assert ret.success is True
    assert ret.properties.calcinfo_nalpha == 2
    assert ret.properties.calcinfo_nbeta == 1


@testing.using_rdkit
def test_rdkit_task():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = qcng.get_molecule("water")
    json_data["driver"] = "gradient"
    json_data["model"] = {"method": "UFF"}
    json_data["keywords"] = {}

    ret = qcng.compute(json_data, "rdkit", raise_error=True)

    assert ret.success is True


@testing.using_rdkit
def test_rdkit_connectivity_error():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = qcng.get_molecule("water").dict()
    json_data["driver"] = "gradient"
    json_data["model"] = {"method": "UFF", "basis": ""}
    json_data["keywords"] = {}
    del json_data["molecule"]["connectivity"]

    ret = qcng.compute(json_data, "rdkit")
    assert ret.success is False
    assert "connectivity" in ret.error.error_message

    with pytest.raises(qcng.exceptions.InputError):
        qcng.compute(json_data, "rdkit", raise_error=True)


@testing.using_torchani
def test_torchani_task():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = qcng.get_molecule("water")
    json_data["driver"] = "gradient"
    json_data["model"] = {"method": "ANI1x", "basis": None}
    json_data["keywords"] = {}

    ret = qcng.compute(json_data, "torchani", raise_error=True)

    assert ret.success is True
    assert ret.driver == "gradient"
