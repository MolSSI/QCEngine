"""
Tests the DQM compute dispatch module
"""

import copy

import numpy as np
import pytest

import qcengine as qcng
from qcelemental.models import Molecule, ResultInput
from qcengine import testing
from qcengine.testing import failure_engine

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
    assert "retries" not in ret["provenance"]


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


@testing.using_mopac
def test_mopac_task():
    json_data = copy.deepcopy(_base_json)
    json_data["molecule"] = qcng.get_molecule("water")
    json_data["driver"] = "gradient"
    json_data["model"] = {"method": "PM6", "basis": None}
    json_data["keywords"] = {}

    ret = qcng.compute(json_data, "mopac", raise_error=True)
    assert ret.extras.keys() >= {"heat_of_formation", "energy_electronic", "dip_vec"}
    energy = pytest.approx(-0.08474117913025125, rel=1.e-5)

    # Check gradient
    ret = qcng.compute(json_data, "mopac", raise_error=True)
    assert ret.extras.keys() >= {"heat_of_formation", "energy_electronic", "dip_vec"}
    assert np.linalg.norm(ret.return_result) == pytest.approx(0.03543560156912385, rel=1.e-4)
    assert ret.properties.return_energy == energy

    # Check energy
    json_data["driver"] = "energy"
    ret = qcng.compute(json_data, "mopac", raise_error=True)
    assert ret.return_result == energy
    assert "== MOPAC DONE ==" in ret.stdout


def test_random_failure_no_retries(failure_engine):

    failure_engine.iter_modes = ["input_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False)
    assert ret.error.error_type == "input_error"
    assert "retries" not in ret.input_data["provenance"].keys()

    failure_engine.iter_modes = ["random_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False)
    assert ret.error.error_type == "random_error"
    assert "retries" not in ret.input_data["provenance"].keys()


def test_random_failure_with_retries(failure_engine):

    failure_engine.iter_modes = ["random_error", "random_error", "random_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False, local_options={"retries": 2})
    assert ret.input_data["provenance"]["retries"] == 2
    assert ret.error.error_type == "random_error"

    failure_engine.iter_modes = ["random_error", "input_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False, local_options={"retries": 4})
    assert ret.input_data["provenance"]["retries"] == 1
    assert ret.error.error_type == "input_error"


def test_random_failure_with_success(failure_engine):

    failure_engine.iter_modes = ["random_error", "pass"]
    failure_engine.ncalls = 0
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False, local_options={"retries": 1})

    assert ret.success, ret.error.error_message
    assert ret.provenance.retries == 1
    assert ret.extras["ncalls"] == 2
