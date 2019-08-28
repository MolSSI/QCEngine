import tempfile
import os

from qcengine import testing
from qcengine import util
from qcengine import cli
from qcengine import get_molecule

from qcelemental.models import ResultInput, OptimizationInput


def test_info():
    """Test for qcengine info"""
    for arg in cli.info_choices:
        args = ["qcengine", "info", arg]
        assert util.execute(args)[0] is True

    args = ["qcengine", "info"]
    assert util.execute(args)[0] is True


@testing.using_psi4
def test_run_psi4():
    """Tests qcengine run with psi4 and JSON input"""
    inp = ResultInput(molecule=get_molecule("hydrogen"),
                      driver="energy",
                      model={"method": "hf", "basis": "6-31G"})

    args = ["qcengine", "run", "psi4", inp.json()]
    assert util.execute(args)[0] is True

    args = ["qcengine", "run", "psi4", "input.json"]
    assert util.execute(args, infiles={"input.json": inp.json()})[0] is True


@testing.using_geometric
@testing.using_psi4
def test_run_procedure():
    """Tests qcengine run-procedure with geometric, psi4, and JSON input"""
    inp = {"schema_name": "qcschema_optimization_input",
           "schema_version": 1,
           "keywords": {
               "coordsys": "tric",
               "maxiter": 100,
               "program": "psi4"
           },
           "input_specification": {
               "schema_name": "qcschema_input",
               "schema_version": 1,
               "driver": "gradient",
               "model": {"method": "HF", "basis": "sto-3g"},
               "keywords": {},
           },
           "initial_molecule": get_molecule("hydrogen")}
    inp = OptimizationInput(**inp)

    args = ["qcengine", "run-procedure", "geometric", inp.json()]
    assert util.execute(args)[0] is True

    args = ["qcengine", "run-procedure", "geometric", "input.json"]
    assert util.execute(args, infiles={"input.json": inp.json()})[0] is True

