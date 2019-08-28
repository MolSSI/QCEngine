import json
import os
import shutil

from qcengine import testing
from qcengine import util
from qcengine import cli
from qcengine import get_molecule

from qcelemental.models import ResultInput, OptimizationInput

cov = shutil.which('coverage')
cov_cmds = []
if cov:
    cov_cmds = [cov, "run", "--append", "--source=" + os.path.dirname(os.path.abspath(__file__))]


def test_info():
    """Test for qcengine info"""
    outputs = []
    for arg in cli.info_choices:
        args = ["qcengine", "info", arg]
        res = util.execute(args)
        assert res[0] is True
        if arg not in {"all", "config"}:  # output of config changes call-to-call depending e.g. on mem available
            outputs.append(res[1]["stdout"])
        util.execute(cov_cmds + args)

    args = ["qcengine", "info"]
    res = util.execute(args)
    assert res[0] is True

    for output in outputs:
        assert output in res[1]["stdout"]
    util.execute(cov_cmds + args)


@testing.using_psi4
def test_run_psi4():
    """Tests qcengine run with psi4 and JSON input"""
    def check_result(result):
        assert result[0] is True
        output = json.loads(result[1]["stdout"])
        assert output["provenance"]["creator"].lower() == "psi4"
        assert output["success"] is True

    inp = ResultInput(molecule=get_molecule("hydrogen"),
                      driver="energy",
                      model={"method": "hf", "basis": "6-31G"})

    args = ["qcengine", "run", "psi4", inp.json()]
    check_result(util.execute(args))
    util.execute(cov_cmds + args)

    args = ["qcengine", "run", "psi4", "input.json"]
    check_result(util.execute(args, infiles={"input.json": inp.json()}))
    util.execute(cov_cmds + args)


@testing.using_geometric
@testing.using_psi4
def test_run_procedure():
    """Tests qcengine run-procedure with geometric, psi4, and JSON input"""
    def check_result(result):
        assert result[0] is True
        output = json.loads(result[1]["stdout"])
        print(output)
        assert output["provenance"]["creator"].lower() == "geometric"
        assert output["success"] is True

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
    check_result(util.execute(args))
    util.execute(cov_cmds + args)

    args = ["qcengine", "run-procedure", "geometric", "input.json"]
    check_result(util.execute(args, infiles={"input.json": inp.json()}))
    util.execute(cov_cmds + args)

