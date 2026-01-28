import json
import os
import subprocess
import sys
from typing import List

import pytest

from qcengine import cli, get_molecule, util
from qcengine.testing import checkver_and_convert, from_v2, schema_versions, using


def run_qcengine_cli(args: List[str], stdin: str = None) -> str:
    """
    Runs qcengine via its CLI

    This method was chosen over the more sophisticated util.execute in order to get pytest-cov to work.
    For the same reason, qcengine is invoked as "python -m qcengine.cli" rather than "qcengine".

    Parameters
    ----------
    args: List[str]
        List of CLI arguments.
    stdin: Optional[str]
        Standard input for the process.

    Returns
    -------
    str
        QCEngine CLI standard output.
    """
    if stdin is not None:
        stdin = stdin.encode("utf-8")

    return subprocess.check_output([sys.executable, "-m", "qcengine"] + args, input=stdin)


def test_no_args():
    """Test for qcengine with no arguments"""
    try:
        run_qcengine_cli([])
    except subprocess.CalledProcessError as e:
        assert e.returncode == 1


@pytest.mark.skipif(
    sys.version_info >= (3, 14), reason="too many programs throw errors on import with Py >=3.14. revisit soon!"
)
def test_info():
    """Test for qcengine info"""
    outputs = []
    for arg in cli.info_choices:
        output = run_qcengine_cli(["info", arg])
        if arg not in {"all", "config"}:  # output of config changes call-to-call depending e.g. on mem available
            outputs.append(output)

    default_output = run_qcengine_cli(["info"])
    for output in outputs:
        assert output in default_output


@using("psi4")
def test_run_psi4(tmp_path, schema_versions, request):
    """Tests qcengine run with psi4 and JSON input"""
    models, retver, _ = schema_versions

    def check_result(stdout):
        output = json.loads(stdout)
        output = checkver_and_convert(output, request.node.name, "post")
        assert output["provenance"]["creator"].lower() == "psi4"
        assert output["success"] is True

    if from_v2(request.node.name):
        inp = models.AtomicInput(
            molecule=models.Molecule(**get_molecule("hydrogen", return_dict=True)),
            specification=models.AtomicSpecification(
                driver="energy",
                model={"method": "hf", "basis": "6-31G"},
            ),
        )
    else:
        inp = models.AtomicInput(
            molecule=models.Molecule(**get_molecule("hydrogen", return_dict=True)),
            driver="energy",
            model={"method": "hf", "basis": "6-31G"},
        )

    inp = checkver_and_convert(inp, request.node.name, "pre")
    args = ["run", "psi4", f"--return-version={retver}", inp.model_dump_json()]
    check_result(run_qcengine_cli(args))

    args = ["run", "psi4", os.path.join(tmp_path, "input.json"), f"--return-version={retver}"]
    with util.disk_files({"input.json": inp.model_dump_json()}, {}, cwd=tmp_path):
        check_result(run_qcengine_cli(args))

    args = ["run", "psi4", f"--return-version={retver}", "-"]
    # model_dump_json() works on v1 or v2 (see above). below tests that json() still works on v1.
    if "as_v1" in request.node.name:
        check_result(run_qcengine_cli(args, stdin=inp.json()))
    else:
        check_result(run_qcengine_cli(args, stdin=inp.model_dump_json()))


@using("geometric")
@using("psi4")
def test_run_procedure(tmp_path, schema_versions, request):
    """Tests qcengine run-procedure with geometric, psi4, and JSON input"""
    models, retver, _ = schema_versions

    def check_result(stdout):
        output = json.loads(stdout)
        output = checkver_and_convert(output, request.node.name, "post", cast_dict_as="OptimizationResult")
        assert output["provenance"]["creator"].lower() == "geometric"
        assert output["success"] is True

    if from_v2(request.node.name):
        inp = {
            "specification": {
                "keywords": {"coordsys": "tric", "maxiter": 100},
                "specification": {
                    "driver": "gradient",
                    "model": {"method": "HF", "basis": "sto-3g"},
                    "keywords": {},
                    "program": "psi4",
                },
            },
            "initial_molecule": models.Molecule(**get_molecule("hydrogen", return_dict=True)),
        }
    else:
        inp = {
            "keywords": {"coordsys": "tric", "maxiter": 100, "program": "psi4"},
            "input_specification": {"driver": "gradient", "model": {"method": "HF", "basis": "sto-3g"}, "keywords": {}},
            "initial_molecule": models.Molecule(**get_molecule("hydrogen", return_dict=True)),
        }
    inp = models.OptimizationInput(**inp)

    inp = checkver_and_convert(inp, request.node.name, "pre")
    if "to_v" in request.node.name:
        args = ["run", "geometric", inp.model_dump_json(), f"--return-version={retver}"]
    else:
        args = ["run-procedure", "geometric", inp.model_dump_json()]
    check_result(run_qcengine_cli(args))

    if "to_v" in request.node.name:
        args = ["run", "geometric", f"--return-version={retver}", os.path.join(tmp_path, "input.json")]
    else:
        args = ["run-procedure", "geometric", os.path.join(tmp_path, "input.json")]
    with util.disk_files({"input.json": inp.model_dump_json()}, {}, cwd=tmp_path):
        check_result(run_qcengine_cli(args))

    # try unified route
    args = ["run", "geometric", inp.model_dump_json(), f"--return-version={retver}"]
    check_result(run_qcengine_cli(args, stdin=inp.model_dump_json()))
