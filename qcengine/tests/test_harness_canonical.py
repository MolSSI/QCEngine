"""
Tests the DQM compute dispatch module
"""
import copy
import sys
import warnings

import msgpack
import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.tests.test_model_results import center_data

import qcengine as qcng
from qcengine.testing import checkver_and_convert, from_v2, has_program, schema_versions, schema_versions2, using

qcsk_bs = {"name": "custom_basis", "center_data": center_data, "atom_map": ["bs_sto3g_h", "bs_sto3g_h"]}

_canonical_methods = [
    ("dftd3", {"method": "b3lyp-d3"}, {}),
    ("qcore", {"method": "pbe", "basis": "6-31G"}, {}),
    ("molpro", {"method": "hf", "basis": "6-31G"}, {}),
    ("mopac", {"method": "PM6"}, {}),
    ("mp2d", {"method": "MP2-DMP2"}, {}),
    ("nwchem", {"method": "hf", "basis": "6-31G"}, {}),
    ("openmm", {"method": "openff-1.0.0", "basis": "smirnoff"}, {}),
    ("psi4", {"method": "hf", "basis": "6-31G"}, {}),
    ("qchem", {"method": "hf", "basis": "6-31G"}, {}),
    ("rdkit", {"method": "UFF"}, {}),
    ("terachem_pbs", {"method": "b3lyp", "basis": "6-31G"}, {}),
    ("torchani", {"method": "ANI1x"}, {}),
    ("turbomole", {"method": "pbe", "basis": "6-31G"}, {}),
    ("xtb", {"method": "GFN2-xTB"}, {}),
    ("adcc", {"method": "adc2", "basis": "6-31G"}, {"n_triplets": 3}),
    ("gcp", {"method": "hf3c"}, {}),
    ("mrchem", {"method": "blyp"}, {"world_prec": 1.0e-3}),
    ("cfour", {"method": "hf", "basis": "6-31G"}, {}),
    ("gamess", {"method": "hf", "basis": "n31"}, {"basis__NGAUSS": 6}),
    ("mctc-gcp", {"method": "dft/sv"}, {}),
    ("mace", {"method": "small"}, {}),
    ("aimnet2", {"method": "b973c"}, {}),
    ("s-dftd3", {"method": "b3lyp-d3"}, {}),
    ("dftd4", {"method": "b3lyp-d4"}, {}),
    # add as programs available
    # ("terachem", {"method": "bad"}),
]

_canonical_methods_qcsk_basis = [
    ("adcc", {"method": "adc2", "basis": qcsk_bs}, {"n_triplets": 3}),
    ("cfour", {"method": "hf", "basis": qcsk_bs}, {}),
    ("gamess", {"method": "hf", "basis": qcsk_bs}, {}),
    ("molpro", {"method": "hf", "basis": qcsk_bs}, {}),
    ("nwchem", {"method": "hf", "basis": qcsk_bs}, {}),
    ("openmm", {"method": "openff-1.0.0", "basis": qcsk_bs}, {}),
    pytest.param("psi4", {"method": "hf", "basis": qcsk_bs}, {}, marks=using("psi4_mp2qcsk")),
    ("qchem", {"method": "hf", "basis": qcsk_bs}, {}),
    ("qcore", {"method": "pbe", "basis": qcsk_bs}, {}),
    ("turbomole", {"method": "pbe", "basis": qcsk_bs}, {}),
]


def _get_molecule(program, molcls):
    if program in ["openmm", "terachem_pbs"]:
        dmol = qcng.get_molecule("water", return_dict=True)
    else:
        dmol = qcng.get_molecule("hydrogen", return_dict=True)
    return molcls(**dmol)


@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_compute_energy(program, model, keywords, schema_versions, request):
    models, retver, _ = schema_versions

    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    molecule = _get_molecule(program, models.Molecule)

    if from_v2(request.node.name):
        inp = models.AtomicInput(
            molecule=molecule, specification=models.AtomicSpecification(driver="energy", model=model, keywords=keywords)
        )
    else:
        inp = models.AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)

    inp = checkver_and_convert(inp, request.node.name, "pre")
    ret = qcng.compute(inp, program, raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success is True
    assert isinstance(ret.return_result, float)
    assert ret.return_result == ret.properties.return_energy

    # official leave this as dict(), not model_dump(), to ensure remains operational
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dret = ret.dict()
    assert dret["return_result"] == dret["properties"]["return_energy"]


@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_compute_gradient(program, model, keywords, schema_versions, request):
    models, retver, _ = schema_versions

    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    molecule = _get_molecule(program, models.Molecule)

    if from_v2(request.node.name):
        inp = models.AtomicInput(
            molecule=molecule,
            specification=models.AtomicSpecification(
                driver="gradient",
                model=model,
                extras={"mytag": "something"},
                keywords=keywords,
            ),
        )
    else:
        inp = models.AtomicInput(
            molecule=molecule,
            driver="gradient",
            model=model,
            extras={"mytag": "something"},
            keywords=keywords,
        )
    if program in ["adcc"]:
        inp = checkver_and_convert(inp, request.node.name, "pre")
        with pytest.raises(qcng.exceptions.InputError) as e:
            qcng.compute(inp, program, raise_error=True, return_version=retver)

        assert "gradient not implemented" in str(e.value)

    else:
        inp = checkver_and_convert(inp, request.node.name, "pre")
        ret = qcng.compute(inp, program, raise_error=True, return_version=retver)
        ret = checkver_and_convert(ret, request.node.name, "post")

        assert ret.success is True
        assert isinstance(ret.return_result, np.ndarray)
        assert len(ret.return_result.shape) == 2
        assert ret.return_result.shape[1] == 3
        if "v2" in request.node.name:
            assert "mytag" in ret.input_data.specification.extras, ret.input_data.specification.extras
            assert "mytag" not in ret.extras, "input extras wrongly present in result"
        else:
            assert "mytag" in ret.extras, ret.extras


@pytest.mark.parametrize("program, model, keywords", _canonical_methods_qcsk_basis)
def test_compute_energy_qcsk_basis(program, model, keywords, schema_versions, request):
    models, retver, _ = schema_versions

    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    molecule = _get_molecule(program, models.Molecule)
    if from_v2(request.node.name):
        inp = models.AtomicInput(
            molecule=molecule, specification=models.AtomicSpecification(driver="energy", model=model, keywords=keywords)
        )
    else:
        inp = models.AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        inp = checkver_and_convert(inp, request.node.name, "pre")
        res = qcng.compute(inp, program, raise_error=True, return_version=retver)
        checkver_and_convert(res, request.node.name, "post")

    assert "QCSchema BasisSet for model.basis not implemented" in str(e.value)


@pytest.mark.parametrize(
    "program, model",
    [
        ("cfour", {"method": "bad"}),
        ("dftd3", {"method": "bad"}),
        ("dftd3", {"method": "b3lyp-d3", "driver": "hessian"}),
        ("qcore", {"method": "bad"}),
        ("gamess", {"method": "bad"}),
        ("mopac", {"method": "bad"}),
        ("mp2d", {"method": "bad"}),
        ("nwchem", {"method": "bad"}),
        ("openmm", {"method": "bad"}),
        ("psi4", {"method": "bad"}),
        ("qchem", {"method": "bad"}),
        ("rdkit", {"method": "bad"}),
        ("terachem_pbs", {"method": "bad"}),
        ("torchani", {"method": "bad"}),
        ("turbomole", {"method": "bad"}),
        ("adcc", {"method": "bad"}),
        ("gcp", {"method": "bad"}),
        ("mrchem", {"method": "bad"}),
        ("mctc-gcp", {"method": "bad"}),
        ("mace", {"method": "bad"})
        # add as programs available
        # ("molpro", {"method": "bad"}),
        # ("terachem", {"method": "bad"}),
        # ("xtb", {"method": "bad"}),
    ],
)
@pytest.mark.parametrize("raiserr", [False, True])
@pytest.mark.parametrize("retdict", [False, True])
@pytest.mark.parametrize("inpdict", [False, True])
def test_compute_bad_models(program, model, schema_versions, request, raiserr, retdict, inpdict):
    models, retver, _ = schema_versions

    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    amodel = copy.deepcopy(model)
    adriver = amodel.pop("driver", "energy")
    if from_v2(request.node.name):
        inp = {
            "molecule": models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
            "specification": {"model": amodel, "driver": adriver},
        }
    else:
        inp = {
            "molecule": models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
            "model": amodel,
            "driver": adriver,
        }
    if not inpdict:
        inp = models.AtomicInput(**inp)
    inp = checkver_and_convert(inp, request.node.name, "pre")

    if raiserr:
        with pytest.raises(qcng.exceptions.InputError) as exc:
            qcng.compute(inp, program, raise_error=raiserr, return_dict=retdict, return_version=retver)
    else:
        ret = qcng.compute(inp, program, raise_error=raiserr, return_dict=retdict, return_version=retver)
        ret = checkver_and_convert(
            ret, request.node.name, "post", vercheck=False, cast_dict_as="FailedOperation"
        )  # TODO release vercheck=F?
        if retdict:
            assert ret["success"] is False, "wrongly successful"
            assert (
                ret["error"]["error_type"] == "input_error"
            ), f"wrong type: {ret['error']['error_type']=} != 'input_error'"
            if "v2" in request.node.name or "to_v1" in request.node.name:
                assert ret["input_data"]["specification"]["model"]["method"] == model["method"], "input not copied over"
            else:
                assert ret["input_data"]["model"]["method"] == model["method"], "input not copied over"
        else:
            assert ret.success is False, "wrongly successful"
            assert isinstance(
                ret, (qcel.models.v1.FailedOperation, qcel.models.v2.FailedOperation, qcel.models._v1v2.FailedOperation)
            ), "wrong class"
            assert ret.error.error_type == "input_error", f"wrong type: {ret.error.error_type=} != 'input_error'"
            if "v2" in request.node.name:
                ret_method = ret.input_data.specification.model.method
            elif "to_v1" in request.node.name:
                ret_method = ret.input_data["specification"]["model"]["method"]
            else:
                ret_method = ret.input_data["model"]["method"]
            assert ret_method == model["method"], "input not copied over"


@using("nwchem")
@pytest.mark.parametrize("envshim", [True, False])
@pytest.mark.parametrize("model", ["Atomic"])
@pytest.mark.parametrize(
    "goodcalc,input_data,return_version,return_dict,xptd",
    [
        # fmt: off
    #           good? in ret  dict?    py313                      py314                      py314+V1V2_SHIM
    pytest.param(True, 1, -1, False,  ("v1.Result",               "v2.FailedOp pyd",         "v2.FailedOp pyd")),          # xptd0
    pytest.param(True, 1,  2, False,  ("v2.Result",               "v2.Result",               "v2.Result")),                # xptd1
    pytest.param(True, 1, -1, True,   ("dict of v1.Result",       "dict of v2.FailedOp pyd", "dict of v1.Result")),        # xptd2
    pytest.param(True, 1,  2, True,   ("dict of v2.Result",       "dict of v2.Result",       "dict of v2.Result")),        # xptd3

    pytest.param(True, 2, -1, False,  ("v2.Result",               "v2.Result",               "v2.Result")),                # xptd4
    pytest.param(True, 2,  1, False , ("v1.Result",               "v2.FailedOp pyd",         "v2.FailedOp pyd")),          # xptd5
    pytest.param(True, 2, -1, True,   ("dict of v2.Result",       "dict of v2.Result",       "dict of v2.Result")),        # xptd6
    pytest.param(True, 2,  1, True,   ("dict of v1.Result",       "dict of v2.FailedOp pyd", "dict of v1.Result")),        # xptd7

    pytest.param(False, 1, -1, False, ("v1.FailedOp mtd",         "v2.FailedOp pyd",         "v2.FailedOp pyd")),          # xptd8
    pytest.param(False, 1,  2, False, ("v2.FailedOp mtd",         "v2.FailedOp mtd",         "v2.FailedOp mtd")),          # xptd9
    pytest.param(False, 1, -1, True,  ("dict of v1.FailedOp mtd", "dict of v2.FailedOp pyd", "dict of v1.FailedOp mtd")),  # xptd10
    pytest.param(False, 1,  2, True,  ("dict of v2.FailedOp mtd", "dict of v2.FailedOp mtd", "dict of v2.FailedOp mtd")),  # xptd11

    pytest.param(False, 2, -1, False, ("v2.FailedOp mtd",         "v2.FailedOp mtd",         "v2.FailedOp mtd")),          # xptd12
    pytest.param(False, 2,  1, False, ("v1.FailedOp mtd",         "v2.FailedOp pyd",         "v2.FailedOp pyd")),          # xptd13
    pytest.param(False, 2, -1, True,  ("dict of v2.FailedOp mtd", "dict of v2.FailedOp mtd", "dict of v2.FailedOp mtd")),  # xptd14
    pytest.param(False, 2,  1, True,  ("dict of v1.FailedOp mtd", "dict of v2.FailedOp pyd", "dict of v1.FailedOp mtd")),
        # xptd15
        # fmt: on
    ],
)
def test_compute_output_table(goodcalc, input_data, return_version, return_dict, xptd, model, envshim, monkeypatch):
    h2 = qcng.get_molecule("hydrogen", return_dict=True)
    method = "hf" if goodcalc else "standard_model_of_particle_physics"

    if model == "Atomic":
        if input_data == 1:
            inp = {"molecule": h2, "driver": "energy", "model": {"method": method, "basis": "sto-3g"}}
        elif input_data == 2:
            inp = {
                "molecule": h2,
                "specification": {"driver": "energy", "model": {"method": method, "basis": "sto-3g"}},
            }

    if sys.version_info >= (3, 14):
        if envshim:
            monkeypatch.setenv("QCNG_USE_V1V2_SHIM", "1")
            result_this_job = xptd[2]
        else:
            result_this_job = xptd[1]
    else:
        result_this_job = xptd[0]

    ret = qcng.compute(inp, "nwchem", raise_error=False, return_dict=return_dict, return_version=return_version)

    if model == "Atomic":
        # pyd_phrase covers _MSG314 from qcng/compute.py and _MSG2 from qcel/models/v1/__init__.py
        # Append ". You can (a)" to catch the former. use ". Use qcel" to catch the latter
        pyd_phrase = "Reason: pydantic.v1 is unavailable on Python 3.14+"
        mtd_phrase = "Method not recognized"

        if result_this_job == "v2.Result":
            assert isinstance(ret, qcel.models.v2.AtomicResult)

        elif result_this_job == "dict of v2.Result":
            assert isinstance(ret, dict)
            assert ret["schema_name"] == "qcschema_atomic_result"
            assert ret["schema_version"] == 2

        elif result_this_job == "v1.Result":
            assert isinstance(ret, qcel.models.v1.AtomicResult)

        elif result_this_job == "dict of v1.Result":
            assert isinstance(ret, dict)
            assert ret["schema_name"] == "qcschema_output"
            assert ret["schema_version"] == 1

        elif result_this_job.startswith("v2.FailedOp"):
            assert isinstance(ret, qcel.models.v2.FailedOperation)
            errmsg = {"v2.FailedOp mtd": mtd_phrase, "v2.FailedOp pyd": pyd_phrase}[result_this_job]
            assert errmsg in ret.error.error_message

        elif result_this_job.startswith("dict of v2.FailedOp"):
            assert isinstance(ret, dict)
            assert ret["schema_name"] == "qcschema_failed_operation"
            assert ret["schema_version"] == 2
            errmsg = {"dict of v2.FailedOp mtd": mtd_phrase, "dict of v2.FailedOp pyd": pyd_phrase}[result_this_job]
            assert errmsg in ret["error"]["error_message"]

        elif result_this_job.startswith("v1.FailedOp"):
            assert isinstance(ret, qcel.models.v1.FailedOperation)
            errmsg = {"v1.FailedOp mtd": mtd_phrase}[result_this_job]
            assert errmsg in ret.error.error_message

        elif result_this_job.startswith("dict of v1.FailedOp"):
            assert isinstance(ret, dict)
            assert "schema_name" not in ret
            assert "schema_version" not in ret
            errmsg = {"dict of v1.FailedOp mtd": mtd_phrase}[result_this_job]
            assert errmsg in ret["error"]["error_message"]

        else:
            assert 0


@pytest.mark.parametrize(
    "program, model",
    [
        ("psi4", {"method": "hf", "driver": "eighth"}),
    ],
)
@pytest.mark.parametrize("raiserr", [False, True])
@pytest.mark.parametrize("retdict", [False, True])
def test_compute_badder_models(program, model, schema_versions2, request, raiserr, retdict):
    models, retver, _ = schema_versions2

    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    amodel = copy.deepcopy(model)
    adriver = amodel.pop("driver", "energy")
    if from_v2(request.node.name):
        inp = {
            "molecule": models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
            "specification": {"model": amodel, "driver": adriver},
        }
    else:
        inp = {
            "molecule": models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
            "model": amodel,
            "driver": adriver,
        }
    # inp = checkver_and_convert(inp, request.node.name, "pre")

    if raiserr:
        with pytest.raises(qcng.exceptions.InputError) as exc:
            qcng.compute(inp, program, raise_error=raiserr, return_dict=retdict, return_version=retver)
    else:
        ret = qcng.compute(inp, program, raise_error=raiserr, return_dict=retdict, return_version=retver)
        if retdict:
            assert ret["success"] is False, "wrongly successful"
            assert (
                ret["error"]["error_type"] == "input_error"
            ), f"wrong type: {ret['error']['error_type']=} != 'input_error'"
            if from_v2(request.node.name):
                assert ret["input_data"]["specification"]["driver"] == "eighth", "input not copied over"
            else:
                assert ret["input_data"]["driver"] == "eighth", "input not copied over"
        else:
            assert ret.success is False, "wrongly successful"
            assert isinstance(ret, (qcel.models.v1.FailedOperation, qcel.models.v2.FailedOperation)), "wrong class"
            assert ret.error.error_type == "input_error", f"wrong type: {ret.error.error_type=} != 'input_error'"
            # note that input_data *always* a dict in this test (even for v2)
            #   since the error is that the AtomicInput model can't be constructed
            if from_v2(request.node.name):
                assert ret.input_data["specification"]["driver"] == "eighth", "input not copied over"
            else:
                assert ret.input_data["driver"] == "eighth", "input not copied over"


def test_psi4_restarts(monkeypatch, schema_versions, request):
    """
    Make sure that a random error is raised which can be restarted if psi4 fails with no error message
    """
    models, retver, _ = schema_versions

    if not has_program("psi4"):
        pytest.skip("Program psi4 not found.")

    # create the psi4 task
    if from_v2(request.node.name):
        inp = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
            specification={"driver": "energy", "model": {"method": "hf", "basis": "6-31G"}},
        )
    else:
        inp = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
            driver="energy",
            model={"method": "hf", "basis": "6-31G"},
        )

    def mock_execute(*args, **kwargs):
        """
        Mock the output of a failed psi4 task with missing error message.
        """

        mock_output = {"sucess": False, "outfiles": {"data.msgpack": msgpack.dumps({"missing": "data"})}}
        return True, mock_output

    monkeypatch.setattr("qcengine.programs.psi4.execute", mock_execute)

    inp = checkver_and_convert(inp, request.node.name, "pre")
    with pytest.raises(qcng.exceptions.RandomError):
        _ = qcng.compute(
            input_data=inp, program="psi4", raise_error=True, task_config={"retries": 0}, return_version=retver
        )
