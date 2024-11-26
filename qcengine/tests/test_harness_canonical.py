"""
Tests the DQM compute dispatch module
"""
import copy
import pprint

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


@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_compute_gradient(program, model, keywords, schema_versions, request):
    models, retver, _ = schema_versions

    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    molecule = _get_molecule(program, models.Molecule)

    if from_v2(request.node.name):
        inp = models.AtomicInput(
            molecule=molecule,
            extras={"mytag": "something"},
            specification=models.AtomicSpecification(
                driver="gradient", model=model, extras={"myspectag": "somethingelse"}, keywords=keywords
            ),
        )
    else:
        inp = models.AtomicInput(
            molecule=molecule,
            driver="gradient",
            model=model,
            extras={"mytag": "something", "myspectag": "somethingelse"},
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
            assert "mytag" in ret.input_data.extras, ret.input_data.extras
            assert "mytag" not in ret.extras, "input extras wrongly present in result"
            # assert "myspectag" in ret.input_data.specification.extras, ret.input_data.specification.extras
            # assert "myspectag" not in ret.extras, "input spec extras wrongly present in result"
        else:
            assert "mytag" in ret.extras, ret.extras
            assert "myspectag" in ret.extras, ret.extras


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
            # TODO if "v2" in request.node.name or "to_v1" in request.node.name:
            # TODO     assert ret["input_data"]["specification"]["model"]["method"] == model["method"], "input not copied over"
            # TODO else:
            # TODO     assert ret["input_data"]["model"]["method"] == model["method"], "input not copied over"
        else:
            assert ret.success is False, "wrongly successful"
            assert isinstance(ret, (qcel.models.v1.FailedOperation, qcel.models.v2.FailedOperation)), "wrong class"
            assert ret.error.error_type == "input_error", f"wrong type: {ret.error.error_type=} != 'input_error'"
            # TODO if "v2" in request.node.name:
            # TODO     ret_method = ret.input_data.specification.model.method
            # TODO elif "to_v1" in request.node.name:
            # TODO     ret_method = ret.input_data["specification"]["model"]["method"]
            # TODO else:
            # TODO     ret_method = ret.input_data["model"]["method"]
            # TODO assert ret_method == model["method"], "input not copied over"


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
