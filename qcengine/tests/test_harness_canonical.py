"""
Tests the DQM compute dispatch module
"""
import copy

import msgpack
import numpy as np
import pytest
from qcelemental.tests.test_model_results import center_data

import qcengine as qcng
from qcengine.testing import checkver_and_convert, has_program, schema_versions, using

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
    models, _ = schema_versions

    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    molecule = _get_molecule(program, models.Molecule)

    inp = models.AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)

    inp = checkver_and_convert(inp, request.node.name, "pre")
    ret = qcng.compute(inp, program, raise_error=True)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success is True
    assert isinstance(ret.return_result, float)


@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_compute_gradient(program, model, keywords, schema_versions, request):
    models, _ = schema_versions

    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    molecule = _get_molecule(program, models.Molecule)

    inp = models.AtomicInput(
        molecule=molecule, driver="gradient", model=model, extras={"mytag": "something"}, keywords=keywords
    )
    if program in ["adcc"]:
        inp = checkver_and_convert(inp, request.node.name, "pre")
        with pytest.raises(qcng.exceptions.InputError) as e:
            qcng.compute(inp, program, raise_error=True)

        assert "gradient not implemented" in str(e.value)

    else:
        inp = checkver_and_convert(inp, request.node.name, "pre")
        ret = qcng.compute(inp, program, raise_error=True)
        ret = checkver_and_convert(ret, request.node.name, "post")

        assert ret.success is True
        assert isinstance(ret.return_result, np.ndarray)
        assert len(ret.return_result.shape) == 2
        assert ret.return_result.shape[1] == 3
        assert "mytag" in ret.extras, ret.extras


@pytest.mark.parametrize("program, model, keywords", _canonical_methods_qcsk_basis)
def test_compute_energy_qcsk_basis(program, model, keywords, schema_versions, request):
    models, _ = schema_versions

    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    molecule = _get_molecule(program, models.Molecule)
    inp = models.AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)

    with pytest.raises(qcng.exceptions.InputError) as e:
        inp = checkver_and_convert(inp, request.node.name, "pre")
        res = qcng.compute(inp, program, raise_error=True)
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
def test_compute_bad_models(program, model, schema_versions, request):
    models, _ = schema_versions

    if not has_program(program):
        pytest.skip("Program '{}' not found.".format(program))

    amodel = copy.deepcopy(model)
    adriver = amodel.pop("driver", "energy")
    inp = models.AtomicInput(molecule=models.Molecule(**qcng.get_molecule("hydrogen",return_dict=True)), driver=adriver, model=amodel)

    inp = checkver_and_convert(inp, request.node.name, "pre")
    with pytest.raises(qcng.exceptions.InputError) as exc:
        ret = qcng.compute(inp, program, raise_error=True)


def test_psi4_restarts(monkeypatch, schema_versions, request):
    """
    Make sure that a random error is raised which can be restarted if psi4 fails with no error message
    """
    models, _ = schema_versions

    if not has_program("psi4"):
        pytest.skip("Program psi4 not found.")

    # create the psi4 task
    inp = models.AtomicInput(molecule=models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)), driver="energy", model={"method": "hf", "basis": "6-31G"})

    def mock_execute(*args, **kwargs):
        """
        Mock the output of a failed psi4 task with missing error message.
        """

        mock_output = {"sucess": False, "outfiles": {"data.msgpack": msgpack.dumps({"missing": "data"})}}
        return True, mock_output

    monkeypatch.setattr("qcengine.programs.psi4.execute", mock_execute)

    inp = checkver_and_convert(inp, request.node.name, "pre")
    with pytest.raises(qcng.exceptions.RandomError):
        _ = qcng.compute(input_data=inp, program="psi4", raise_error=True, task_config={"retries": 0})
