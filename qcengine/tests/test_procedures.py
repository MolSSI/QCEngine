"""
Tests the DQM compute dispatch module
"""

import warnings

import numpy as np
import pytest
import qcelemental as qcel

import qcengine as qcng
from qcengine.testing import checkver_and_convert, failure_engine, from_v2, schema_versions, using


@pytest.fixture(scope="function")
def input_data(request):
    if from_v2(request.node.name):
        return {
            "specification": {
                "keywords": {"coordsys": "tric", "maxiter": 100},
                "specification": {"driver": "gradient", "model": None, "keywords": {}, "program": None},
            },
            "initial_molecule": None,
        }
    else:
        return {
            "keywords": {"coordsys": "tric", "maxiter": 100, "program": None},
            "input_specification": {"driver": "gradient", "model": None, "keywords": {}},
            "initial_molecule": None,
        }


@using("psi4")
@pytest.mark.parametrize("ncores", [1, 4])
@pytest.mark.parametrize(
    "optimizer",
    [
        pytest.param("geometric", marks=using("geometric")),
        pytest.param("optking", marks=using("optking")),
        pytest.param("berny", marks=using("berny")),
        pytest.param("nwchemdriver", marks=using("nwchem")),
    ],
)
def test_geometric_psi4(input_data, optimizer, ncores, schema_versions, request):
    models, retver, _ = schema_versions

    input_data["initial_molecule"] = models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True))
    if optimizer == "nwchemdriver":
        grad_program = "nwchem"
        grad_kw = {}
    else:
        grad_program = "psi4"
        grad_kw = {"scf_properties": ["wiberg_lowdin_indices"]}

    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
        input_data["specification"]["specification"]["keywords"] = grad_kw
        input_data["specification"]["specification"]["program"] = grad_program
        input_data["specification"]["specification"]["extras"] = {"myqctag": "hello psi4"}
        input_data["specification"]["extras"] = {"myopttag": "hello qcengine"}
        input_data["specification"]["protocols"] = {"trajectory_results": "all"}

    else:
        input_data["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
        input_data["input_specification"]["keywords"] = grad_kw
        input_data["keywords"]["program"] = grad_program
        input_data["input_specification"]["extras"] = {"myqctag": "hello psi4"}
        input_data["extras"] = {"myopttag": "hello qcengine"}

    input_data = models.OptimizationInput(**input_data)

    task_config = {
        "ncores": ncores,
    }

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        ret = qcng.compute_procedure(
            input_data, optimizer, raise_error=True, task_config=task_config, return_version=retver
        )
    ret = checkver_and_convert(ret, request.node.name, "post")
    trajs_tgt = ret.trajectory_results if "v2" in request.node.name else ret.trajectory

    assert 10 > len(trajs_tgt) > 1

    assert pytest.approx(ret.final_molecule.measure([0, 1]), 1.0e-4) == 1.3459150737
    assert ret.provenance.creator.lower() == optimizer
    assert trajs_tgt[0].provenance.creator.lower() == grad_program

    # Note: thread passing is semi-implemented in optking and subject to other env factors
    # if optimizer == "optking" and ncores != 1:
    #     with pytest.raises(AssertionError):
    #         # Optking not passing threads to psi4
    #         assert trajs_tgt[0].provenance.nthreads == ncores
    if optimizer in ["optking", "nwchemdriver"]:
        pass
        # threading not yet implemented in NWChemDriver
    else:
        assert trajs_tgt[0].provenance.nthreads == ncores

    # Check keywords passing
    if optimizer != "nwchemdriver":
        for single in trajs_tgt:
            if "v2" in request.node.name:
                assert "scf_properties" in single.input_data.specification.keywords
            else:
                assert "scf_properties" in single.keywords
            assert (
                "WIBERG_LOWDIN_INDICES" in single.extras["qcvars"] or "WIBERG LOWDIN INDICES" in single.extras["qcvars"]
            )
            # TODO: old WIBERG qcvar used underscore; new one uses space. covering bases here but remove someday

    # Check extras passing
    if "v2" in request.node.name:
        assert (
            "myqctag" in ret.input_data.specification.specification.extras
        ), ret.input_data.specification.specification.extras.keys()
        if optimizer != "geometric":
            # geometric: below can't happen until optimizers call gradients with v2
            assert "myqctag" not in ret.trajectory_results[0].extras, "input extras wrongly present in sub-result"
        assert "myopttag" in ret.input_data.specification.extras, ret.input_data.specification.extras.keys()
        assert "myopttag" not in ret.extras, "input extras wrongly present in top-result"
    else:
        if optimizer != "optking":
            assert "myqctag" in ret.trajectory[0].extras, ret.trajectory[0].extras.keys()
            assert "myopttag" in ret.extras, ret.extras.keys()


@using("psi4")
@pytest.mark.parametrize(
    "optimizer",
    [
        pytest.param("geometric", marks=using("geometric")),
        pytest.param("optking", marks=using("optking")),
        pytest.param("berny", marks=using("berny")),
    ],
)
def test_geometric_local_options(input_data, schema_versions, request, optimizer):
    models, retver, _ = schema_versions

    input_data["initial_molecule"] = models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True))

    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
        input_data["specification"]["specification"]["program"] = "psi4"
        input_data["specification"]["protocols"] = {"trajectory_results": "final"}
    else:
        input_data["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
        input_data["keywords"]["program"] = "psi4"

    input_data = models.OptimizationInput(**input_data)

    # Set some extremely large number to test
    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, optimizer, raise_error=True, task_config={"memory": "5000"}, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    trajs_tgt = ret.trajectory_results if "v2" in request.node.name else ret.trajectory
    atin_tgt = ret.input_data.specification.specification if "v2" in request.node.name else ret.input_specification

    if optimizer == "optking":
        with pytest.raises(AssertionError):
            # Optking not passing memory to psi4
            assert pytest.approx(trajs_tgt[0].provenance.memory, 1) == 4900
    else:
        assert pytest.approx(trajs_tgt[0].provenance.memory, 1) == 4900

    # Make sure we cleaned up
    assert "_qcengine_local_config" not in atin_tgt.extras
    assert "_qcengine_local_config" not in trajs_tgt[0].extras
    if "v2" in request.node.name:
        assert "_qcengine_local_config" not in trajs_tgt[0].input_data.specification.extras


@pytest.mark.parametrize(
    "optimizer,converged",
    [
        pytest.param("geometric", "Converged!", marks=using("geometric")),
        pytest.param("optking", "Convergence check returned True", marks=using("optking")),
        pytest.param("berny", "All criteria matched", marks=using("berny")),
    ],
)
@pytest.mark.parametrize(
    "gradprog,gradmodel",
    [
        pytest.param("rdkit", {"method": "UFF", "basis": ""}, marks=using("rdkit")),
        pytest.param("psi4", {"method": "HF", "basis": "sto-3g"}, marks=using("psi4")),
    ],
)
def test_optimizer_stdout(optimizer, gradprog, gradmodel, converged, input_data, schema_versions, request):
    models, retver, _ = schema_versions

    input_data["initial_molecule"] = models.Molecule(**qcng.get_molecule("water", return_dict=True))

    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = gradmodel
        input_data["specification"]["specification"]["program"] = gradprog
        input_data["specification"]["specification"]["protocols"] = {"stdout": False}
        input_data["specification"]["protocols"] = {"trajectory_results": "all"}
    else:
        input_data["input_specification"]["model"] = gradmodel
        input_data["keywords"]["program"] = gradprog
        # no way to turn off gradient stdout in v1

    input_data = models.OptimizationInput(**input_data)

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, optimizer, raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success is True
    assert converged in ret.stdout

    trajs_tgt = ret.trajectory_results if "v2" in request.node.name else ret.trajectory
    assert ret.provenance.creator.lower() == optimizer
    assert trajs_tgt[0].provenance.creator.lower() == gradprog
    # rdkit does not have stdout, while psi4 does and is suppressed by protocol in v2 but can't be thru v1
    if from_v2(request.node.name) or gradprog == "rdkit":
        assert trajs_tgt[0].stdout is None
    else:
        assert trajs_tgt[0].stdout is not None

    if "v2" in request.node.name:
        assert ret.properties.optimization_iterations in [3, 4, 5]
        assert len(ret.trajectory_properties) == ret.properties.optimization_iterations
        assert np.allclose(ret.properties.return_gradient, ret.trajectory_results[-1].return_result)
        assert ret.trajectory_results[-1].properties.return_energy == pytest.approx(
            ret.properties.return_energy, 1.0e-4
        )
        if gradprog == "psi4":
            assert ret.properties.nuclear_repulsion_energy == pytest.approx(8.906, 1.0e-3)
            assert ret.properties.return_energy == pytest.approx(-74.96599, 1.0e-4)
        elif gradprog == "rdkit":
            assert ret.properties.nuclear_repulsion_energy == pytest.approx(8.888, 1.0e-3)
            assert ret.properties.return_energy == pytest.approx(0, 1.0e-2)  # 0 seems wrong


@pytest.mark.parametrize(
    "optimizer",
    [
        pytest.param("geometric", marks=using("geometric")),
        pytest.param("optking", marks=using("optking")),
        pytest.param("berny", marks=using("berny")),
    ],
)
@pytest.mark.parametrize(
    "gradprog,gradmodel",
    [
        pytest.param("rdkit", {"method": "UFF", "basis": ""}, marks=using("rdkit")),
        # pytest.param("psi4", {"method": "HF", "basis": "sto-3g"}, marks=using("psi4")),
    ],
)
@pytest.mark.parametrize("traj_ptcl", ["none", "final", "all", "default"])
def test_optimizer_protocols(optimizer, gradprog, gradmodel, input_data, schema_versions, request, traj_ptcl):
    models, retver, _ = schema_versions

    input_data["initial_molecule"] = models.Molecule(**qcng.get_molecule("water", return_dict=True))

    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = gradmodel
        input_data["specification"]["specification"]["program"] = gradprog
        if traj_ptcl != "default":
            input_data["specification"]["protocols"] = {"trajectory_results": traj_ptcl}
    else:
        input_data["input_specification"]["model"] = gradmodel
        input_data["keywords"]["program"] = gradprog
        if traj_ptcl != "default":
            input_data["protocols"] = {"trajectory": traj_ptcl}

    input_data = models.OptimizationInput(**input_data)

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, optimizer, raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success is True

    trajs_tgt = ret.trajectory_results if "v2" in request.node.name else ret.trajectory
    assert ret.provenance.creator.lower() == optimizer

    if "v2" in request.node.name:
        assert len(ret.trajectory_properties) > 2

    if traj_ptcl == "none":
        assert len(trajs_tgt) == 0
    elif traj_ptcl == "final":
        assert len(trajs_tgt) == 1
    elif traj_ptcl == "all":
        assert len(trajs_tgt) > 2
    elif traj_ptcl == "default":
        if from_v2(request.node.name):
            assert len(trajs_tgt) == 0
        else:
            assert len(trajs_tgt) > 2


@using("psi4")
@using("berny")
def test_berny_failed_gradient_computation(input_data, schema_versions, request):
    models, retver, _ = schema_versions

    input_data["initial_molecule"] = models.Molecule(**qcng.get_molecule("water", return_dict=True))
    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
        input_data["specification"]["specification"]["keywords"] = {"badpsi4key": "badpsi4value"}
        input_data["specification"]["specification"]["program"] = "psi4"
    else:
        input_data["input_specification"]["model"] = {"method": "HF", "basis": "sto-3g"}
        input_data["input_specification"]["keywords"] = {"badpsi4key": "badpsi4value"}
        input_data["keywords"]["program"] = "psi4"

    input_data = models.OptimizationInput(**input_data)

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "berny", raise_error=False, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    assert isinstance(ret, (qcel.models.v1.FailedOperation, qcel.models.v2.FailedOperation))
    assert ret.success is False
    assert ret.error.error_type == qcng.exceptions.InputError.error_type


@using("geometric")
@using("rdkit")
def test_geometric_rdkit_error(input_data, schema_versions, request):
    models, retver, _ = schema_versions

    input_data["initial_molecule"] = models.Molecule(**qcng.get_molecule("water", return_dict=True)).copy(
        exclude={"connectivity_"}
    )
    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = {"method": "UFF", "basis": ""}
        input_data["specification"]["specification"]["program"] = "rdkit"
    else:
        input_data["input_specification"]["model"] = {"method": "UFF", "basis": ""}
        input_data["keywords"]["program"] = "rdkit"

    input_data = models.OptimizationInput(**input_data)

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "geometric", return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    assert ret.success is False
    assert isinstance(ret.error.error_message, str)


@using("rdkit")
@pytest.mark.parametrize(
    "optimizer",
    [
        pytest.param("geometric", marks=using("geometric")),
        pytest.param("optking", marks=using("optking")),
        pytest.param("berny", marks=using("berny")),
        pytest.param("nwchemdriver", marks=using("nwchem")),
    ],
)
def test_optimization_protocols(optimizer, input_data, schema_versions, request):
    models, retver, _ = schema_versions

    if optimizer == "nwchemdriver":
        grad_program = "nwchem"
        grad_model = {"method": "HF", "basis": "sto-3g"}
    else:
        grad_program = "rdkit"
        grad_model = {"method": "UFF"}

    input_data["initial_molecule"] = models.Molecule(**qcng.get_molecule("water", return_dict=True))
    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = grad_model
        input_data["specification"]["specification"]["program"] = grad_program
        input_data["specification"]["protocols"] = {"trajectory_results": "initial_and_final"}
    else:
        input_data["input_specification"]["model"] = grad_model
        input_data["keywords"]["program"] = grad_program
        input_data["protocols"] = {"trajectory": "initial_and_final"}

    input_data = models.OptimizationInput(**input_data)

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, optimizer, raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    trajs_tgt = ret.trajectory_results if "v2" in request.node.name else ret.trajectory
    initmol_tgt = ret.input_data.initial_molecule if "v2" in request.node.name else ret.initial_molecule

    assert ret.success, ret.error.error_message
    assert len(trajs_tgt) == 2
    assert (
        pytest.approx(initmol_tgt.nuclear_repulsion_energy(), 1.0e-4)
        == trajs_tgt[0].molecule.nuclear_repulsion_energy()
    )
    if optimizer != "nwchemdriver":
        # in internal-opt mode, fix_/orient is not called, so the hash will be different
        assert initmol_tgt.get_hash() == trajs_tgt[0].molecule.get_hash()
    assert (
        pytest.approx(ret.final_molecule.nuclear_repulsion_energy(), 1.0e-4)
        == trajs_tgt[1].molecule.nuclear_repulsion_energy()
    )
    if optimizer not in ["nwchemdriver", "optking"]:
        # optking would pass if: optking/optimize.py: #computer.update_geometry(o_molsys.geom)
        assert ret.final_molecule.get_hash() == trajs_tgt[1].molecule.get_hash()


@using("geometric")
def test_geometric_retries(failure_engine, input_data, schema_versions, request):
    import geometric

    tric_ver = geometric.__version__
    models, retver, _ = schema_versions

    failure_engine.iter_modes = ["random_error", "pass", "random_error", "random_error", "pass"]  # Iter 1  # Iter 2
    failure_engine.iter_modes.extend(["pass"] * 20)

    input_data["initial_molecule"] = {
        "symbols": ["He", "He"],
        "geometry": [0, 0, 0, 0, 0, failure_engine.start_distance],
    }

    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = {"method": "something"}
        input_data["specification"]["specification"]["program"] = failure_engine.name
        input_data["specification"]["keywords"][
            "coordsys"
        ] = "cart"  # needed by geometric v1.0 to play nicely with failure_engine
        input_data["specification"]["protocols"] = {"trajectory_results": "all"}
    else:
        input_data["input_specification"]["model"] = {"method": "something"}
        input_data["keywords"]["program"] = failure_engine.name
        input_data["keywords"]["coordsys"] = "cart"  # needed by geometric v1.0 to play nicely with failure_engine

    input_data = models.OptimizationInput(**input_data)

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "geometric", task_config={"ncores": 13}, raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    trajs_tgt = ret.trajectory_results if "v2" in request.node.name else ret.trajectory

    assert ret.success is True
    assert trajs_tgt[0].provenance.retries == 1
    assert trajs_tgt[0].provenance.ncores == 13
    assert trajs_tgt[1].provenance.retries == 2
    assert trajs_tgt[1].provenance.ncores == 13
    assert "retries" not in trajs_tgt[2].provenance.model_dump()

    # Ensure we still fail
    failure_engine.iter_modes = ["random_error", "pass", "random_error", "random_error", "pass"]  # Iter 1  # Iter 2

    ret = qcng.compute(input_data, "geometric", task_config={"ncores": 13, "retries": 1}, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    assert ret.success is False
    assert ret.input_data["trajectory"][0]["provenance"]["retries"] == 1
    if tric_ver == "1.1":
        # bad! temp until https://github.com/leeping/geomeTRIC/pull/222 available
        assert len(ret.input_data["trajectory"]) == 1
    else:
        assert len(ret.input_data["trajectory"]) == 2


@using("geometric")
@pytest.mark.parametrize(
    "program, model, bench",
    [
        pytest.param(
            "rdkit", {"method": "UFF"}, [1.87130923886072, 2.959448636243545, 104.5099642579023], marks=using("rdkit")
        ),
        pytest.param(
            "rdkit",
            {"method": "mmff94"},
            [1.8310842343589573, 2.884612338953529, 103.93822919865106],
            marks=using("rdkit"),
        ),
        pytest.param(
            "rdkit",
            {"method": "MMFF94s"},
            [1.8310842343589573, 2.884612338953529, 103.93822919865106],
            marks=using("rdkit"),
        ),
        pytest.param(
            "torchani",
            {"method": "ANI1x"},
            [1.82581873750194, 2.866376526793269, 103.4332610730292],
            marks=using("torchani"),
        ),
        pytest.param(
            "mopac",
            {"method": "PM6"},
            [1.793052302291527, 2.893333237502448, 107.57254391453196],
            marks=using("mopac"),
        ),
        pytest.param(
            "openmm",
            {"method": "openff-1.0.0", "basis": "smirnoff"},
            [1.8344994291796748, 3.010099477501204, 110.25177977849998],
            marks=using("openmm"),
        ),
        pytest.param(
            "openmm",
            {"method": "openff_unconstrained-1.0.0", "basis": "smirnoff"},
            [1.8344994291195869, 3.0100994772976124, 110.25259556886984],
            marks=using("openmm"),
        ),
        pytest.param(
            "openmm",
            {"method": "smirnoff99Frosst-1.1.0", "basis": "smirnoff"},
            [1.814137087600702, 3.025566213038376, 112.9999999990053],
            marks=using("openmm"),
        ),
        pytest.param(
            "qcore",
            {"method": "GFN1"},
            [1.8104763949897031, 2.9132449420655213, 107.13403040879244],
            marks=using("qcore"),
        ),
    ],
)
def test_geometric_generic(input_data, program, model, bench, schema_versions, request):
    models, retver, _ = schema_versions

    input_data["initial_molecule"] = models.Molecule(**qcng.get_molecule("water", return_dict=True))
    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = model
        input_data["specification"]["specification"]["program"] = program
        input_data["specification"]["specification"]["extras"] = {
            "_secret_tags": {"mysecret_tag": "data1"}  # pragma: allowlist secret
        }
        input_data["specification"]["protocols"] = {"trajectory_results": "final"}
    else:
        input_data["input_specification"]["model"] = model
        input_data["keywords"]["program"] = program
        input_data["input_specification"]["extras"] = {
            "_secret_tags": {"mysecret_tag": "data1"}  # pragma: allowlist secret
        }

    input_data = checkver_and_convert(input_data, request.node.name, "pre", cast_dict_as="OptimizationInput")
    ret = qcng.compute(input_data, "geometric", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    trajs_tgt = ret.trajectory_results if "v2" in request.node.name else ret.trajectory

    assert ret.success is True
    assert "Converged!" in ret.stdout

    r01, r02, r12, a102 = ret.final_molecule.measure([[0, 1], [0, 2], [1, 2], [1, 0, 2]])

    assert pytest.approx(r01, 1.0e-4) == bench[0]
    assert pytest.approx(r02, 1.0e-4) == bench[0]
    assert pytest.approx(r12, 5.0e-4) == bench[1]
    assert pytest.approx(a102, 5.0e-4) == bench[2]

    assert "_secret_tags" in trajs_tgt[0].extras
    assert "data1" == trajs_tgt[0].extras["_secret_tags"]["mysecret_tag"], trajs_tgt[0].extras["_secret_tags"]
    # the _secret_tags shows up in Res.extras, not in Res.input_data.extras b/c geometric running v1 internally


@using("nwchem")
@pytest.mark.parametrize("linopt", [0, 1])
def test_nwchem_relax(linopt, schema_versions, request):
    models, retver, _ = schema_versions

    # Make the input file
    if from_v2(request.node.name):
        input_data = {
            "specification": {
                "specification": {
                    "model": {"method": "HF", "basis": "sto-3g"},
                    "keywords": {"set__driver:linopt": linopt},
                    "driver": "gradient",
                },
                "protocols": {"trajectory_results": "all"},
            },
            "initial_molecule": models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
        }
    else:
        input_data = {
            "input_specification": {
                "model": {"method": "HF", "basis": "sto-3g"},
                "keywords": {"set__driver:linopt": linopt},
            },
            "initial_molecule": models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
        }
    input_data = models.OptimizationInput(**input_data)

    # Run the relaxation
    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "nwchemdriver", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    trajs_tgt = ret.trajectory_results if "v2" in request.node.name else ret.trajectory

    assert 10 > len(trajs_tgt) > 1

    assert pytest.approx(ret.final_molecule.measure([0, 1]), 1.0e-4) == 1.3459150737


@using("nwchem")
def test_nwchem_restart(tmpdir, schema_versions, request):
    models, retver, _ = schema_versions

    # Make the input file
    if from_v2(request.node.name):
        input_data = {
            "specification": {
                "specification": {
                    "model": {"method": "HF", "basis": "sto-3g"},
                    "keywords": {"driver__maxiter": 2, "set__driver:linopt": 0},
                    "driver": "gradient",  # newly req'd by model in v2
                    "extras": {"allow_restarts": True},
                },
            },
            "initial_molecule": models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
        }
    else:
        input_data = {
            "input_specification": {
                "model": {"method": "HF", "basis": "sto-3g"},
                "keywords": {"driver__maxiter": 2, "set__driver:linopt": 0},
                "extras": {"allow_restarts": True},
            },
            "initial_molecule": models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True)),
        }
    input_data = models.OptimizationInput(**input_data)

    # Run an initial step, which should not converge
    local_opts = {"scratch_messy": True, "scratch_directory": str(tmpdir)}

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "nwchemdriver", task_config=local_opts, raise_error=False, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)
    assert not ret.success

    # Run it again, which should converge
    new_ret = qcng.compute(input_data, "nwchemdriver", task_config=local_opts, raise_error=True, return_version=retver)
    new_ret = checkver_and_convert(new_ret, request.node.name, "post")
    assert new_ret.success


@using("rdkit")
@using("torsiondrive")
@pytest.mark.parametrize("scan_ptcl", ["none", "all", "lowest"])
def test_torsiondrive_generic(schema_versions, request, scan_ptcl):
    models, retver, _ = schema_versions

    if from_v2(request.node.name):
        input_data = models.TorsionDriveInput(
            initial_molecule=[models.Molecule(**qcng.get_molecule("ethane", return_dict=True))] * 2,
            specification=models.TorsionDriveSpecification(
                keywords=models.TorsionDriveKeywords(dihedrals=[(2, 0, 1, 5)], grid_spacing=[180]),
                protocols=models.TorsionDriveProtocols(scan_results=scan_ptcl),
                specification=models.OptimizationSpecification(
                    protocols=models.OptimizationProtocols(trajectory_results="all"),
                    program="geomeTRIC",
                    keywords={
                        "coordsys": "hdlc",
                        "maxiter": 500,
                    },
                    specification=models.AtomicSpecification(
                        program="rdkit", driver=models.DriverEnum.gradient, model=models.Model(method="UFF", basis=None)
                    ),
                ),
            ),
        )
    else:
        input_data = models.TorsionDriveInput(
            keywords=models.TDKeywords(dihedrals=[(2, 0, 1, 5)], grid_spacing=[180]),
            input_specification=models.QCInputSpecification(
                driver=models.DriverEnum.gradient, model=models.Model(method="UFF", basis=None)
            ),
            initial_molecule=[models.Molecule(**qcng.get_molecule("ethane", return_dict=True))] * 2,
            optimization_spec=models.OptimizationSpecification(
                procedure="geomeTRIC",
                keywords={
                    "coordsys": "hdlc",
                    "maxiter": 500,
                    "program": "rdkit",
                },
            ),
        )

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "torsiondrive", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    if "_v2" not in request.node.name:
        assert ret.error is None
    assert ret.success

    expected_grid_ids = {"180", "0"}
    opthist_tgt = ret.scan_results if "v2" in request.node.name else ret.optimization_history

    if not (from_v2(request.node.name) and scan_ptcl == "none"):
        assert {*opthist_tgt} == expected_grid_ids

    if "v2" in request.node.name:
        assert {*ret.scan_properties} == expected_grid_ids
    else:
        assert {*ret.final_energies} == expected_grid_ids
    assert {*ret.final_molecules} == expected_grid_ids

    assert (
        pytest.approx(ret.final_molecules["180"].measure([2, 0, 1, 5]), abs=1.0e-2) == 180.0
        or pytest.approx(ret.final_molecules["180"].measure([2, 0, 1, 5]), abs=1.0e-2) == -180.0
    )
    assert pytest.approx(ret.final_molecules["0"].measure([2, 0, 1, 5]), abs=1.0e-2) == 0.0

    assert ret.provenance.creator.lower() == "torsiondrive"

    if from_v2(request.node.name):
        # properly "to_v1" should be always `== 4` but data lost in v2 by default
        if scan_ptcl == "none":
            assert not opthist_tgt
        else:
            assert len(opthist_tgt["180"]) == {"lowest": 1, "all": 4}[scan_ptcl]
    else:
        assert len(opthist_tgt["180"]) == 4

    if not (from_v2(request.node.name) and scan_ptcl == "none"):
        assert opthist_tgt["180"][0].provenance.creator.lower() == "geometric"
    if "v2" in request.node.name:
        if scan_ptcl != "none":
            assert opthist_tgt["180"][0].trajectory_results[0].provenance.creator.lower() == "rdkit"
            assert opthist_tgt["180"][0].trajectory_results[0].schema_version == 2
    else:
        if not ("to_v1" in request.node.name and scan_ptcl == "none"):
            assert opthist_tgt["180"][0].trajectory[0].provenance.creator.lower() == "rdkit"
            assert opthist_tgt["180"][0].trajectory[0].schema_version == 1

    assert ret.stdout == "All optimizations converged at lowest energy. Job Finished!\n"


@using("mace")
@using("torsiondrive")
def test_torsiondrive_extra_constraints(schema_versions, request):
    models, retver, _ = schema_versions

    keywords = {
        "coordsys": "dlc",
        "constraints": {
            "set": [
                {
                    "type": "dihedral",  # hold a dihedral through the other C-C bond fixed
                    "indices": (0, 1, 2, 10),
                    "value": 0.0,
                }
            ]
        },
    }

    if from_v2(request.node.name):
        input_data = models.TorsionDriveInput(
            initial_molecule=[models.Molecule(**qcng.get_molecule("propane", return_dict=True))],
            specification=models.TorsionDriveSpecification(
                keywords=models.TorsionDriveKeywords(dihedrals=[(3, 0, 1, 2)], grid_spacing=[180]),
                specification=models.OptimizationSpecification(
                    program="geomeTRIC",
                    keywords=keywords,
                    specification=models.AtomicSpecification(
                        # use mace as it does not have convergence issues like UFF
                        program="mace",
                        driver=models.DriverEnum.gradient,
                        model=models.Model(method="small", basis=None),
                    ),
                    protocols=models.OptimizationProtocols(trajectory_results="all"),
                ),
                protocols=models.TorsionDriveProtocols(scan_results="all"),
            ),
        )
    else:
        input_data = models.TorsionDriveInput(
            keywords=models.TDKeywords(dihedrals=[(3, 0, 1, 2)], grid_spacing=[180]),
            input_specification=models.QCInputSpecification(
                driver=models.DriverEnum.gradient, model=models.Model(method="small", basis=None)
            ),
            initial_molecule=[models.Molecule(**qcng.get_molecule("propane", return_dict=True))],
            optimization_spec=models.OptimizationSpecification(
                procedure="geomeTRIC",
                keywords={
                    **keywords,
                    # use mace as it does not have convergence issues like UFF
                    "program": "mace",
                },
            ),
        )

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "torsiondrive", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    if "_v2" not in request.node.name:
        assert ret.error is None
    assert ret.success

    expected_grid_ids = {"180", "0"}
    opthist_tgt = ret.scan_results if "v2" in request.node.name else ret.optimization_history

    assert {*opthist_tgt} == expected_grid_ids
    if "v2" in request.node.name:
        assert {*ret.scan_properties} == expected_grid_ids
    else:
        assert {*ret.final_energies} == expected_grid_ids
    assert {*ret.final_molecules} == expected_grid_ids
    if "v2" in request.node.name:
        assert ret.properties.calcinfo_ngrid == 2

    assert (
        pytest.approx(ret.final_molecules["180"].measure([3, 0, 1, 2]), abs=1.0e-2) == 180.0
        or pytest.approx(ret.final_molecules["180"].measure([3, 0, 1, 2]), abs=1.0e-2) == -180.0
    )
    assert pytest.approx(ret.final_molecules["180"].measure([0, 1, 2, 10]), abs=1.0e-2) == 0.0
    assert pytest.approx(ret.final_molecules["0"].measure([3, 0, 1, 2]), abs=1.0e-2) == 0.0

    assert ret.provenance.creator.lower() == "torsiondrive"
    assert opthist_tgt["180"][0].provenance.creator.lower() == "geometric"
    if "v2" in request.node.name:
        assert opthist_tgt["180"][0].trajectory_results[0].provenance.creator.lower() == "mace"
    else:
        assert opthist_tgt["180"][0].trajectory[0].provenance.creator.lower() == "mace"

    assert "Using MACE-OFF23 MODEL for MACECalculator" in ret.stdout
    assert "All optimizations converged at lowest energy. Job Finished!\n" in ret.stdout


@using("mrchem")
@pytest.mark.parametrize(
    "optimizer",
    [
        pytest.param("geometric", marks=using("geometric")),
        pytest.param("optking", marks=using("optking")),
        pytest.param("berny", marks=using("berny")),
    ],
)
def test_optimization_mrchem(input_data, optimizer, schema_versions, request):
    models, retver, _ = schema_versions

    input_data["initial_molecule"] = models.Molecule(**qcng.get_molecule("hydrogen", return_dict=True))
    if from_v2(request.node.name):
        input_data["specification"]["specification"]["model"] = {"method": "HF"}
        input_data["specification"]["specification"]["keywords"] = {"world_prec": 1.0e-4}
        input_data["specification"]["specification"]["program"] = "mrchem"
        input_data["specification"]["protocols"] = {"trajectory_results": "final"}  # to test provenance
    else:
        input_data["input_specification"]["model"] = {"method": "HF"}
        input_data["input_specification"]["keywords"] = {"world_prec": 1.0e-4}
        input_data["keywords"]["program"] = "mrchem"

    input_data = models.OptimizationInput(**input_data)

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, optimizer, raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    trajs_tgt = ret.trajectory_results if "v2" in request.node.name else ret.trajectory  # for looking at grad jobs
    props_tgt = ret.trajectory_properties if "v2" in request.node.name else ret.energies  # for counting opt iterations

    assert 10 > len(props_tgt) > 1
    if "v2" in request.node.name:
        assert 10 > ret.properties.optimization_iterations > 1

    assert pytest.approx(ret.final_molecule.measure([0, 1]), 1.0e-3) == 1.3860734486984705
    assert ret.provenance.creator.lower() == optimizer
    assert trajs_tgt[0].provenance.creator.lower() == "mrchem"
