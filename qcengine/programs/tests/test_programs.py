"""
Tests the DQM compute dispatch module
"""


import numpy as np
import pytest
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import checkver_and_convert, failure_engine, from_v2, schema_versions, using


def test_missing_key(schema_versions, request):
    _, retver, _ = schema_versions

    ret = qcng.compute({"hello": "hi"}, "bleh", return_version=retver)
    if "as_v2" not in request.node.name:
        # input provides no clue to v2/v1, so checkver fails
        ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    assert ret.success is False
    assert "hello" in ret.input_data


def test_missing_key_raises(schema_versions, request):
    with pytest.raises(qcng.exceptions.InputError):
        ret = qcng.compute({"hello": "hi"}, "bleh", raise_error=True)


@using("psi4")
def test_psi4_task(schema_versions, request):
    models, retver, _ = schema_versions

    if from_v2(request.node.name):
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "specification": {
                "driver": "energy",
                "model": {"method": "SCF", "basis": "sto-3g"},
                "keywords": {"scf_type": "df"},
            },
        }
    else:
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "driver": "energy",
            "model": {"method": "SCF", "basis": "sto-3g"},
            "keywords": {"scf_type": "df"},
        }

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "psi4", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    if "v2" in request.node.name:
        assert ret.input_data.specification.driver == "energy"
    else:
        assert ret.driver == "energy"
    assert "Final Energy" in ret.stdout

    prov_keys = {"cpu", "hostname", "username", "wall_time"}
    assert ret.provenance.model_dump().keys() >= prov_keys
    assert "retries" not in ret.provenance.model_dump()

    assert ret.success is True


@using("psi4")
@using("gcp")
def test_psi4_hf3c_task(schema_versions, request):
    models, retver, _ = schema_versions

    if from_v2(request.node.name):
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "specification": {
                "driver": "energy",
                "model": {"method": "HF3c"},
                "keywords": {"scf_type": "df"},
            },
        }
    else:
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "driver": "energy",
            "model": {"method": "HF3c"},
            "keywords": {"scf_type": "df"},
        }

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "psi4", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success is True
    if "v2" in request.node.name:
        # prior to 0.50, None, now ""
        assert not ret.input_data.specification.model.basis
    else:
        assert not ret.model.basis


@using("psi4_runqcsk")
def test_psi4_interactive_task(schema_versions, request):
    models, retver, _ = schema_versions

    if from_v2(request.node.name):
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "specification": {
                "driver": "energy",
                "model": {"method": "SCF", "basis": "sto-3g"},
                "keywords": {"scf_type": "df"},
                "extras": {"psiapi": True},
            },
        }
    else:
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "driver": "energy",
            "model": {"method": "SCF", "basis": "sto-3g"},
            "keywords": {"scf_type": "df"},
            "extras": {"psiapi": True},
        }

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "psi4", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert "Final Energy" in ret.stdout
    assert ret.success
    is_psiapi_evaluated = ret.extras.pop("psiapi_evaluated", False)
    assert is_psiapi_evaluated


@using("psi4_runqcsk")
def test_psi4_wavefunction_task(schema_versions, request):
    models, retver, _ = schema_versions

    if from_v2(request.node.name):
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "specification": {
                "driver": "energy",
                "model": {"method": "SCF", "basis": "sto-3g"},
                "protocols": {"wavefunction": "orbitals_and_eigenvalues"},
                "keywords": {"scf_type": "df"},
            },
        }
    else:
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "driver": "energy",
            "model": {"method": "SCF", "basis": "sto-3g"},
            "keywords": {"scf_type": "df"},
            "protocols": {"wavefunction": "orbitals_and_eigenvalues"},
        }

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "psi4", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success, ret.error.error_message
    assert ret.wavefunction.scf_orbitals_a.shape == (7, 7)


@using("psi4")
def test_psi4_internal_failure(schema_versions, request):
    models, retver, _ = schema_versions

    mol = models.Molecule.from_data(
        """0 3
     O    0.000000000000     0.000000000000    -0.068516245955
    """
    )

    if from_v2(request.node.name):
        psi4_task = {
            "molecule": mol,
            "specification": {
                "driver": "energy",
                "model": {"method": "ccsd", "basis": "6-31g"},
                "keywords": {"reference": "rhf"},
            },
        }
    else:
        psi4_task = {
            "molecule": mol,
            "driver": "energy",
            "model": {"method": "ccsd", "basis": "6-31g"},
            "keywords": {"reference": "rhf"},
        }
    with pytest.raises(qcng.exceptions.InputError) as exc:
        psi4_task = checkver_and_convert(psi4_task, request.node.name, "pre")
        ret = qcng.compute(psi4_task, "psi4", raise_error=True, return_version=retver)

    assert "reference is only" in str(exc.value)


@using("psi4")
def test_psi4_ref_switch(schema_versions, request):
    models, retver, _ = schema_versions

    if from_v2(request.node.name):
        inp = models.AtomicInput(
            **{
                "molecule": {"symbols": ["Li"], "geometry": [0, 0, 0], "molecular_multiplicity": 2},
                "specification": {
                    "driver": "energy",
                    "model": {"method": "B3LYP", "basis": "sto-3g"},
                    "keywords": {"scf_type": "df"},
                },
            }
        )
    else:
        inp = models.AtomicInput(
            **{
                "molecule": {"symbols": ["Li"], "geometry": [0, 0, 0], "molecular_multiplicity": 2},
                "driver": "energy",
                "model": {"method": "B3LYP", "basis": "sto-3g"},
                "keywords": {"scf_type": "df"},
            }
        )

    inp = checkver_and_convert(inp, request.node.name, "pre")
    ret = qcng.compute(inp, "psi4", raise_error=True, return_dict=False, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success is True
    assert ret.properties.calcinfo_nalpha == 2
    assert ret.properties.calcinfo_nbeta == 1


@using("rdkit")
@pytest.mark.parametrize("method", ["UFF", "MMFF94", "MMFF94s"])
def test_rdkit_task(method, schema_versions, request):
    models, retver, _ = schema_versions

    if from_v2(request.node.name):
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "specification": {"driver": "gradient", "model": {"method": method}, "keywords": {}},
        }
    else:
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "driver": "gradient",
            "model": {"method": method},
            "keywords": {},
        }

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "rdkit", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success is True


@using("rdkit")
def test_rdkit_connectivity_error(schema_versions, request):
    models, retver, _ = schema_versions

    if from_v2(request.node.name):
        input_data = {
            "molecule": qcng.get_molecule("water", return_dict=True),
            "specification": {"driver": "energy", "model": {"method": "UFF", "basis": ""}, "keywords": {}},
        }
    else:
        input_data = {
            "molecule": qcng.get_molecule("water", return_dict=True),
            "driver": "energy",
            "model": {"method": "UFF", "basis": ""},
            "keywords": {},
        }
    del input_data["molecule"]["connectivity"]

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "rdkit", return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    assert ret.success is False
    assert "connectivity" in ret.error.error_message

    with pytest.raises(qcng.exceptions.InputError):
        qcng.compute(input_data, "rdkit", raise_error=True, return_version=retver)


@using("torchani")
def test_torchani_task(schema_versions, request):
    models, retver, _ = schema_versions

    if from_v2(request.node.name):
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "specification": {"driver": "gradient", "model": {"method": "ANI1x", "basis": None}, "keywords": {}},
        }
    else:
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "driver": "gradient",
            "model": {"method": "ANI1x", "basis": None},
            "keywords": {},
        }

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "torchani", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success is True
    if "v2" in request.node.name:
        assert ret.input_data.specification.driver == "gradient"
    else:
        assert ret.driver == "gradient"


@using("mopac")
def test_mopac_task(schema_versions, request):
    _, retver, _ = schema_versions

    if from_v2(request.node.name):
        input_data = {
            "molecule": qcng.get_molecule("water", return_dict=True),
            "specification": {
                "driver": "gradient",
                "model": {"method": "PM6", "basis": None},
                "keywords": {"pulay": False},
            },
        }
    else:
        input_data = {
            "molecule": qcng.get_molecule("water", return_dict=True),
            "driver": "gradient",
            "model": {"method": "PM6", "basis": None},
            "keywords": {"pulay": False},
        }

    input_data1 = checkver_and_convert(input_data.copy(), request.node.name, "pre")
    ret = qcng.compute(input_data1, "mopac", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.extras.keys() >= {"heat_of_formation", "dip_vec"}
    energy = pytest.approx(-0.08474117913025125, rel=1.0e-5)

    # Check gradient
    input_data2 = checkver_and_convert(input_data.copy(), request.node.name, "pre")
    ret = qcng.compute(input_data2, "mopac", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.extras.keys() >= {"heat_of_formation", "dip_vec"}
    assert np.linalg.norm(ret.return_result) == pytest.approx(0.03543560156912385, rel=3.0e-4)
    assert ret.properties.return_energy == energy

    # Check energy
    input_data3 = input_data.copy()
    if from_v2(request.node.name):
        input_data3["specification"]["driver"] = "energy"
    else:
        input_data3["driver"] = "energy"
    input_data3 = checkver_and_convert(input_data3, request.node.name, "pre")
    ret = qcng.compute(input_data3, "mopac", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.return_result == energy
    assert "== MOPAC DONE ==" in ret.stdout


def test_random_failure_no_retries_input(failure_engine, schema_versions, request):
    _, retver, _ = schema_versions

    failure_engine.iter_modes = ["input_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    assert ret.error.error_type == "input_error"
    provenance_keys = (
        ret.input_data.provenance.model_dump().keys()
        if ("_v2" in request.node.name)
        else ret.input_data["provenance"].keys()
    )
    assert "retries" not in provenance_keys


def test_random_failure_no_retries_random(failure_engine, schema_versions, request):
    _, retver, _ = schema_versions

    failure_engine.iter_modes = ["random_error"]
    ret = qcng.compute(failure_engine.get_job(), failure_engine.name, raise_error=False, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    assert ret.error.error_type == "random_error"
    provenance_keys = (
        ret.input_data.provenance.model_dump().keys()
        if ("_v2" in request.node.name)
        else ret.input_data["provenance"].keys()
    )
    assert "retries" not in provenance_keys


def test_random_failure_with_retries(failure_engine, schema_versions, request):
    _, retver, _ = schema_versions

    failure_engine.iter_modes = ["random_error", "random_error", "random_error"]
    ret = qcng.compute(
        failure_engine.get_job(),
        failure_engine.name,
        raise_error=False,
        task_config={"retries": 2},
        return_version=retver,
    )
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    retries = (
        ret.input_data.provenance.retries if ("_v2" in request.node.name) else ret.input_data["provenance"]["retries"]
    )
    assert retries == 2
    assert ret.error.error_type == "random_error"

    failure_engine.iter_modes = ["random_error", "input_error"]
    ret = qcng.compute(
        failure_engine.get_job(),
        failure_engine.name,
        raise_error=False,
        task_config={"retries": 4},
        return_version=retver,
    )
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    retries = (
        ret.input_data.provenance.retries if ("_v2" in request.node.name) else ret.input_data["provenance"]["retries"]
    )
    assert retries == 1
    assert ret.error.error_type == "input_error"


def test_random_failure_with_success(failure_engine, schema_versions, request):
    _, retver, _ = schema_versions

    failure_engine.iter_modes = ["random_error", "pass"]
    failure_engine.ncalls = 0
    ret = qcng.compute(
        failure_engine.get_job(),
        failure_engine.name,
        raise_error=False,
        task_config={"retries": 1},
        return_version=retver,
    )
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success, ret.error.error_message
    assert ret.provenance.retries == 1
    assert ret.extras["ncalls"] == 2


@using("openmm")
def test_openmm_task_smirnoff(schema_versions, request):
    models, retver, _ = schema_versions

    from qcengine.programs.openmm import OpenMMHarness

    if from_v2(request.node.name):
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "specification": {
                "driver": "energy",
                "model": {"method": "openff-1.0.0", "basis": "smirnoff"},
                "keywords": {},
            },
        }
    else:
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "driver": "energy",
            "model": {"method": "openff-1.0.0", "basis": "smirnoff"},
            "keywords": {},
        }

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "openmm", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    cachelength = len(OpenMMHarness._CACHE)

    assert cachelength > 0
    assert ret.success is True

    ret = qcng.compute(input_data, "openmm", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    # ensure cache has not grown
    assert len(OpenMMHarness._CACHE) == cachelength
    assert ret.success is True


@pytest.mark.skip("`basis` must be explicitly specified at this time")
@using("openmm")
def test_openmm_task_url_basis(schema_versions, request):
    models, retver, _ = schema_versions

    from qcengine.programs.openmm import OpenMMHarness

    if from_v2(request.node.name):
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "specification": {
                "driver": "energy",
                "model": {
                    "method": "openmm",
                    "basis": "openff-1.0.0",
                    "url": "https://raw.githubusercontent.com/openforcefield/openff-forcefields/1.0.0/openforcefields/offxml/openff-1.0.0.offxml",
                },
                "keywords": {},
            },
        }
    else:
        input_data = {
            "molecule": models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            "driver": "energy",
            "model": {
                "method": "openmm",
                "basis": "openff-1.0.0",
                "url": "https://raw.githubusercontent.com/openforcefield/openff-forcefields/1.0.0/openforcefields/offxml/openff-1.0.0.offxml",
            },
            "keywords": {},
        }

    input_data = checkver_and_convert(input_data, request.node.name, "pre")
    ret = qcng.compute(input_data, "openmm", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    cachelength = len(OpenMMHarness._CACHE)

    assert cachelength > 0
    assert ret.success is True

    ret = qcng.compute(input_data, "openmm", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    # ensure cache has not grown
    assert len(OpenMMHarness._CACHE) == cachelength
    assert ret.success is True


@using("openmm")
def test_openmm_cmiles_gradient(schema_versions, request):
    models, retver, _ = schema_versions

    program = "openmm"

    water = models.Molecule(**qcng.get_molecule("water", return_dict=True))

    water_dict = water.model_dump()
    # add water cmiles to the molecule
    water_dict["extras"] = {"cmiles": {"canonical_isomeric_explicit_hydrogen_mapped_smiles": "[H:2][O:1][H:3]"}}

    molecule = models.Molecule.from_data(water_dict)

    model = {"method": "openff-1.0.0", "basis": "smirnoff"}

    if from_v2(request.node.name):
        inp = models.AtomicInput(molecule=molecule, specification={"driver": "gradient", "model": model})
    else:
        inp = models.AtomicInput(molecule=molecule, driver="gradient", model=model)

    inp = checkver_and_convert(inp, request.node.name, "pre")
    ret = qcng.compute(inp, program, raise_error=False, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert ret.success is True


@using("openmm")
def test_openmm_cmiles_gradient_nomatch(schema_versions, request):
    models, retver, _ = schema_versions

    program = "openmm"

    water = models.Molecule(**qcng.get_molecule("water", return_dict=True))

    water_dict = water.model_dump()
    # add ethane cmiles to the molecule
    water_dict["extras"] = {
        "cmiles": {
            "canonical_isomeric_explicit_hydrogen_mapped_smiles": "[H:3][C:1]([H:4])([H:5])[C:2]([H:6])([H:7])[H:8]"
        }
    }

    molecule = models.Molecule.from_data(water_dict)

    model = {"method": "openff-1.0.0", "basis": "smirnoff"}

    if from_v2(request.node.name):
        inp = models.AtomicInput(molecule=molecule, specification={"driver": "gradient", "model": model})
    else:
        inp = models.AtomicInput(molecule=molecule, driver="gradient", model=model)

    inp = checkver_and_convert(inp, request.node.name, "pre")
    ret = qcng.compute(inp, program, raise_error=False, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post", vercheck=False)

    # if we correctly find the cmiles this should fail as the molecule and cmiles are different
    assert ret.success is False
    assert (
        "molecule.add_conformer given input of the wrong shape: Given (3, 3), expected (8, 3)"
        in ret.error.error_message
    )


@using("openmm")
@pytest.mark.parametrize(
    "gaff_settings",
    [
        pytest.param(({}, None, 0.0013904199062156914), id="gaff no keywords"),
        pytest.param(({"constraints": "ALLBONDS"}, None, 8.108238580315493e-05), id="constraints allbonds"),
        pytest.param(({"nonbondedMethod": "LjPmE"}, None, 0.0013904199062156914), id="nonbonded ljpme"),
        pytest.param(
            ({"nonbondedMethod": "PME", "constraints": "Hbonds"}, None, 8.108238580315493e-05),
            id="nonbonded pme constraints hbonds",
        ),
        pytest.param(({"constraints": "badmethod"}, ValueError, 0), id="incorrect constraint"),
        pytest.param(({"nonbondedMethod": "badmethod"}, ValueError, 0), id="incorrect nonbondedmethod"),
    ],
)
def test_openmm_gaff_keywords(gaff_settings, schema_versions, request):
    """
    Test the different running settings with gaff.
    """
    models, retver, _ = schema_versions

    program = "openmm"
    water = models.Molecule(**qcng.get_molecule("water", return_dict=True))

    water_dict = water.model_dump()
    # add water cmiles to the molecule
    water_dict["extras"] = {"cmiles": {"canonical_isomeric_explicit_hydrogen_mapped_smiles": "[H:2][O:1][H:3]"}}

    molecule = models.Molecule.from_data(water_dict)
    keywords, error, expected_result = gaff_settings
    model = {"method": "gaff-2.1", "basis": "antechamber"}
    if from_v2(request.node.name):
        inp = models.AtomicInput(
            molecule=molecule, specification={"driver": "energy", "model": model, "keywords": keywords}
        )
    else:
        inp = models.AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)
    if error is not None:
        with pytest.raises(error):
            inp = checkver_and_convert(inp, request.node.name, "pre")
            _ = qcng.compute(inp, program, raise_error=True, return_version=retver)
    else:
        inp = checkver_and_convert(inp, request.node.name, "pre")
        ret = qcng.compute(inp, program, raise_error=False, return_version=retver)
        ret = checkver_and_convert(ret, request.node.name, "post")

        assert ret.success is True
        assert ret.return_result == pytest.approx(expected_result, rel=1e-6)


@using("mace")
def test_mace_energy(schema_versions, request):
    """
    Test calculating the energy with mace
    """
    models, retver, _ = schema_versions
    water = models.Molecule(**qcng.get_molecule("water", return_dict=True))
    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule=water, specification={"driver": "energy", "model": {"method": "small", "basis": None}}
        )
    else:
        atomic_input = models.AtomicInput(molecule=water, model={"method": "small", "basis": None}, driver="energy")

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    result = qcng.compute(atomic_input, "mace", return_version=retver)
    result = checkver_and_convert(result, request.node.name, "post")

    assert result.success
    assert pytest.approx(result.return_result) == -76.47683956098838


@using("mace")
def test_mace_gradient(schema_versions, request):
    """
    Test calculating the gradient with mace
    """
    models, retver, _ = schema_versions

    water = models.Molecule(**qcng.get_molecule("water", return_dict=True))
    expected_result = np.array(
        [
            [0.0, -2.1590400539385646e-18, -0.04178551770271103],
            [0.0, -0.029712483642769006, 0.020892758851355515],
            [0.0, 0.029712483642769006, 0.020892758851355518],
        ]
    )

    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule=water, specification={"driver": "gradient", "model": {"method": "small", "basis": None}}
        )
    else:
        atomic_input = models.AtomicInput(molecule=water, model={"method": "small", "basis": None}, driver="gradient")

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    result = qcng.compute(atomic_input, "mace", return_version=retver)
    result = checkver_and_convert(result, request.node.name, "post")

    assert result.success
    assert pytest.approx(result.return_result) == expected_result


@using("aimnet2")
@pytest.mark.parametrize(
    "model, expected_energy",
    [
        pytest.param("b973c", -76.39604306960972, id="b973c"),
        pytest.param("wb97m-d3", -76.47412023758551, id="wb97m-d3"),
    ],
)
def test_aimnet2_energy(model, expected_energy, schema_versions, request):
    """Test computing the energies of water with two aimnet2 models."""
    models, retver, _ = schema_versions

    water = models.Molecule(**qcng.get_molecule("water", return_dict=True))
    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule=water, specification={"driver": "energy", "model": {"method": model, "basis": None}}
        )
    else:
        atomic_input = models.AtomicInput(molecule=water, model={"method": model, "basis": None}, driver="energy")

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    result = qcng.compute(atomic_input, "aimnet2", return_version=retver)
    result = checkver_and_convert(result, request.node.name, "post")

    assert result.success
    assert pytest.approx(result.return_result) == expected_energy
    assert "charges" in result.extras["aimnet2"]
    assert "ensemble_charges_std" in result.extras["aimnet2"]
    assert "ensemble_forces_std" in result.extras["aimnet2"]


@using("aimnet2")
def test_aimnet2_gradient(schema_versions, request):
    """Test computing the gradient of water using one aimnet2 model."""
    models, retver, _ = schema_versions

    water = models.Molecule(**qcng.get_molecule("water", return_dict=True))
    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule=water, specification={"driver": "gradient", "model": {"method": "wb97m-d3", "basis": None}}
        )
    else:
        atomic_input = models.AtomicInput(
            molecule=water, model={"method": "wb97m-d3", "basis": None}, driver="gradient"
        )

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    result = qcng.compute(atomic_input, "aimnet2", return_version=retver)
    result = checkver_and_convert(result, request.node.name, "post")

    assert result.success
    # make sure the gradient is now the return result
    assert np.allclose(
        result.return_result,
        np.array(
            [
                [-0.0, 2.6080331227973375e-09, -0.04097248986363411],
                [-0.0, -0.029529934749007225, 0.020486244931817055],
                [-0.0, 0.029529931023716927, 0.020486244931817055],
            ]
        ),
    )
    assert pytest.approx(result.properties.return_energy) == -76.47412023758551
    # make sure the other properties were also saved
    assert "charges" in result.extras["aimnet2"]


@using("psi4")
def test_psi4_properties_driver(schema_versions, request):
    models, retver, _ = schema_versions

    import numpy as np

    json_data = {
        "schema_name": "qcschema_input",
        "schema_version": 1,
        "molecule": {
            "geometry": [
                0.0,
                0.0,
                -0.1294769411935893,
                0.0,
                -1.494187339479985,
                1.0274465079245698,
                0.0,
                1.494187339479985,
                1.0274465079245698,
            ],
            "symbols": ["O", "H", "H"],
            "fix_com": True,
            "fix_orientation": True,
        },
        "driver": "properties",
        "model": {
            "method": "HF",
            "basis": "6-31G",
            "properties": [
                "dipole",
                "quadrupole",
                "mulliken_charges",
                "lowdin_charges",
                "wiberg_lowdin_indices",
                "mayer_indices",
            ],
        },
        "keywords": {"scf_type": "df", "mp2_type": "df", "e_convergence": 9},
    }
    if from_v2(request.node.name):
        json_data["schema_name"] = "qcschema_atomic_input"
        json_data["schema_version"] = 2
        json_data["specification"] = {
            "driver": json_data.pop("driver"),
            "model": json_data.pop("model"),
            "keywords": json_data.pop("keywords"),
        }

    # Write expected output (dipole & quadrupole in au)
    expected_return_result = {
        "dipole": np.array([0.0, 0.0, 1.04037263345]).reshape((3,)),
        "quadrupole": np.array([-5.42788218076, 0.0, 0.0, 0.0, -3.07521129293, 0.0, 0.0, 0.0, -4.36605314966]).reshape(
            (3, 3)
        ),
        "mulliken_charges": [-0.7967275695997689, 0.3983637847998823, 0.3983637847998822],
        "lowdin_charges": [-0.5945105406840803, 0.29725527034203636, 0.29725527034203636],
        "wiberg_lowdin_indices": np.array(
            [
                0.0,
                0.9237385044125344,
                0.9237385044125329,
                0.9237385044125344,
                0.0,
                0.006992650019441531,
                0.9237385044125329,
                0.006992650019441531,
                0.0,
            ]
        ).reshape((3, 3)),
        "mayer_indices": np.array(
            [
                0.0,
                0.802064044935596,
                0.8020640449355959,
                0.802064044935596,
                0.0,
                0.003020025778524219,
                0.8020640449355959,
                0.003020025778524219,
                0.0,
            ]
        ).reshape((3, 3)),
    }

    expected_properties = {
        "calcinfo_nbasis": 13,
        "calcinfo_nmo": 13,
        "calcinfo_nalpha": 5,
        "calcinfo_nbeta": 5,
        "calcinfo_natom": 3,
        "scf_one_electron_energy": -122.27509111304202,
        "scf_two_electron_energy": 37.49348718008625,
        "nuclear_repulsion_energy": 8.80146206062943,
        "scf_total_energy": -75.98014187232634,
        "return_energy": -75.98014187232634,
    }

    json_data = checkver_and_convert(json_data, request.node.name, "pre")
    json_ret = qcng.compute(json_data, "psi4", return_version=retver)
    json_ret = checkver_and_convert(json_ret, request.node.name, "post")

    assert json_ret.success
    xptd_schema_name = "qcschema_atomic_result" if "v2" in request.node.name else "qcschema_output"
    assert json_ret.schema_name == xptd_schema_name
    for k in expected_return_result.keys():
        assert compare_values(expected_return_result[k], json_ret.return_result[k], atol=1.0e-5)
    for k in expected_properties.keys():
        assert compare_values(expected_properties[k], getattr(json_ret.properties, k), atol=1.0e-5)
