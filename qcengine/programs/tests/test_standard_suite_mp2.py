import pytest

import qcelemental as qcel
import qcengine as qcng
from qcelemental.testing import compare_values
from qcengine import testing


@pytest.fixture
def h2o():
    smol = """
 # R=0.958 A=104.5
 H                  0.000000000000     1.431430901356     0.984293362719
 O                  0.000000000000     0.000000000000    -0.124038860300
 H                  0.000000000000    -1.431430901356     0.984293362719
 units au
"""
    return qcel.models.Molecule.from_data(smol)


@pytest.fixture
def nh2():
    smol = """
 # R=1.008 #A=105.0
 0 2
 N   0.000000000000000   0.000000000000000  -0.145912918634892
 H   0.000000000000000  -1.511214298139000   1.013682596946108
 H   0.000000000000000   1.511214298139000   1.013682596946108
 units au
"""
    return qcel.models.Molecule.from_data(smol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param("cfour", "aug-pvdz", {"scf_conv": 12}, marks=testing.using_cfour),
        pytest.param("cfour", "aug-pvdz", {}, marks=testing.using_cfour),
        pytest.param("gamess", "accd", {"mp2__nacore": 0, "contrl__ispher": 1}, marks=testing.using_gamess),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True}, marks=testing.using_nwchem),
        pytest.param(
            "nwchem", "aug-cc-pvdz", {"basis__spherical": True, "qc_module": "tce"}, marks=testing.using_nwchem
        ),
        pytest.param("psi4", "aug-cc-pvdz", {"mp2_type": "conv"}, marks=testing.using_psi4),
        pytest.param("qchem", "aug-cc-pvdz", {"N_FROZEN_CORE": 0}, marks=testing.using_qchem),
        # TODO Molpro has frozen-core on by default. For this to pass need keyword frozen_core = False
        # pytest.param('molpro', 'aug-cc-pvdz', {}, marks=testing.using_molpro),
    ],
)
def test_sp_mp2_rhf_full(program, basis, keywords, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point MP2/adz on water

    """
    resi = {"molecule": h2o, "driver": "energy", "model": {"method": "mp2", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    mp2_tot = -76.2632792578

    atol = 1.0e-6
    assert compare_values(mp2_tot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param(
            "cfour",
            "aug-pvdz",
            {
                "reference": "uhf",
                "occupation": [[3, 1, 1, 0], [3, 0, 1, 0]],
                "dropmo": [1],
                "scf_conv": 12,
                "cc_conv": 12,
            },
            marks=testing.using_cfour,
        ),
        pytest.param("cfour", "aug-pvdz", {"reference": "uhf", "dropmo": 1}, marks=testing.using_cfour),
        pytest.param("gamess", "accd", {"contrl__ispher": 1, "contrl__scftyp": "uhf"}, marks=testing.using_gamess),
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"basis__spherical": True, "qc_module": "tce", "scf__uhf": True, "tce__freeze": 1},
            marks=testing.using_nwchem,
        ),
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"basis__spherical": True, "scf__uhf": True, "mp2__freeze": 1},
            marks=testing.using_nwchem,
        ),
        pytest.param(
            "psi4",
            "aug-cc-pvdz",
            {"reference": "uhf", "freeze_core": True, "mp2_type": "conv"},
            marks=testing.using_psi4,
        ),
        pytest.param("qchem", "aug-cc-pvdz", {"N_frozen_CORE": "fC"}, marks=testing.using_qchem),
        # TODO Molpro needs a new keyword for unrestricted MP2 (otherwise RMP2 by default) and needs symmetry c1
        # pytest.param('molpro', 'aug-cc-pvdz', {"reference": "unrestricted"}, marks=testing.using_molpro),
    ],
)
def test_sp_mp2_uhf_fc(program, basis, keywords, nh2):
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "mp2", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    mp2_tot = -55.727565606601

    atol = 1.0e-6
    assert compare_values(mp2_tot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords,errmsg",
    [pytest.param("nwchem", "aug-cc-pvdz", {"scf__rohf": True}, "unknown SCFTYPE", marks=testing.using_nwchem)],
)
def test_sp_mp2_rohf_full_error(program, basis, keywords, nh2, errmsg):
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "mp2", "basis": basis}, "keywords": keywords}

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert errmsg in str(e.value)
