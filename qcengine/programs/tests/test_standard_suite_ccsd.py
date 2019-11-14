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
        pytest.param("cfour", "aug-pvdz", {"SCF_CONV": 12, "CC_CONV": 12}, marks=testing.using_cfour),
        pytest.param("cfour", "aug-pvdz", {}, marks=testing.using_cfour),
        pytest.param("gamess", "accd", {"ccinp__ncore": 0, "contrl__ispher": 1}, marks=testing.using_gamess),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True}, marks=testing.using_nwchem),
        pytest.param(
            "nwchem", "aug-cc-pvdz", {"basis__spherical": True, "qc_module": "tce"}, marks=testing.using_nwchem
        ),
        pytest.param("psi4", "aug-cc-pvdz", {}, marks=testing.using_psi4),
        # pytest.param("qchem", "aug-cc-pvdz", {"N_FROZEN_CORE": 0}, marks=testing.using_qchem),
        # TODO Molpro has frozen-core on by default. For this to pass need new keyword frozen_core = False
        # pytest.param('molpro', 'aug-cc-pvdz', {}, marks=testing.using_molpro),
    ],
)
def test_sp_ccsd_rhf_full(program, basis, keywords, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD/qz2p on water

    """
    resi = {"molecule": h2o, "driver": "energy", "model": {"method": "ccsd", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    # scftot = -76.0413815332
    # mp2tot = -76.2632792578
    # mp2corl = -0.2218977246
    # ccsdcorl = -0.2294105794
    ccsdtot = -76.2707921127

    atol = 1.0e-6
    # assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
    # if not (mtd == 'nwc-ccsd' and keywords.get('qc_module', 'nein').lower() == 'tce'):
    #    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    # assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD corl', atol=atol)
    assert compare_values(ccsdtot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords,errmsg",
    [
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"ccsd__freeze": 1, "scf__uhf": True},
            "ccsd: nopen is not zero",
            marks=testing.using_nwchem,
        ),
        pytest.param(
            "gamess",
            "accd",
            {"contrl__scftyp": "uhf"},
            "CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF",
            marks=testing.using_gamess,
        ),
    ],
)
def test_sp_ccsd_uhf_fc_error(program, basis, keywords, nh2, errmsg):
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "ccsd", "basis": basis}, "keywords": keywords}

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert errmsg in str(e.value)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param(
            "cfour",
            "AUG-PVDZ",
            {"refERENCE": "ROhf", "OCCUPATION": [[3, 1, 1, 0], [3, 0, 1, 0]], "SCF_CONV": 12, "CC_CONV": 12},
            marks=testing.using_cfour,
        ),
        pytest.param("cfour", "AUG-PVDZ", {"REFERENCE": "ROHF"}, marks=testing.using_cfour),
        pytest.param(
            "gamess",
            "ACCD",
            {"CONTRL__ISPHER": 1, "CONTRL__SCFTYP": "ROHF", "CCINP__NCORE": 0},
            marks=testing.using_gamess,
        ),
        pytest.param(
            "nwchem",
            "AUG-CC-PVDZ",
            {"BASIS__SPHERICAL": True, "QC_MODULE": "TCE", "SCF__ROHF": True},
            marks=testing.using_nwchem,
        ),
        pytest.param("psi4", "AUG-CC-PVDZ", {"REFERENCE": "ROHF"}, marks=testing.using_psi4),
        # pytest.param("qchem", "AUG-CC-PVDZ", {"UNRESTRICTED": False}, marks=testing.using_qchem),
        # TODO Molpro has frozen-core on by default. For this to pass need new keyword frozen_core = False
        # pytest.param('molpro', 'aug-cc-pvdz', {}, marks=testing.using_molpro),
    ],
)
def test_sp_ccsd_rohf_full(program, basis, keywords, nh2):
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "ccsd", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    # scftot = -55.570724348574
    # ssccsdcorl = -0.0339827
    # osccsdcorl = -0.1442533
    # ccsdcorl = -0.178236032911
    ccsdtot = -55.748960381485

    atol = 1.0e-6
    # assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + 'SCF', atol=atol)
    # if not (method in ['gms-ccsd', 'nwc-ccsd']):
    #    # cfour isn't splitting out the singles from OS. and psi4 isn't incl the singles in either so SS + OS != corl
    #    # and maybe singles moved from SS to OS in Cfour btwn 2010 and 2014 versions (change in ref)
    #    # assert compare_values(osccsdcorl, qcdb.variable('ccsd opposite-spin correlation energy'), tnm() + ' CCSD OS corl', atol=atol)
    #    assert compare_values(ssccsdcorl, qcdb.variable('ccsd same-spin correlation energy'), tnm() + ' CCSD SS corl', atol=atol)
    # assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD corl', atol=atol)
    # assert compare_values(ccsdtot, qcdb.variable('ccsd total energy'), tnm() + ' CCSD', atol=atol)
    # assert compare_values(ccsdcorl, qcdb.variable('current correlation energy'), tnm() + ' Current corl', atol=atol)
    assert compare_values(ccsdtot, res["return_result"], atol=atol)
