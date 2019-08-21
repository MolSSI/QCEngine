import copy

import pytest

import qcelemental as qcel
from qcelemental.testing import compare_values
import qcengine as qcng
from qcengine import testing


_base_json = {"schema_name": "qcschema_input", "schema_version": 1}


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
units au
"""
    return qcel.models.Molecule.from_data(smol)


@pytest.mark.parametrize('prog,bas,opts', [
    #pytest.param('cfour', 'aug-pvdz', {'cfour_SCF_CONV': 12, 'cfour_CC_CONV': 12}, marks=testing.using_cfour),
    #pytest.param('cfour', {}, marks=testing.using_cfour),
    #pytest.param('nwchem', {}, marks=testing.using_nwchem),
    #pytest.param('nwchem', {'qc_module': 'tce'}, marks=testing.using_nwchem),
    pytest.param('psi4', 'aug-cc-pvdz', {}, marks=testing.using_psi4),
    pytest.param('gamess', 'accd', {'ccinp__ncore': 0, 'contrl__ispher': 1}, marks=testing.using_gamess),
])
def test_sp_ccsd_rhf_full(prog, bas, opts, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD/qz2p on water

    """
    resi = copy.deepcopy(_base_json)
    resi["molecule"] = h2o
    resi["driver"] = "energy"
    resi["model"] = {"method": "ccsd", "basis": bas}
    resi["keywords"] = opts

    res = qcng.compute(resi, prog, raise_error=True, return_dict=True)
#    import pprint
#    pprint.pprint(res)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    scftot = -76.0413815332
    mp2tot = -76.2632792578
    mp2corl = -0.2218977246
    ccsdcorl = -0.2294105794
    ccsdtot = -76.2707921127

    atol = 1.e-6
    #assert compare_values(scftot, qcdb.variable('scf total energy'), tnm() + ' SCF', atol=atol)
    #if not (mtd == 'nwc-ccsd' and opts.get('qc_module', 'nein').lower() == 'tce'):
    #    assert compare_values(mp2tot, qcdb.variable('mp2 total energy'), tnm() + ' MP2', atol=atol)
    #assert compare_values(ccsdcorl, qcdb.variable('ccsd correlation energy'), tnm() + ' CCSD corl', atol=atol)
    assert compare_values(ccsdtot, res["return_result"], atol=atol)

