import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import checkver_and_convert, from_v2, schema_versions, using


@pytest.fixture
def h2o_data():
    return """
 # R=0.958 A=104.5
 H                  0.000000000000     1.431430901356     0.984293362719
 O                  0.000000000000     0.000000000000    -0.124038860300
 H                  0.000000000000    -1.431430901356     0.984293362719
 units au
"""


@pytest.fixture
def nh2_data():
    return """
 # R=1.008 #A=105.0
 0 2
 N   0.000000000000000   0.000000000000000  -0.145912918634892
 H   0.000000000000000  -1.511214298139000   1.013682596946108
 H   0.000000000000000   1.511214298139000   1.013682596946108
 units au
"""


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param("cfour", "aug-pvdz", {"scf_conv": 12, "cc_conv": 12}, marks=using("cfour")),
        pytest.param("cfour", "aug-pvdz", {}, marks=using("cfour")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True}, marks=using("nwchem")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "qc_module": "tce"}, marks=using("nwchem")),
        pytest.param("psi4", "aug-cc-pvdz", {}, marks=using("psi4")),
        pytest.param("gamess", "accd", {"ccinp__ncore": 0, "contrl__ispher": 1}, marks=using("gamess")),
    ],
)
def test_sp_ccsd_t_rhf_full(program, basis, keywords, h2o_data, schema_versions, request):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD(T)/adz on water

    """
    models, retver, _ = schema_versions
    h2o = models.Molecule.from_data(h2o_data)

    if from_v2(request.node.name):
        resi = {
            "molecule": h2o,
            "specification": {"driver": "energy", "model": {"method": "ccsd(t)", "basis": basis}, "keywords": keywords},
        }
    else:
        resi = {
            "molecule": h2o,
            "driver": "energy",
            "model": {"method": "ccsd(t)", "basis": basis},
            "keywords": keywords,
        }

    resi = checkver_and_convert(resi, request.node.name, "pre")
    res = qcng.compute(resi, program, raise_error=True, return_dict=True, return_version=retver)
    res = checkver_and_convert(res, request.node.name, "post")

    if "v2" in request.node.name:
        assert res["input_data"]["specification"]["driver"] == "energy"
    else:
        assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    ccsd_t_tot = -76.276030676767

    atol = 1.0e-6
    assert compare_values(ccsd_t_tot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords,errmsg",
    [
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"ccsd__freeze": 1, "scf__uhf": True},
            "ccsd: nopen is not zero",
            marks=using("nwchem"),
        ),
        pytest.param(
            "gamess",
            "accd",
            {"contrl__scftyp": "uhf"},
            "CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF",
            marks=using("gamess"),
        ),
    ],
)
def test_sp_ccsd_t_uhf_fc_error(program, basis, keywords, nh2_data, errmsg, schema_versions, request):
    models, retver, _ = schema_versions
    nh2 = models.Molecule.from_data(nh2_data)

    if from_v2(request.node.name):
        resi = {
            "molecule": nh2,
            "specification": {"driver": "energy", "model": {"method": "ccsd(t)", "basis": basis}, "keywords": keywords},
        }
    else:
        resi = {
            "molecule": nh2,
            "driver": "energy",
            "model": {"method": "ccsd(t)", "basis": basis},
            "keywords": keywords,
        }

    with pytest.raises(qcng.exceptions.InputError) as e:
        resi = checkver_and_convert(resi, request.node.name, "pre")
        qcng.compute(resi, program, raise_error=True, return_dict=True, return_version=retver)

    assert errmsg in str(e.value)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param(
            "cfour",
            "aug-pvdz",
            {"reference": "rohf", "occupation": [[3, 1, 1, 0], [3, 0, 1, 0]], "scf_conv": 12, "cc_conv": 12},
            marks=using("cfour"),
        ),
        pytest.param("cfour", "aug-pvdz", {"reference": "rohf"}, marks=using("cfour")),
        # pytest.param('nwchem', 'aug-cc-pvdz', {'basis__spherical': True, 'qc_module': 'tce', 'scf__rohf': True}, marks=using("nwchem")),
        pytest.param("psi4", "aug-cc-pvdz", {"reference": "rohf"}, marks=using("psi4")),
    ],
)
def test_sp_ccsd_t_rohf_full(program, basis, keywords, nh2_data, schema_versions, request):
    models, retver, _ = schema_versions
    nh2 = models.Molecule.from_data(nh2_data)

    if from_v2(request.node.name):
        resi = {
            "molecule": nh2,
            "specification": {"driver": "energy", "model": {"method": "ccsd(t)", "basis": basis}, "keywords": keywords},
        }
    else:
        resi = {
            "molecule": nh2,
            "driver": "energy",
            "model": {"method": "ccsd(t)", "basis": basis},
            "keywords": keywords,
        }

    resi = checkver_and_convert(resi, request.node.name, "pre")
    res = qcng.compute(resi, program, raise_error=True, return_version=retver)
    res = checkver_and_convert(res, request.node.name, "post")
    res = res.model_dump()

    if "v2" in request.node.name:
        assert res["input_data"]["specification"]["driver"] == "energy"
    else:
        assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    ccsd_t_tot = -55.752861467462

    atol = 1.0e-6
    assert compare_values(ccsd_t_tot, res["return_result"], atol=atol)
