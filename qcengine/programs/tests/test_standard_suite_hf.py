import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import checkver_and_convert, from_v2, schema_versions, schema_versions5, using


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
 symmetry c1
"""


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param("cfour", "cC-pvdZ", {"scf_conv": 12}, marks=using("cfour")),  # test basis handling, not results
        pytest.param("cfour", "aug-pvdz", {"scf_conv": 12}, marks=using("cfour")),
        pytest.param("cfour", "aug-pvdz", {}, marks=using("cfour")),
        pytest.param(
            "qcore",
            "aug-cc-pVDZ",
            {"coulomb_method": "direct_4idx", "exchange_method": "direct_4idx"},
            marks=using("qcore"),
        ),
        pytest.param("gamess", "accd", {"contrl__ispher": 1}, marks=using("gamess")),
        pytest.param("molpro", "aug-cc-pvdz", {}, marks=using("molpro")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True}, marks=using("nwchem")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "qc_module": "tce"}, marks=using("nwchem")),
        pytest.param("psi4", "aug-cc-pvdz", {"scf_type": "direct"}, marks=using("psi4")),
        pytest.param("qchem", "aug-cc-pvdz", {}, marks=using("qchem")),
        pytest.param("turbomole", "aug-cc-pVDZ", {}, marks=using("turbomole")),
        pytest.param("terachem_pbs", "aug-cc-pvdz", {}, marks=using("terachem_pbs")),
    ],
)
def test_sp_hf_rhf_v1v2shim(program, basis, keywords, h2o_data, schema_versions5, request, monkeypatch):
    """cfour/sp-rhf-hf/input.dat
    #! single point HF/adz on water

    """
    models, retver, _ = schema_versions5
    h2o = qcel.models.v2.Molecule.from_data(h2o_data).model_dump()
    if not from_v2(request.node.name):
        h2o["schema_version"] = 2

    # Note that this test is identical to test_sp_hf_rhf for <py314 (but with fewer checks)

    if from_v2(request.node.name):
        resi = {
            "molecule": h2o,
            "specification": {"driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords},
        }
    else:
        resi = {"molecule": h2o, "driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords}

    monkeypatch.setenv("QCNG_USE_V1V2_SHIM", "1")
    # resi = checkver_and_convert(resi, request.node.name, "pre")
    res = qcng.compute(resi, program, raise_error=True, return_dict=True, return_version=retver)
    import pprint

    pprint.pprint(res, width=200)
    # res = checkver_and_convert(res, request.node.name, "post")

    if "_v2" in request.node.name:
        assert res["input_data"]["specification"]["driver"] == "energy"
    else:
        assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    if basis == "cC-pvdZ":
        return

    # aug-cc-pvdz
    scf_tot = -76.0413815332

    atol = 1.0e-6
    assert compare_values(scf_tot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param("cfour", "cC-pvdZ", {"scf_conv": 12}, marks=using("cfour")),  # test basis handling, not results
        pytest.param("cfour", "aug-pvdz", {"scf_conv": 12}, marks=using("cfour")),
        pytest.param("cfour", "aug-pvdz", {}, marks=using("cfour")),
        pytest.param(
            "qcore",
            "aug-cc-pVDZ",
            {"coulomb_method": "direct_4idx", "exchange_method": "direct_4idx"},
            marks=using("qcore"),
        ),
        pytest.param("gamess", "accd", {"contrl__ispher": 1}, marks=using("gamess")),
        pytest.param("molpro", "aug-cc-pvdz", {}, marks=using("molpro")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True}, marks=using("nwchem")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "qc_module": "tce"}, marks=using("nwchem")),
        pytest.param("psi4", "aug-cc-pvdz", {"scf_type": "direct"}, marks=using("psi4")),
        pytest.param("qchem", "aug-cc-pvdz", {}, marks=using("qchem")),
        pytest.param("turbomole", "aug-cc-pVDZ", {}, marks=using("turbomole")),
        pytest.param("terachem_pbs", "aug-cc-pvdz", {}, marks=using("terachem_pbs")),
    ],
)
def test_sp_hf_rhf(program, basis, keywords, h2o_data, schema_versions, request):
    """cfour/sp-rhf-hf/input.dat
    #! single point HF/adz on water

    """
    models, retver, _ = schema_versions
    h2o = models.Molecule.from_data(h2o_data)

    if from_v2(request.node.name):
        resi = {
            "molecule": h2o,
            "specification": {"driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords},
        }
    else:
        resi = {"molecule": h2o, "driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords}

    resi = checkver_and_convert(resi, request.node.name, "pre")
    res = qcng.compute(resi, program, raise_error=True, return_dict=True, return_version=retver)
    res = checkver_and_convert(res, request.node.name, "post")

    if "v2" in request.node.name:
        assert res["input_data"]["specification"]["driver"] == "energy"
    else:
        assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    if basis == "cC-pvdZ":
        return

    # aug-cc-pvdz
    scf_tot = -76.0413815332

    atol = 1.0e-6
    assert compare_values(scf_tot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param(
            "cfour",
            "aug-pvdz",
            {"reference": "uhf", "occupation": [[3, 1, 1, 0], [3, 0, 1, 0]], "scf_conv": 12},
            marks=using("cfour"),
        ),
        pytest.param("cfour", "aug-pvdz", {"reference": "uhf"}, marks=using("cfour")),
        pytest.param(
            "qcore",
            "aug-cc-pVDZ",
            {"ansatz": "u", "coulomb_method": "direct_4idx", "exchange_method": "direct_4idx"},
            marks=using("qcore"),
        ),
        pytest.param("gamess", "accd", {"contrl__ispher": 1, "contrl__scftyp": "uhf"}, marks=using("gamess")),
        pytest.param("molpro", "aug-cc-pvdz", {"reference": "unrestricted"}, marks=using("molpro")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "scf__uhf": True}, marks=using("nwchem")),
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"basis__spherical": True, "qc_module": "tce", "scf__uhf": True},
            marks=using("nwchem"),
        ),
        pytest.param("psi4", "aug-cc-pvdz", {"reference": "uhf", "scf_type": "direct"}, marks=using("psi4")),
        pytest.param("qchem", "aug-cc-pvdz", {}, marks=using("qchem")),
        pytest.param("turbomole", "aug-cc-pVDZ", {}, marks=using("turbomole")),
    ],
)
def test_sp_hf_uhf(program, basis, keywords, nh2_data, schema_versions, request):
    models, retver, _ = schema_versions

    nh2 = models.Molecule.from_data(nh2_data)
    if from_v2(request.node.name):
        resi = {
            "molecule": nh2,
            "specification": {"driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords},
        }
    else:
        resi = {"molecule": nh2, "driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords}
    resi = models.AtomicInput(**resi)

    resi = checkver_and_convert(resi, request.node.name, "pre")
    res = qcng.compute(resi, program, raise_error=True, return_dict=False, return_version=retver)
    res = checkver_and_convert(res, request.node.name, "post")

    assert res.success is True
    res = res.model_dump()
    if "v2" in request.node.name:
        assert res["input_data"]["specification"]["driver"] == "energy"
    else:
        assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    scf_tot = -55.57513805253009

    atol = 1.0e-6
    assert compare_values(scf_tot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param(
            "cfour",
            "aug-pvdz",
            {"reference": "rohf", "occupation": [[3, 1, 1, 0], [3, 0, 1, 0]], "scf_conv": 12},
            marks=using("cfour"),
        ),
        pytest.param("cfour", "aug-pvdz", {"reference": "rohf"}, marks=using("cfour")),
        pytest.param("gamess", "accd", {"contrl__ispher": 1, "contrl__scftyp": "rohf"}, marks=using("gamess")),
        pytest.param("molpro", "aug-cc-pvdz", {}, marks=using("molpro")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "scf__rohf": True}, marks=using("nwchem")),
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"basis__spherical": True, "qc_module": "tce", "scf__rohf": True},
            marks=using("nwchem"),
        ),
        pytest.param("psi4", "aug-cc-pvdz", {"reference": "rohf", "scf_type": "direct"}, marks=using("psi4")),
        pytest.param("qchem", "aug-cc-pvdz", {"UNRESTRICTED": False}, marks=using("qchem")),
    ],
)
def test_sp_hf_rohf(program, basis, keywords, nh2_data, schema_versions, request):
    models, retver, models_out = schema_versions

    nh2 = models.Molecule.from_data(nh2_data)
    if from_v2(request.node.name):
        resi = {
            "molecule": nh2,
            "specification": {"driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords},
        }
    else:
        resi = {"molecule": nh2, "driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords}

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
    scf_tot = -55.570724348574

    atol = 1.0e-6
    assert compare_values(scf_tot, res["return_result"], atol=atol)
