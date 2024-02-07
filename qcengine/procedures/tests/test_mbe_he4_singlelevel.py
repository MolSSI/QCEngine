import pprint

import pytest

from qcelemental import constants
from qcelemental.models import Molecule
from qcelemental.models.procedures_manybody import AtomicSpecification, ManyBodyKeywords, ManyBodyInput
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using

def skprop(qcvar):
    return qcng.procedures.manybody.qcvars_to_manybodyproperties[qcvar]


he4_refs_conv = {
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":           -11.530668717083888,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":           -11.522467757090013,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":           -11.522702864080149,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":           -11.522639870651439,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":       0.008200959993875045,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":       0.007965853003739198,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":       0.008028846432448944,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":           0.008200959993875045,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":          -0.00023510699013584713,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":           6.299342870974556e-05,

        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.530668717083888,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.522851206178828,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         -11.523095269671348,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.523038093664368,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     0.007817510905059777,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":     0.0075734474125397355,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     0.007630623419519367,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         0.007817510905059777,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        -0.00024406349252004134,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":         5.717600697963121e-05,

        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.530668717083888,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.52244892169719,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         -11.52268452228489,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.522621528856181,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     0.00821979538669737,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":     0.007984194798996924,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     0.00804718822770667,
        "VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         0.00821979538669737,
        "VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        -0.00023560058770044634,
        "VMFC-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":         6.299342870974556e-05,
}

he4_refs_df = {
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":          -11.530751941948,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":          -11.522403579651,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":          -11.522640167467,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":          -11.522576639404,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":      0.008348362297,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":      0.008111774481,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":      0.008175302544,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":          0.008348362297,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":         -0.000236587816,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":          0.000063528063,

        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":        -11.530751941948,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":        -11.522760073327,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":        -11.523005411447,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":        -11.522948420000,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":    0.007991868621,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":    0.007746530501,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":    0.007803521948,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":        0.007991868621,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":       -0.000245338120,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":        0.000056991448,

        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY":        -11.530751941948,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY":        -11.522390319401,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY":        -11.522627256726,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY":        -11.522563728663,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":    0.008361622547,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":    0.008124685222,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":    0.008188213285,
        "VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":        0.008361622547,
        "VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":       -0.000236937325,
        "VMFC-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":        0.000063528063,
    }

sumdict = {
    "4b-all": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b-nocp-rtd-sio": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b-nocp-sio": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b-cp-rtd-sio": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b-cp-sio": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b-nocp-rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b-nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b-cp-rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b-cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "3b-nocp-rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b-nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b-cp-rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b-cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "2b-nocp-rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b-nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b-cp-rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b-cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "1b-nocp-rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "zero",  #"NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b-nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "zero",  #"NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b-cp-rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "zero",  #"CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b-cp": {
        "CP-CORRECTED INTERACTION ENERGY": "zero",  #"CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
# TODO table defines the general qcvar as 0 even if 1-body qcvar not available. continue?
}


@pytest.fixture
def he_tetramer():
    a2 = 2 / constants.bohr2angstroms
    return Molecule(symbols=["He", "He", "He", "He"], fragments=[[0], [1], [2], [3]], geometry=[0, 0, 0, 0, 0, a2, 0, a2, 0, 0, a2, a2])


@pytest.mark.parametrize("program,basis,keywords", [
    pytest.param("cfour", "aug-pvdz", {"frozen_core": False}, id="cfour-conv", marks=using("cfour")),
    pytest.param("gamess", "accd", {"contrl__ispher": 1, "mp2__nacore": 0}, id="gamess-conv", marks=using("gamess")),
    pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "scf__thresh": 1.0e-8, "mp2__freeze": False}, id="nwchem-conv", marks=using("nwchem")),
    pytest.param("psi4", "aug-cc-pvdz", {"e_convergence": 1.e-10, "d_convergence": 1.e-10, "scf_type": "pk", "mp2_type": "conv"}, id="psi4-conv", marks=using("psi4")),
    pytest.param("psi4", "aug-cc-pvdz", {"e_convergence": 1.e-10, "d_convergence": 1.e-10}, id="psi4-df", marks=using("psi4")),
])
@pytest.mark.parametrize("mbe_keywords,anskey,bodykeys,calcinfo_nmbe", [
    pytest.param(
        {"bsse_type": ["nocp", "cp", "vmfc"]},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv],
        65,
        id="4b-all"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "supersystem_ie_only": True},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
        5,
        id="4b-nocp-rtd-sio"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "supersystem_ie_only": True},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
        5,
        id="4b-nocp-sio"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "supersystem_ie_only": True},
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
        9,
        id="4b-cp-rtd-sio"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "supersystem_ie_only": True},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k) and "TOTAL ENERGY" not in k)],
        5,
        id="4b-cp-sio"),
## TODO add vmfc. 3b nmbe=50
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-"))],
        15,
        id="4b-nocp-rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-"))],
        15,
        id="4b-nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True},
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-"))],
        19,
        id="4b-cp-rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "TOTAL ENERGY" not in k)],
        15,
        id="4b-cp"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 3},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k))],
        14,
        id="3b-nocp-rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 3},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and "4-BODY" not in k)],
        14,
        id="3b-nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 3},
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "4-BODY" not in k)],
        18,  # bugfix: was 28
        id="3b-cp-rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 3},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "4-BODY" not in k and "TOTAL ENERGY" not in k)],
        14,
        id="3b-cp"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 2},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        10,
        id="2b-nocp-rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 2},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        10,
        id="2b-nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 2},
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        14,
        id="2b-cp-rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 2},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("4-BODY" not in k) and ("3-BODY" not in k) and "TOTAL ENERGY" not in k)],
        10,
        id="2b-cp"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 1},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("1-BODY" in k))],
        4,
        id="1b-nocp-rtd"),
# TODO fix 1b for rtd=F
#    pytest.param(
#        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 1},
#        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
#        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("1-BODY" in k))],
#        10,
#        id="1b-nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 1},
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("1-BODY" in k))],
        4,
        id="1b-cp-rtd"),
#    pytest.param(
#        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 1},
#        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
#        [k for k in he4_refs_conv if (k.startswith("CP-") and ("1-BODY" in k) and "TOTAL ENERGY" not in k)],
#        4,
#        id="1b-cp"),
])
def test_nbody_he4_single(program, basis, keywords, mbe_keywords, anskey, bodykeys, calcinfo_nmbe, he_tetramer, request):
    #! MP2/aug-cc-pvDZ many body energies of an arbitrary Helium complex,
    #   addressing 4-body formulas
    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, ...)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    if program == "gamess":
        with pytest.raises(ValueError) as exe:
            qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
        assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
        pytest.xfail("GAMESS can't do ghosts")

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print(f"SSSSSSS {request.node.name}")
    pprint.pprint(ret.model_dump(), width=200)

    _inner = request.node.name.split("[")[1].split("]")[0]
    kwdsln, progln = "-".join(_inner.split("-")[:-2]), "-".join(_inner.split("-")[-2:])
    refs = he4_refs_df if progln == "psi4-df" else he4_refs_conv
    ans = refs[anskey]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if qcv in bodykeys:
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[z] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv in sumdict["4b-all"]:
        skp = skprop(qcv)
        if qcv in sumdict[kwdsln]:
            refkey = sumdict[kwdsln][qcv]
            ref = 0.0 if refkey == "zero" else refs[refkey]
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[c] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[d] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[y] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[g] ret")
    assert ret.properties.calcinfo_nmbe == calcinfo_nmbe, f"{ret.properties.calcinfo_nmbe=} != {calcinfo_nmbe}"
