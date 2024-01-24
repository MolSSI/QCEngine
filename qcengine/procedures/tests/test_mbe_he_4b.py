import pprint

import pytest

from qcelemental import constants
from qcelemental.models import Molecule
from qcelemental.models.procedures_layered import AtomicSpecification, ManyBodyKeywords, ManyBodyInput
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


@pytest.fixture
def he_tetramer():
    a2 = 2 / constants.bohr2angstroms
    return Molecule(symbols=["He", "He", "He", "He"], fragments=[[0], [1], [2], [3]], geometry=[0, 0, 0, 0, 0, a2, 0, a2, 0, 0, a2, a2])


_qcprog_lanes = [
    pytest.param("cfour", "aug-pvdz", {"frozen_core": False}, id="cfour-conv", marks=using("cfour")),
    pytest.param("gamess", "accd", {"contrl__ispher": 1, "mp2__nacore": 0}, id="gamess-conv", marks=using("gamess")),
    pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "scf__thresh": 1.0e-8, "mp2__freeze": False}, id="nwchem-conv", marks=using("nwchem")),
    pytest.param("psi4", "aug-cc-pvdz", {"e_convergence": 1.e-10, "d_convergence": 1.e-10, "scf_type": "pk", "mp2_type": "conv"}, id="psi4-conv", marks=using("psi4")),
    pytest.param("psi4", "aug-cc-pvdz", {"e_convergence": 1.e-10, "d_convergence": 1.e-10}, id="psi4-df", marks=using("psi4")),
    ]


@pytest.mark.parametrize("program,basis,keywords", _qcprog_lanes)
def test_nbody_he_4b(program, basis, keywords, he_tetramer, request):
    #! MP2/aug-cc-pvDZ many body energies of an arbitrary Helium complex,
    #   addressing 4-body formulas

    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, bsse_type=["nocp",
    #                 "cp", "vmfc"], return_wfn=True)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_keywords = ManyBodyKeywords(bsse_type=["nocp", "cp", "vmfc"])
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    if program == "gamess":
        with pytest.raises(ValueError) as exe:
            qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
        assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
        pytest.xfail("GAMESS can't do ghosts")

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print("SSSSSSS [1] All,4b")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_df if "psi4-df" in request.node.name else he4_refs_conv
    ans = refs["NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY"]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[1a] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[1b] skprop {skp}")

    for qcv, ref in {
        "NOCP-CORRECTED TOTAL ENERGY": refs["NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY"],
        "NOCP-CORRECTED INTERACTION ENERGY": ans,
        "CP-CORRECTED TOTAL ENERGY": refs["CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY"],
        "CP-CORRECTED INTERACTION ENERGY": refs["CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY"],
        "VMFC-CORRECTED TOTAL ENERGY": refs["VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY"],
        "VMFC-CORRECTED INTERACTION ENERGY": refs["VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY"],
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[1c] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[1d] skprop {skp}")

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[1e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[1f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[1g] ret")
    assert ret.properties.calcinfo_nmbe == 65, f"{ret.properties.calcinfo_nmbe=} != 65"


@pytest.mark.parametrize("program,basis,keywords", _qcprog_lanes)
def test_nbody_he_3b_nocp(program, basis, keywords, he_tetramer, request):
    #! MP2/aug-cc-pvDZ many body energies of an arbitrary Helium complex,
    #   addressing 4-body formulas, 3-body NOCP subset

    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, bsse_type="nocp",
    #                  return_wfn=True, max_nbody=3)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_keywords = ManyBodyKeywords(bsse_type=["nocp"], max_nbody=3)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)
    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print("SSSSSSS [2] NOCP,3b")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_df if "psi4-df" in request.node.name else he4_refs_conv
    ans = refs["NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if "4-BODY" not in qcv and qcv.startswith("NOCP-"):
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[2a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[2b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[2] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "NOCP-CORRECTED TOTAL ENERGY": refs["NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY"],
        "NOCP-CORRECTED INTERACTION ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[2c] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[2d] skprop {skp}")

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[2e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[2f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[1g] ret")
    assert ret.properties.calcinfo_nmbe == 14, f"{ret.properties.calcinfo_nmbe=} != 14"


@pytest.mark.parametrize("program,basis,keywords", _qcprog_lanes)
def test_nbody_he_3b_cp_rtd(program, basis, keywords, he_tetramer, request):
    #! MP2/aug-cc-pvDZ many body energies of an arbitrary Helium complex,
    #   addressing 4-body formulas, 3-body CP subset, total data

    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, bsse_type="cp",
    #                  return_wfn=True, max_nbody=3, return_total_data=True)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_keywords = ManyBodyKeywords(bsse_type=["cp"], max_nbody=3, return_total_data=True)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    if program == "gamess":
        with pytest.raises(ValueError) as exe:
            qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
        assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
        pytest.xfail("GAMESS can't do ghosts")

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print("SSSSSSS [3] CP,3b,rtd=T")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_df if "psi4-df" in request.node.name else he4_refs_conv
    ans = refs["CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY"]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if "4-BODY" not in qcv and qcv.startswith("CP-"):
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[3a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[3b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[3] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CP-CORRECTED TOTAL ENERGY": ans,
        "CP-CORRECTED INTERACTION ENERGY": refs["CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"],
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[3c] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[3d] skprop {skp}")

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[3e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[3f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[1g] ret")
    assert ret.properties.calcinfo_nmbe == 18, f"{ret.properties.calcinfo_nmbe=} != 18"  # bugfix: was 28


@pytest.mark.parametrize("program,basis,keywords", _qcprog_lanes)
def test_nbody_he_3b_cp(program, basis, keywords, he_tetramer, request):
    #! MP2/aug-cc-pvDZ many body energies of an arbitrary Helium complex,
    #   addressing 4-body formulas, 3-body CP subset

    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, bsse_type="cp",
    #                  return_wfn=True, max_nbody=3, return_total_data=False)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_keywords = ManyBodyKeywords(bsse_type=["cp"], max_nbody=3, return_total_data=False)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    if program == "gamess":
        with pytest.raises(ValueError) as exe:
            qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
        assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
        pytest.xfail("GAMESS can't do ghosts")

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print("SSSSSSS [4] CP,3b,rtd=F")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_df if "psi4-df" in request.node.name else he4_refs_conv
    ans = refs["CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if "4-BODY" not in qcv and qcv.startswith("CP-") and "TOTAL ENERGY" not in qcv:
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[4a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[4b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[4] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CP-CORRECTED INTERACTION ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[4c] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[4d] skprop {skp}")

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[4e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[4f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[1g] ret")
    assert ret.properties.calcinfo_nmbe == 14, f"{ret.properties.calcinfo_nmbe=} != 14"


@pytest.mark.parametrize("program,basis,keywords", _qcprog_lanes)
def test_nbody_he_3b_vmfc(program, basis, keywords, he_tetramer, request):
    #! MP2/aug-cc-pvDZ many body energies of an arbitrary Helium complex,
    #   addressing 4-body formulas, 3-body VMFC subset

    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, bsse_type="vmfc",
    #                  return_wfn=True, max_nbody=3)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_keywords = ManyBodyKeywords(bsse_type=["vmfc"], max_nbody=3)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    if program == "gamess":
        with pytest.raises(ValueError) as exe:
            qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
        assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
        pytest.xfail("GAMESS can't do ghosts")

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print("SSSSSSS [5] VMFC,3b")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_df if "psi4-df" in request.node.name else he4_refs_conv
    ans = refs["VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY"]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if "4-BODY" not in qcv and qcv.startswith("VMFC-"):
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[5a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[5b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[5] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "VMFC-CORRECTED TOTAL ENERGY": refs["VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY"],
        "VMFC-CORRECTED INTERACTION ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[5c] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[5d] skprop {skp}")

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[5e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[5f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[1g] ret")
    assert ret.properties.calcinfo_nmbe == 50, f"{ret.properties.calcinfo_nmbe=} != 50"


@pytest.mark.parametrize("program,basis,keywords", _qcprog_lanes)
def test_nbody_he_4b_nocp_scm(program, basis, keywords, he_tetramer, request):
    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, bsse_type="nocp",
    #                  return_total_data=False, short_circuit_mbe=True, return_wfn=True)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_keywords = ManyBodyKeywords(bsse_type=["nocp"], return_total_data=False, short_circuit_mbe=True)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print("SSSSSSS [9] NOCP,4b,rtd=F,scm=T")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_df if "psi4-df" in request.node.name else he4_refs_conv
    ans = refs["NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY"]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if ("THROUGH 4-BODY" in qcv or "THROUGH 1-BODY" in qcv) and qcv.startswith("NOCP-"):
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[9a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[9b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[9] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "NOCP-CORRECTED INTERACTION ENERGY": ans,
        "NOCP-CORRECTED TOTAL ENERGY": refs["NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY"],
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[9c] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[9d] skprop {skp}")

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[9e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[9f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[1g] ret")
    assert ret.properties.calcinfo_nmbe == 5, f"{ret.properties.calcinfo_nmbe=} != 5"


@pytest.mark.parametrize("program,basis,keywords", _qcprog_lanes)
def test_nbody_he_4b_nocp_rtd_scm(program, basis, keywords, he_tetramer, request):
    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, bsse_type="nocp",
    #                  return_total_data=True, short_circuit_mbe=True, return_wfn=True)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_keywords = ManyBodyKeywords(bsse_type="nocp", return_total_data=True, short_circuit_mbe=True)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print("SSSSSSS [10] NOCP,4b,rtd=T,scm=T")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_df if "psi4-df" in request.node.name else he4_refs_conv
    ans = refs["NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY"]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if ("THROUGH 4-BODY" in qcv or "THROUGH 1-BODY" in qcv) and qcv.startswith("NOCP-"):
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[10a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[10b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[10] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "NOCP-CORRECTED TOTAL ENERGY": ans,
        "NOCP-CORRECTED INTERACTION ENERGY": refs["NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY"],
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[10c] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[10d] skprop {skp}")

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[10e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[10f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[1g] ret")
    assert ret.properties.calcinfo_nmbe == 5, f"{ret.properties.calcinfo_nmbe=} != 5"


@pytest.mark.parametrize("program,basis,keywords", _qcprog_lanes)
def test_nbody_he_4b_cp_scm(program, basis, keywords, he_tetramer, request):
    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, bsse_type="cp",
    #                  return_total_data=False, short_circuit_mbe=True, return_wfn=True)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_keywords = ManyBodyKeywords(bsse_type=["cp"], return_total_data=False, short_circuit_mbe=True)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    if program == "gamess":
        with pytest.raises(ValueError) as exe:
            qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
        assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
        pytest.xfail("GAMESS can't do ghosts")

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print("SSSSSSS [11] CP,4b,rtd=F,scm=T")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_df if "psi4-df" in request.node.name else he4_refs_conv
    ans = refs["CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY"]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if "INTERACTION ENERGY THROUGH 4-BODY" in qcv and qcv.startswith("CP-"):
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[11a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[11b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[11] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CP-CORRECTED INTERACTION ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[11c] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[11d] skprop {skp}")

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[11e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[11f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[11g] ret")
    assert ret.properties.calcinfo_nmbe == 5, f"{ret.properties.calcinfo_nmbe=} != 5"


@pytest.mark.parametrize("program,basis,keywords", _qcprog_lanes)
def test_nbody_he_4b_cp_rtd_scm(program, basis, keywords, he_tetramer, request):
    # e, wfn = energy('MP2/aug-cc-pVDZ', molecule=he_tetramer, bsse_type="cp",
    #                  return_total_data=True, short_circuit_mbe=True, return_wfn=True)

    atomic_spec = AtomicSpecification(model={"method": "mp2", "basis": basis}, program=program, driver="energy", keywords=keywords)
    mbe_keywords = ManyBodyKeywords(bsse_type="cp", return_total_data=True, short_circuit_mbe=True)
    mbe_model = ManyBodyInput(specification={"specification": atomic_spec, "keywords": mbe_keywords, "driver": "energy"}, molecule=he_tetramer)

    if program == "gamess":
        with pytest.raises(ValueError) as exe:
            qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
        assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
        pytest.xfail("GAMESS can't do ghosts")

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print("SSSSSSS [12] CP,4b,rtd=T,scm=T")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_df if "psi4-df" in request.node.name else he4_refs_conv
    ans = refs["CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY"]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if ("THROUGH 4-BODY" in qcv or "THROUGH 1-BODY" in qcv) and qcv.startswith("CP-"):
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[12a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[12b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[12] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CP-CORRECTED TOTAL ENERGY": ans,
        "CP-CORRECTED INTERACTION ENERGY": refs["CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY"],
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[12c] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[12d] skprop {skp}")

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[12e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[12f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[12g] ret")
    assert ret.properties.calcinfo_nmbe == 9, f"{ret.properties.calcinfo_nmbe=} != 9"
