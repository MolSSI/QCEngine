import pprint

import pytest

from qcelemental import constants
from qcelemental.models import Molecule
from qcelemental.testing import compare, compare_values

import qcengine as qcng
from qcengine.testing import using


@pytest.fixture
def he_tetramer():
    a2 = 2 / constants.bohr2angstroms
    return Molecule(
        symbols=["He", "He", "He", "He"],
        fragments=[[0], [1], [2], [3]],
        geometry=[0, 0, 0, 0, 0, a2, 0, a2, 0, 0, a2, a2],
    )


@using("qcmanybody")
@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param("cfour", "aug-pvdz", {"frozen_core": False}, id="cfour_conv", marks=using("cfour")),
        pytest.param(
            "gamess", "accd", {"contrl__ispher": 1, "mp2__nacore": 0}, id="gamess_conv", marks=using("gamess")
        ),
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"basis__spherical": True, "scf__thresh": 1.0e-8, "mp2__freeze": False},
            id="nwchem_conv",
            marks=using("nwchem"),
        ),
        pytest.param(
            "psi4",
            "aug-cc-pvdz",
            {"e_convergence": 1.0e-10, "d_convergence": 1.0e-10, "scf_type": "pk", "mp2_type": "conv"},
            id="psi4_conv",
            marks=using("psi4"),
        ),
    ],
)
@pytest.mark.parametrize(
    "mbe_keywords,anskey,calcinfo_nmbe",
    [
        pytest.param(
            {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 3},
            "nocp_corrected_interaction_energy_through_3_body",
            14,
            id="3b_nocp",
        ),
    ],
)
def test_nbody_he4_single(program, basis, keywords, mbe_keywords, anskey, calcinfo_nmbe, he_tetramer, request):
    from qcmanybody.models import AtomicSpecification, ManyBodyInput

    atomic_spec = AtomicSpecification(
        model={"method": "mp2", "basis": basis},
        program=program,
        driver="energy",
        keywords=keywords,
        protocols={"stdout": False},
    )
    mbe_model = ManyBodyInput(
        specification={
            "specification": atomic_spec,
            "keywords": mbe_keywords,
            "driver": "energy",
            "protocols": {"component_results": "all"},
        },
        molecule=he_tetramer,
    )

    ret = qcng.compute_procedure(mbe_model, "qcmanybody", raise_error=True)
    pprint.pprint(ret.dict(), width=200)

    assert ret.extras == {}, f"[w] extras wrongly present: {ret.extras.keys()}"

    refs = {
        "nocp_corrected_total_energy_through_1_body": -11.530668717083888,
        "nocp_corrected_total_energy_through_2_body": -11.522851206178828,
        "nocp_corrected_total_energy_through_3_body": -11.523095269671348,
        "nocp_corrected_total_energy_through_4_body": -11.523038093664368,
        "nocp_corrected_interaction_energy_through_1_body": 0.0,
        "nocp_corrected_interaction_energy_through_2_body": 0.007817510905059777,
        "nocp_corrected_interaction_energy_through_3_body": 0.0075734474125397355,
        "nocp_corrected_interaction_energy_through_4_body": 0.007630623419519367,
        "nocp_corrected_2_body_contribution_to_energy": 0.007817510905059777,
        "nocp_corrected_3_body_contribution_to_energy": -0.00024406349252004134,
        "nocp_corrected_4_body_contribution_to_energy": 5.717600697963121e-05,
    }

    ans = refs[anskey]
    ref_nmbe = calcinfo_nmbe
    atol = 1.0e-8

    for skp, ref in refs.items():
        if not "4_body" in skp:
            assert compare_values(ref, getattr(ret.properties, skp), atol=atol, label=f"[b] skprop {skp}")
        else:
            assert getattr(ret.properties, skp) is None

    assert compare_values(
        refs["nocp_corrected_total_energy_through_3_body"],
        ret.properties.nocp_corrected_total_energy,
        atol=atol,
        label=f"[d] skprop tot",
    )
    assert compare_values(
        refs["nocp_corrected_interaction_energy_through_3_body"],
        ret.properties.nocp_corrected_interaction_energy,
        atol=atol,
        label=f"[d] skprop IE",
    )

    assert compare_values(ans, ret.properties.return_energy, atol=atol, label=f"[f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=atol, label=f"[g] ret")

    assert ret.properties.calcinfo_nmbe == ref_nmbe, f"[i] {ret.properties.calcinfo_nmbe} != {ref_nmbe}"
    assert (
        len(ret.component_results) == ref_nmbe
    ), f"[k] {len(ret.component_results)} != {ref_nmbe}; mbe protocol did not take"
    if ref_nmbe > 0:
        an_atres = next(iter(ret.component_results.values()))
        assert an_atres.stdout is None, f"[l] atomic protocol did not take"


@using("qcmanybody")
@pytest.mark.parametrize(
    "qcprog",
    [
        pytest.param("cfour", marks=using("cfour")),
        pytest.param("gamess", marks=using("gamess")),
        pytest.param("nwchem", marks=using("nwchem")),
        pytest.param("psi4", marks=using("psi4")),
    ],
)
def test_bsse_ene_tu6_cp_ne2(qcprog):
    """
    from https://github.com/psi4/psi4/blob/master/tests/tu6-cp-ne2/input.dat
    Example potential energy surface scan and CP-correction for Ne2
    """
    from qcmanybody.models import ManyBodyInput

    tu6_ie_scan = {2.5: 0.757717, 3.0: 0.015685, 4.0: -0.016266}  # Ang: kcal/mol IE

    keywords = {
        "cfour": {"frozen_core": True},
        "gamess": {"contrl__ispher": 1},
        "nwchem": {"ccsd__freeze__atomic": True, "basis__spherical": True},
        "psi4": {"freeze_core": True},
    }
    basis = {
        "cfour": "aug-pvdz",
        "gamess": "accd",
        "nwchem": "aug-cc-pvdz",
        "psi4": "aug-cc-pvdz",
    }

    mbe_data = {
        "specification": {
            "specification": {
                "model": {
                    "method": "ccsd(t)",
                    "basis": basis[qcprog],
                },
                "driver": "energy",
                "program": qcprog,
                "keywords": keywords[qcprog],
            },
            "keywords": {
                "bsse_type": "cp",
            },
            "driver": "energy",
        },
        "molecule": None,
    }

    for R in tu6_ie_scan:
        nene = Molecule(
            symbols=["Ne", "Ne"], fragments=[[0], [1]], geometry=[0, 0, 0, 0, 0, R / constants.bohr2angstroms]
        )
        mbe_data["molecule"] = nene

        mbe_model = ManyBodyInput(**mbe_data)
        if qcprog == "gamess":
            with pytest.raises(RuntimeError) as exe:
                qcng.compute_procedure(mbe_model, "qcmanybody", raise_error=True)
            assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
            pytest.xfail("GAMESS can't do ghosts")

        ret = qcng.compute_procedure(mbe_model, "qcmanybody", raise_error=True)
        pprint.pprint(ret.dict(), width=200)

        assert compare_values(
            tu6_ie_scan[R], ret.return_result * constants.hartree2kcalmol, atol=1.0e-4, label=f"CP-CCSD(T) [{R:3.1f}]"
        )
        assert compare(3, ret.properties.calcinfo_nmbe, label="nmbe")
        assert compare(1, ret.properties.calcinfo_nmc, label="nmc")
        assert compare(2, ret.properties.calcinfo_nfr, label="nfr")
        assert compare(2, ret.properties.calcinfo_natom, label="nat")


@using("qcmanybody")
def test_mbe_error():
    from qcmanybody.models import ManyBodyInput

    mbe_data = {
        "specification": {
            "specification": {
                "model": {
                    "method": "nonsense",
                    "basis": "nonsense",
                },
                "driver": "energy",
                "program": "cms",
            },
            "keywords": {},
            "driver": "energy",
        },
        "molecule": None,
    }

    nene = Molecule(
        symbols=["Ne", "Ne"], fragments=[[0], [1]], geometry=[0, 0, 0, 0, 0, 3.0 / constants.bohr2angstroms]
    )
    mbe_data["molecule"] = nene

    mbe_model = ManyBodyInput(**mbe_data)

    # test 1
    with pytest.raises(RuntimeError) as exc:
        qcng.compute_procedure(mbe_model, "qcmanybody", raise_error=True)

    assert "Program cms is not registered to QCEngine" in str(exc.value)

    # test 2
    ret = qcng.compute_procedure(mbe_model, "qcmanybody")
    assert ret.success is False
    assert "Program cms is not registered to QCEngine" in ret.error.error_message


@using("psi4")
@using("qcmanybody")
@pytest.mark.parametrize(
    "optimizer,bsse_type,sio",
    [
        pytest.param("optking", "none", None, marks=using("optking")),
        # pytest.param("genoptking", "none", None, marks=using("optking_genopt")),
        # pytest.param("genoptking", "nocp", True, marks=using("optking_genopt")),
        # pytest.param("genoptking", "cp", False, marks=using("optking_genopt")),
        pytest.param("geometric", "none", None, marks=using("geometric")),
        # pytest.param("gengeometric", "none", None, marks=using("geometric_genopt")),
        # pytest.param("gengeometric", "nocp", False, marks=using("geometric_genopt")),
        # pytest.param("gengeometric", "cp", True, marks=using("geometric_genopt")),
    ],
)
def test_bsse_opt_hf_trimer(optimizer, bsse_type, sio):

    initial_molecule = Molecule.from_data(
        """
F         -0.04288        2.78905        0.00000
H          0.59079        2.03435        0.00000
--
F         -1.94320       -0.70822        0.00000
H         -1.60642        0.21789       -0.00000
--
F          2.03569       -0.60531       -0.00000
H          1.06527       -0.77673        0.00000
units ang
"""
    )

    at_spec = {
        # schema_name needed for differentiation in genopt
        "schema_name": "qcschema_input",
        "model": {
            "method": "hf",
            "basis": "6-31g",
        },
        "keywords": {
            "scf_type": "df",
        },
    }

    mbe_spec = {
        # schema_name needed for differentiation in genopt
        "schema_name": "qcschema_manybodyspecification",
        "specification": {
            "model": {
                "method": "hf",
                "basis": "6-31g",
            },
            "driver": "energy",
            "program": "psi4",
            "keywords": {},
            "protocols": {
                "stdout": False,
            },
            "extras": {
                "psiapi": True,
            },
        },
        "keywords": {
            "bsse_type": bsse_type,
            "supersystem_ie_only": sio,
        },
        "driver": "energy",
        "protocols": {
            "component_results": "all",
        },
    }

    opt_data = {
        "initial_molecule": initial_molecule,
        "input_specification": at_spec if (bsse_type == "none") else mbe_spec,
        "keywords": {
            "program": "psi4",
            "g_convergence": "nwchem_loose",
        },
        "protocols": {
            "trajectory": "initial_and_final",
        },
    }
    # from qcmanybody.models.generalized_optimization import GeneralizedOptimizationInput
    # opt_data = GeneralizedOptimizationInput(**opt_data)

    ret = qcng.compute_procedure(opt_data, optimizer, raise_error=True)

    print("FFFFFFFFFF")
    pprint.pprint(ret.dict(), width=200)

    r_fh_hb_xptd = {
        "none": 2.18 / constants.bohr2angstroms,
        "nocp": 2.18 / constants.bohr2angstroms,
        "cp": 2.27 / constants.bohr2angstroms,
    }[bsse_type]
    r_fh_computed = ret.final_molecule.measure([1, 3])
    assert (
        pytest.approx(r_fh_computed, 1.0e-2) == r_fh_hb_xptd
    ), f"hydrogen bond length computed ({r_fh_computed}) != expected ({r_fh_hb_xptd})"
    assert (
        len(ret.trajectory) == 2
    ), f"trajectory protocol did not take. len(ret.trajectory)={len(ret.trajectory)} != 2 (initial_and_final)"
    if bsse_type != "none":
        xptd_nmbe = {
            ("nocp", False): 7,
            ("nocp", True): 4,
            ("cp", False): 10,
            ("cp", True): 7,
        }[(bsse_type, sio)]
        assert xptd_nmbe == len(ret.trajectory[-1].component_results), f"mbe protocol did not take"
        assert (
            ret.trajectory[-1].component_results['["(auto)", [1, 2, 3], [1, 2, 3]]'].stdout is None
        ), f"atomic protocol did not take"


# TODO Add this back with genopt. test name in geomeTRIC is test_lif_bsse
# def test_bsse_opt_lif_dimer(optimizer, opt_keywords, bsse_type, qcprog, qc_keywords):
