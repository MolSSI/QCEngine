import pprint

import pytest

from qcelemental import constants
from qcelemental.models import Molecule
from qcelemental.testing import compare, compare_values

import qcengine as qcng
from qcengine.testing import using


@using("qcmanybody")
@pytest.mark.parametrize("qcprog", [
    pytest.param("cfour", marks=using("cfour")),
    pytest.param("gamess", marks=using("gamess")),
    pytest.param("nwchem", marks=using("nwchem")),
    pytest.param("psi4", marks=using("psi4")),
])
def test_tu6_cp_ne2(qcprog):
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
                    #TODO error handling cfour "basis": "aug-cc-pvdz",
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
        # TODO fix_symmetry='c1' propagate
        nene = Molecule(symbols=["Ne", "Ne"], fragments=[[0], [1]], geometry=[0, 0, 0, 0, 0, R / constants.bohr2angstroms])
        mbe_data["molecule"] = nene

        mbe_model = ManyBodyInput(**mbe_data)
        print("IIIIIII")
        pprint.pprint(mbe_model.dict(), width=200)
        if qcprog == "gamess":
            with pytest.raises(ValueError) as exe:
                qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
            assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
            pytest.xfail("GAMESS can't do ghosts")

        ret = qcng.compute_procedure(mbe_model, "qcmanybody", raise_error=True)
        # print("SSSSSSS")
        # pprint.pprint(ret.dict(), width=200)

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
                    "method": "nonsense",
                    "basis": "nonsense",
                },
                "driver": "energy",
                "program": "cms",
            },
            "keywords": {
            },
            "driver": "energy",
        },
        "molecule": None,
    }

    nene = Molecule(symbols=["Ne", "Ne"], fragments=[[0], [1]], geometry=[0, 0, 0, 0, 0, 3.0 / constants.bohr2angstroms])
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
    "optimizer,bsse_type",
    [
        pytest.param("optking", "none", marks=using("optking")),
        pytest.param("genoptking", "none", marks=using("optking")),
        pytest.param("genoptking", "nocp", marks=using("optking")),
        pytest.param("genoptking", "cp", marks=using("optking")),
    ],
)
def test_optimization_qcmanybody(optimizer, bsse_type):
    from qcmanybody.models.generalized_optimization import GeneralizedOptimizationInput

    initial_molecule = Molecule.from_data("""
F         -0.04288        2.78905        0.00000
H          0.59079        2.03435        0.00000
--
F         -1.94320       -0.70822        0.00000
H         -1.60642        0.21789       -0.00000
--
F          2.03569       -0.60531       -0.00000
H          1.06527       -0.77673        0.00000
units ang
""")

    at_spec = {
        "model": {
            "method": "hf",
            "basis": "6-31g",
        },
        "keywords": {
            "scf_type": "df",
        },
    }

    mbe_spec = {
         "specification": {
             "model": {
                 "method": "hf",
                 "basis": "6-31g",
             },
             "driver": "energy",
             "program": "psi4",
             "keywords": {},
         },
         "keywords": {
             "bsse_type": bsse_type,
         },
         "driver": "energy",
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
    opt_data = GeneralizedOptimizationInput(**opt_data)

    ret = qcng.compute_procedure(opt_data, optimizer, raise_error=True)
    import pprint
    pprint.pprint(ret.dict())

    r_fh_hb = {
        "none": 2.18 /constants.bohr2angstroms,
        "nocp": 2.18 /constants.bohr2angstroms,
        "cp": 2.27 / constants.bohr2angstroms,
    }
    r_fh_computed = ret.final_molecule.measure([1, 3])
    assert pytest.approx(r_fh_computed, 1.0e-2) == r_fh_hb[bsse_type]
