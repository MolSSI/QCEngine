import pprint

import pytest

from qcelemental import constants
from qcelemental.models import Molecule
from qcelemental.models.procedures_layered import ManyBodyInput
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using


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
        if qcprog == "gamess":
            with pytest.raises(ValueError) as exe:
                qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
            assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
            pytest.xfail("GAMESS can't do ghosts")

        ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
        # print("SSSSSSS")
        # pprint.pprint(ret.model_dump(), width=200)

        assert compare_values(
            tu6_ie_scan[R], ret.return_result * constants.hartree2kcalmol, atol=1.0e-4, label=f"CP-CCSD(T) [{R:3.1f}]"
        )
