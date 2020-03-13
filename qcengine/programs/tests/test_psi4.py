import qcengine as qcng
import qcelemental as qcel
from qcengine.testing import using


@using("psi4")
def test_3c():
    mol = qcel.models.Molecule(geometry=[[0, 0, 0], [0, 1.5, 0], [0, 0, 1.5]],
                               symbols=["O", "H", "H"],
                               connectivity=[[0, 1, 1], [0, 2, 1]])
    computation = {
        "molecule": mol,
        "driver": "energy",
        "model": {"method": "B3LYP", "basis": "6-31g"}
    }
    ret = qcng.compute(computation, "psi4")
    assert ret.success is True
