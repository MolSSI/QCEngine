import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using


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


@using("madness")
@pytest.mark.parametrize(
    "program,basis,keywords",
    [pytest.param("madness", None, {"dft__k": 7, "dft__aobasis": "sto-3g", "dft__econv": 1.0000e-05}),],
)
@using("madness")
def test_mad_hf(program, basis, keywords, h2o):
    resi = {"molecule": h2o, "driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)
    print(res["stdout"])

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # k=7
    scf_tot = -76.06718632

    atol = 1.0e-5
    assert compare_values(scf_tot, res["return_result"], atol=atol)
