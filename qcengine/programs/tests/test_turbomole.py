import numpy as np
import pytest
import qcelemental
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using


@pytest.fixture
def h2o():
    mol = qcelemental.models.Molecule.from_data(
        """
            O 0.000000000000     0.000000000000    -0.068516245955
            H 0.000000000000    -0.790689888800     0.543701278274
            H 0.000000000000     0.790689888800     0.543701278274
    """
    )
    return mol


@pytest.mark.parametrize(
    "method, keywords, ref_energy",
    [
        pytest.param("hf", {}, -75.95536954370, marks=using("turbomole")),
        pytest.param("pbe0", {"grid": "m5"}, -76.27371135900, marks=using("turbomole")),
        pytest.param("ricc2", {}, -76.1603807755, marks=using("turbomole")),
        pytest.param("rimp2", {}, -76.1593614075, marks=using("turbomole")),
    ],
)
def test_turbomole_energy(method, keywords, ref_energy, h2o):
    resi = {"molecule": h2o, "driver": "energy", "model": {"method": method, "basis": "def2-SVP"}, "keywords": keywords}

    res = qcng.compute(resi, "turbomole", raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert res["success"] is True

    assert compare_values(ref_energy, res["return_result"])


@pytest.mark.parametrize(
    "method, keywords, ref_norm",
    [
        pytest.param("hf", {}, 0.099340, marks=using("turbomole")),
        pytest.param("pbe0", {"grid": "m5"}, 0.060631, marks=using("turbomole")),
        pytest.param("ricc2", {}, 0.059378, marks=using("turbomole")),
        pytest.param("rimp2", {}, 0.061576, marks=using("turbomole")),
    ],
)
def test_turbomole_gradient(method, keywords, ref_norm, h2o):
    resi = {
        "molecule": h2o,
        "driver": "gradient",
        "model": {"method": method, "basis": "def2-SVP"},
        "keywords": keywords,
    }

    res = qcng.compute(resi, "turbomole", raise_error=True)

    assert res.driver == "gradient"
    assert res.success is True

    grad = res.return_result
    grad_norm = np.linalg.norm(grad)
    assert compare_values(ref_norm, grad_norm)


@using("turbomole")
def test_turbomole_ri_dsp(h2o):
    resi = {
        "molecule": h2o,
        "driver": "energy",
        "model": {"method": "b-p", "basis": "def2-SVP"},
        "keywords": {"ri": True, "d3bj": True},
    }

    res = qcng.compute(resi, "turbomole", raise_error=True)

    assert res.driver == "energy"
    assert res.success is True

    energy = res.return_result
    ref_energy = -76.36275642866
    assert compare_values(ref_energy, energy)
