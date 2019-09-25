import qcelemental
from qcelemental.testing import compare_values
import qcengine
from qcengine import testing
import pytest


@pytest.fixture
def h2o():
    mol = qcelemental.models.Molecule.from_data("""
            O 0.000000000000     0.000000000000    -0.068516245955
            H 0.000000000000    -0.790689888800     0.543701278274
            H 0.000000000000     0.790689888800     0.543701278274
    """)
    return mol


@pytest.mark.parametrize(
    "method, keywords, ref_energy",
    [
        pytest.param('hf', {}, -75.95536954370, marks=testing.using_turbomole),
        pytest.param('pbe0', {"grid": "m5"}, -76.27371135900, marks=testing.using_turbomole),
        pytest.param('ricc2', {}, -76.1603807755, marks=testing.using_turbomole),
        pytest.param('rimp2', {}, -76.1593614075, marks=testing.using_turbomole),
    ])  # yapf: disable
def test_turbomole_energies(method, keywords, ref_energy, h2o):
    resi = {
        "molecule": h2o,
        "driver": "energy",
        "model": {
            "method": method,
            "basis": "def2-SVP",
        },
        "keywords": keywords,
    }

    res = qcengine.compute(resi, "turbomole", raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert res["success"] is True

    assert compare_values(ref_energy, res["return_result"])
