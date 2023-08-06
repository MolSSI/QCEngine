import numpy as np
import pytest
import qcelemental
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.programs.turbomole.harvester import parse_hessian
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


@pytest.fixture
def h2o_ricc2_def2svp():
    """NumForce calls only make sense for stationary points. So this
    geometry was optimized at the ricc2/def2-svp level of theory and
    can be used to run NumForce with ricc2."""

    mol = qcelemental.models.Molecule.from_data(
        """
        O     0.0000000    0.0000000   -0.0835835
        H     0.7501772    0.0000000    0.5210589
        H    -0.7501772    0.0000000    0.5210589
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
        pytest.param(
            "hf",
            {"scf_conv": 4, "scf_iters": 1},
            -75.95536954370,
            marks=[using("turbomole"), pytest.mark.xfail(raises=AssertionError, strict=True)],
        ),
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
        pytest.param("pbe0", {"grid": "m5"}, 0.0606266, marks=using("turbomole")),
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
    assert res.properties.return_energy

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


def assert_hessian(H, ref_eigvals, ref_size):
    w, v = np.linalg.eigh(H)
    last_eigvals = w[-3:]
    # Hessian must be symmetric
    np.testing.assert_allclose(H, H.T)
    # Check eigenvalues
    np.testing.assert_allclose(last_eigvals, ref_eigvals)
    # Hessian must be of shape (3N x 3N)
    assert H.shape == (ref_size, ref_size)


@using("turbomole")
@pytest.mark.parametrize(
    "method, keywords, ref_eigvals",
    [
        ("hf", {}, (2.00771683e-01, 7.77977644e-01, 9.91091318e-01)),
        ("pbe0", {"grid": "m5"}, (1.72092719e-01, 7.38603449e-01, 9.73783598e-01)),
        ("b-p", {"grid": "m5", "ri": True}, (1.59729409e-01, 7.21364827e-01, 9.63399519e-01)),
    ],
)
def test_turbomole_hessian(method, keywords, ref_eigvals, h2o):
    resi = {
        "molecule": h2o,
        "driver": "hessian",
        "model": {
            "method": method,
            "basis": "def2-SVP",
        },
        "keywords": keywords,
    }

    res = qcng.compute(resi, "turbomole", raise_error=True)
    H = res.return_result
    size = h2o.geometry.size

    assert res.driver == "hessian"
    assert res.success is True
    assert res.properties.return_energy
    assert_hessian(H, ref_eigvals, size)


@using("turbomole")
@pytest.mark.parametrize(
    "method, keywords, ref_eigvals",
    [
        ("ricc2", {}, (1.65405531e-01, 9.63690706e-01, 1.24676634e00)),
    ],
)
def test_turbomole_num_hessian(method, keywords, ref_eigvals, h2o_ricc2_def2svp):
    resi = {
        "molecule": h2o_ricc2_def2svp,
        "driver": "hessian",
        "model": {
            "method": method,
            "basis": "def2-SVP",
        },
        "keywords": keywords,
    }

    res = qcng.compute(resi, "turbomole", raise_error=True)
    H = res.return_result

    size = h2o_ricc2_def2svp.geometry.size

    assert res.driver == "hessian"
    assert res.success is True
    assert res.properties.return_energy
    assert_hessian(H, ref_eigvals, size)


@pytest.fixture
def h2o_nprhessian():
    return """$nprhessian
     1  1   0.6142699252  -0.0000000000   0.0000000000  -0.3071349626  -0.2479448514
     1  2  -0.0000000000  -0.3071349626   0.2479448514  -0.0000000000
     2  1  -0.0000000000   0.4365036678   0.0000000000  -0.1885017686  -0.2182518339
     2  2  -0.0000000000   0.1885017686  -0.2182518339   0.0000000000
     3  1   0.0000000000   0.0000000000  -0.0000524175  -0.0000000000  -0.0000000000
     3  2   0.0000262088  -0.0000000000   0.0000000000   0.0000262088
     4  1  -0.3071349626  -0.1885017686  -0.0000000000   0.3389423895   0.2182233100
     4  2   0.0000000000  -0.0318074269  -0.0297215414  -0.0000000000
     5  1  -0.2479448514  -0.2182518339  -0.0000000000   0.2182233100   0.2092172237
     5  2   0.0000000000   0.0297215414   0.0090346102   0.0000000000
     6  1  -0.0000000000  -0.0000000000   0.0000262088   0.0000000000   0.0000000000
     6  2  -0.0000125560  -0.0000000000   0.0000000000  -0.0000136528
     7  1  -0.3071349626   0.1885017686  -0.0000000000  -0.0318074269   0.0297215414
     7  2  -0.0000000000   0.3389423895  -0.2182233100   0.0000000000
     8  1   0.2479448514  -0.2182518339   0.0000000000  -0.0297215414   0.0090346102
     8  2   0.0000000000  -0.2182233100   0.2092172237  -0.0000000000
     9  1  -0.0000000000   0.0000000000   0.0000262088  -0.0000000000   0.0000000000
     9  2  -0.0000136528   0.0000000000  -0.0000000000  -0.0000125560
    $end"""


def test_turbomole_parse_hessian(h2o_nprhessian):
    """Test parsing of unproject Turbomole Hessian for water."""
    hessian = parse_hessian(h2o_nprhessian)
    assert hessian.shape == (9, 9)
    eigvals, _ = np.linalg.eigh(hessian)
    assert eigvals[-1] == pytest.approx(1.12157030e00)
