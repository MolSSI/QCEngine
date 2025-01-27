import pprint
import re

import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using

# point charge coords in psi4 and cfour must be in bohr. nwchem can be ang or bohr.
h2o_extern_bohr_qxyz = [
    [-0.834, [3.11659683, 0.0, -4.45223936]],
    [0.417, [1.02944157, 0.0, -7.18088642]],
    [0.417, [1.02944157, 0.0, -1.72359229]],
]
h2o_extern_bohr_xyzq = [[*pt[1], pt[0]] for pt in h2o_extern_bohr_qxyz]
h2o_extern_ang_xyzq = np.array(h2o_extern_bohr_xyzq)
# h2o_extern_ang_xyzq = [1.649, 0.0, -2.356, 0.545, 0.0, -3.800, 0.545, 0.0, -0.912]
h2o_extern_ang_xyzq[:, [0, 1, 2]] *= qcel.constants.bohr2angstroms
h2o_extern_ang_xyzq = h2o_extern_ang_xyzq.tolist()
h2o_extern_bohr_Nxyzq = [[len(h2o_extern_bohr_xyzq)], *h2o_extern_bohr_xyzq]

ne2_extern_bohr_qxyz = [[1.05, [0.0, 0.94486306, 0.94486306]]]
ne2_extern_bohr_xyzq = [[*pt[1], pt[0]] for pt in ne2_extern_bohr_qxyz]
ne2_extern_ang_xyzq = np.array(ne2_extern_bohr_xyzq)
ne2_extern_ang_xyzq[:, [0, 1, 2]] *= qcel.constants.bohr2angstroms
ne2_extern_ang_xyzq = ne2_extern_ang_xyzq.tolist()
ne2_extern_bohr_Nxyzq = [[len(ne2_extern_bohr_xyzq)], *ne2_extern_bohr_xyzq]


@pytest.mark.parametrize("driver", ["energy", "gradient"])
@pytest.mark.parametrize(
    "variation,program,basis,keywords",
    [
        pytest.param("h2o_plain_df", "psi4", "6-31G*", {}, marks=using("psi4")),
        pytest.param("h2o_plain_conv", "psi4", "6-31G*", {"scf_type": "direct"}, marks=using("psi4")),
        pytest.param(
            "h2o_ee_df",
            "psi4",
            "6-31G*",
            {"function_kwargs": {"external_potentials": h2o_extern_bohr_qxyz}},
            marks=using("psi4"),
        ),
        pytest.param(
            "h2o_ee_conv",
            "psi4",
            "6-31G*",
            {"scf_type": "direct", "function_kwargs": {"external_potentials": h2o_extern_bohr_qxyz}},
            marks=using("psi4"),
        ),
        pytest.param("h2o_plain_conv", "nwchem", "6-31G*", {}, marks=using("nwchem")),
        pytest.param(
            "h2o_ee_conv",
            "nwchem",
            "6-31G*",
            {
                "bq__units": "au",
                "bq__data": h2o_extern_bohr_xyzq,
                "geometry__nocenter": True,
                "geometry__autosym": False,
            },
            marks=using("nwchem"),
        ),
        pytest.param(
            "h2o_ee_conv",
            "nwchem",
            "6-31G*",
            {"bq__data": h2o_extern_ang_xyzq, "geometry__nocenter": True, "geometry__autosym": False},
            marks=using("nwchem"),
        ),
        pytest.param("h2o_plain_conv", "cfour", "6-31g*", {"spherical": 0, "scf_conv": 12}, marks=using("cfour")),
        pytest.param(
            "h2o_ee_conv",
            "cfour",
            "6-31g*",
            {"spherical": 0, "scf_conv": 12, "extern_pot": True, "%extern_pot*": h2o_extern_bohr_Nxyzq},
            marks=using("cfour"),
        ),
        # with ghost atom
        pytest.param("ne2_plain_conv", "psi4", "6-31G", {"scf_type": "direct"}, marks=using("psi4")),
        pytest.param(
            "ne2_ee_conv",
            "psi4",
            "6-31G",
            {"scf_type": "direct", "function_kwargs": {"external_potentials": ne2_extern_bohr_qxyz}},
            marks=using("psi4"),
        ),
        pytest.param("ne2_plain_conv", "nwchem", "6-31G", {}, marks=using("nwchem")),
        pytest.param(
            "ne2_ee_conv",
            "nwchem",
            "6-31G",
            {
                "bq__units": "au",
                "bq__data": ne2_extern_bohr_xyzq,
                "geometry__nocenter": True,
                "geometry__autosym": False,
            },
            marks=using("nwchem"),
        ),
        pytest.param(
            "ne2_ee_conv",
            "nwchem",
            "6-31G",
            {"bq__data": ne2_extern_ang_xyzq, "geometry__nocenter": True, "geometry__autosym": False},
            marks=using("nwchem"),
        ),
        pytest.param("ne2_plain_conv", "cfour", "6-31g", {"spherical": 0, "scf_conv": 12}, marks=using("cfour")),
        pytest.param(
            "ne2_ee_conv",
            "cfour",
            "6-31g",
            {"spherical": 0, "scf_conv": 12, "extern_pot": True, "%extern_pot*": ne2_extern_bohr_Nxyzq},
            marks=using("cfour"),
        ),
    ],
)
def test_simple_external_charges(driver, variation, program, basis, keywords, request):

    h2o_bohr = [-1.47172438, 0.0, 2.14046066, -1.25984639, 1.44393784, 3.22442268, -1.25984639, -1.44393784, 3.22442079]
    mol_ang = [-0.778803, 0.000000, 1.132683, -0.666682, 0.764099, 1.706291, -0.666682, -0.764099, 1.706290]
    ne2_bohr = [0.0, 0.0, 1.8897261254578281, 0.0, 0.0, 0.0]
    ne2_ang = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]

    water = qcel.models.Molecule(geometry=h2o_bohr, symbols=["O", "H", "H"], fix_com=True, fix_orientation=True)
    ne2 = qcel.models.Molecule(
        geometry=ne2_bohr, symbols=["Ne", "Ne"], real=[False, True], fix_com=True, fix_orientation=True
    )
    mol = water if variation.startswith("h2o_") else ne2

    atin = {"molecule": mol, "driver": driver, "model": {"method": "hf", "basis": basis}, "keywords": keywords}

    # fmt: off
    ans = {
        # plain are values w/o embedded point charges
        "h2o_plain_df": {  # copied from psi4/extern1 test
            "energy": -76.010274923509,
            "gradient": [
        -3.61281240e-03, -4.93441226e-07, -1.84830343e-02,
         1.80644240e-03,  1.27352333e-02,  9.24171070e-03,
         1.80637000e-03, -1.27347399e-02,  9.24132361e-03],
        },
        "h2o_plain_conv": {
            "energy": -76.01030124,
            "gradient": [
        -0.00360958806,  -0.000000493573, -0.018466538001,
         0.001804830243,  0.012733248901,  0.009233462596,
         0.001804757818, -0.012732755328,  0.009233075406],
        },
        "h2o_ee_df": {  # copied from psi4/extern1 test
            "energy": -76.0194112285529968,
            "gradient": [
        -6.82055778e-03, -4.91935494e-07, -1.37359563e-02,
         1.97327564e-03,  1.21536580e-02,  8.95300956e-03,
         1.97320511e-03, -1.21531644e-02,  8.95261922e-03],
        },
        "h2o_ee_conv": {
            "energy": -76.01943687,
            "gradient": [
        -0.00681662153,  -0.000000492037, -0.01372019165,
         0.001971465483,  0.012151813875,  0.008945014046,
         0.001971394934, -0.012151320165,  0.008944623631],
        },
        "ne2_plain_conv": {
            "energy": -128.4754448306251,
            "gradient": [
         0.0, 0.0,  1.61194524e-04,
         0.0, 0.0, -1.61194524e-04],
        },
        "ne2_ee_conv": {
            "energy": -128.2723042630023,
            "gradient": [
         0.0, -2.28636935e-02,  3.64658269e-02,
         0.0,  6.15278512e-01,  6.01676379e-01],
        },
    }
    # fmt: on
    thisans = {
        "energy": ans[variation]["energy"],
        "gradient": ans[variation]["gradient"],
    }

    atres = qcng.compute(atin, program, raise_error=True, return_dict=True)
    assert atres["success"] is True
    pprint.pprint(atres, width=200)

    atol = 2.0e-6

    assert compare_values(thisans["energy"], atres["properties"]["return_energy"], atol=atol, label="ene")
    assert compare_values(thisans[driver], atres["return_result"], atol=atol, label="return")
