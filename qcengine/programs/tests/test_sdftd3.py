"""
Testing for the DFT-D3 program harness.

Most of the tests use mindless molecules for the diversity of element species
to test as much different interactions as possible.
"""

import pprint

import numpy as np
import pytest
import qcelemental as qcel

import qcengine as qcng
from qcengine.testing import checkver_and_convert, from_v2, schema_versions, uusing


@uusing("s-dftd3")
def test_dftd3_task_b97m_m01(schema_versions, request):
    import dftd3

    models, retver, _ = schema_versions

    thr = 1.0e-8

    return_result = -0.05879001214961249

    if qcel.util.parse_version(dftd3.__version__) >= qcel.util.parse_version("1.3.0"):
        atm_correction = 6.48732559e-05
        return_result -= atm_correction

    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("mindless-01", return_dict=True)),
            specification={"model": {"method": "b97m"}, "driver": "energy"},
        )
    else:
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("mindless-01", return_dict=True)),
            model={"method": "b97m"},
            driver="energy",
        )

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    atomic_result = qcng.compute(atomic_input, "s-dftd3", return_version=retver)
    atomic_result = checkver_and_convert(atomic_result, request.node.name, "post")

    print(atomic_result.return_result)
    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@uusing("s-dftd3")
@pytest.mark.parametrize("atm", ["nil", True, False])
@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"return_result": -0.024863196328457682, "level_hint": "d3bj"}),
        pytest.param({"return_result": -0.011922614341086182, "level_hint": "d3zero"}),
        pytest.param({"return_result": -0.054219891333201896, "level_hint": "d3mbj"}),
        pytest.param({"return_result": -0.087100927341320660, "level_hint": "d3mzero"}),
        pytest.param({"return_result": -0.014849063343059512, "level_hint": "d3op"}),
    ],
    ids=["d3bj", "d3zero", "d3mbj", "d3mzero", "d3op"],
)
def test_dftd3_task_pbe_m02(inp, schema_versions, request, atm):
    models, retver, _ = schema_versions

    # return to 1.0e-8 after https://github.com/MolSSI/QCEngine/issues/370
    thr = 1.0e-7

    return_result = inp["return_result"]

    # s-dftd3 v1.3 changed the 3-body default to False to match the executable interface.
    #   Toggle is keywords["params_tweaks"]["atm"] which isn't handled by the harness.
    #   Here, demonstrating pass-through to s-dftd3 1st-class QCSchema interface. Requires repeating method from model.
    atm_correction = 8.94882569e-05
    if atm == "nil":
        import dftd3

        xtra_kw = {}
        if qcel.util.parse_version(dftd3.__version__) >= qcel.util.parse_version("1.3.0"):
            return_result -= atm_correction
    else:
        xtra_kw = {"params_tweaks": {"method": "pbe", "atm": atm}}
        if atm is False:
            return_result -= atm_correction

    spec = {
        "model": {"method": "pbe"},
        "keywords": {"level_hint": inp["level_hint"], **xtra_kw},
        "driver": "energy",
    }

    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("mindless-02", return_dict=True)),
            specification=spec,
        )
    else:
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("mindless-02", return_dict=True)),
            **spec,
        )

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    atomic_result = qcng.compute(atomic_input, "s-dftd3", return_version=retver)
    atomic_result = checkver_and_convert(atomic_result, request.node.name, "post")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@uusing("s-dftd3")
def test_dftd3_task_tpss_m02(schema_versions, request):
    models, retver, _ = schema_versions

    thr = 1.0e-8

    return_result = np.array(
        [
            [-1.38459808e-04, -3.39780053e-04, +6.63854384e-06],
            [-8.49987303e-04, -1.75644402e-05, +6.92859669e-04],
            [-1.87640802e-05, +3.48459127e-04, +3.58802840e-04],
            [-1.94903751e-04, -1.27629117e-04, -6.58938185e-04],
            [+8.81075012e-05, +6.41604548e-04, +2.74678977e-04],
            [-1.64367324e-04, -2.16197141e-04, -3.84148192e-04],
            [-1.12458629e-04, -2.56329049e-04, +4.98979702e-05],
            [-1.52971321e-05, +3.02144394e-05, +1.36451189e-04],
            [+6.66230167e-04, -8.99527309e-05, +5.62299622e-04],
            [+3.09620572e-05, -1.73604390e-05, -4.31385583e-04],
            [-1.30132201e-04, -6.02965004e-04, +2.06807277e-04],
            [-6.34318034e-04, +3.32680767e-04, -5.16509587e-04],
            [-1.90415533e-04, +4.13192949e-04, -1.76537818e-04],
            [+2.65591351e-04, -1.50550857e-04, +1.39890439e-04],
            [-3.56346859e-05, -3.56851554e-04, +1.11332319e-04],
            [+1.43384741e-03, +4.09028554e-04, -3.72139482e-04],
        ]
    )

    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("mindless-02", return_dict=True)),
            specification={
                "model": {"method": ""},
                "keywords": {
                    "level_hint": "d3mbj",
                    "params_tweaks": {
                        "s8": 1.76596355,
                        "a1": 0.42822303,
                        "a2": 4.54257102,
                    },
                },
                "driver": "gradient",
            },
        )
    else:
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("mindless-02", return_dict=True)),
            model={"method": ""},
            keywords={
                "level_hint": "d3mbj",
                "params_tweaks": {
                    "s8": 1.76596355,
                    "a1": 0.42822303,
                    "a2": 4.54257102,
                },
            },
            driver="gradient",
        )

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    atomic_result = qcng.compute(atomic_input, "s-dftd3", return_version=retver)
    atomic_result = checkver_and_convert(atomic_result, request.node.name, "post")

    assert atomic_result.success
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@uusing("s-dftd3")
@pytest.mark.parametrize("atm", ["nil", True, False])
def test_dftd3_task_r2scan_m03(schema_versions, request, atm):
    models, retver, _ = schema_versions

    thr = 1.0e-8

    return_result = {
        # s9 = 0.0
        False: np.array(
            [
                [-9.35558196e-06, -7.34067064e-06, -5.91238608e-05],
                [2.53304403e-05, -5.73248858e-05, 6.46059022e-05],
                [-2.85980097e-05, 1.93908141e-05, 3.67399397e-05],
                [-7.30494324e-06, -5.48095855e-05, 3.76154049e-05],
                [1.07471133e-05, 3.03306726e-05, -2.46209798e-05],
                [-6.64928435e-05, 1.44200457e-05, -3.64523447e-05],
                [-7.88327389e-06, -3.03833956e-05, -4.75808422e-05],
                [-1.90420772e-06, -4.31459317e-05, -2.02338533e-05],
                [-3.27860042e-05, 3.68088235e-05, 2.55372601e-05],
                [4.41178314e-06, -1.11489999e-05, -4.32651442e-05],
                [7.61947667e-05, 2.64655830e-05, 2.10390818e-05],
                [-3.38080082e-05, 5.29154102e-05, -1.44995850e-04],
                [-3.85124417e-05, 4.82248838e-05, 3.75602755e-05],
                [1.01670679e-05, -3.64563841e-05, 4.13102437e-05],
                [1.02281493e-04, -2.98543913e-05, 8.22232550e-05],
                [-2.48735050e-06, 4.19080117e-05, 2.96415119e-05],
            ]
        ),
        # s9 = 1.0
        True: np.array(
            [
                [-9.15721221e-06, -8.18139252e-06, -6.04628002e-05],
                [+2.63353857e-05, -5.74347326e-05, +6.43925910e-05],
                [-2.99422439e-05, +1.88531187e-05, +3.80203831e-05],
                [-1.05590494e-05, -6.00459729e-05, +3.93481945e-05],
                [+1.59600548e-05, +3.66166973e-05, -2.69939628e-05],
                [-7.21060928e-05, +1.35991320e-05, -3.71000739e-05],
                [-9.01933781e-06, -3.41989101e-05, -4.92317946e-05],
                [-2.38625512e-06, -4.42678339e-05, -1.95513968e-05],
                [-3.09663159e-05, +3.41418638e-05, +2.51926884e-05],
                [+5.60572318e-06, -1.06356845e-05, -3.91159008e-05],
                [+7.73090102e-05, +2.75681060e-05, +2.02933984e-05],
                [-3.23109274e-05, +5.39319105e-05, -1.44309772e-04],
                [-4.05785623e-05, +5.03251549e-05, +3.92193348e-05],
                [+8.46365844e-06, -3.35282588e-05, +3.76937643e-05],
                [+1.04861677e-04, -2.99444252e-05, +8.20822571e-05],
                [-1.50951292e-06, +4.32012272e-05, +3.05230899e-05],
            ]
        ),
    }

    if atm == "nil":
        import dftd3

        xtra_kw = {}
        return_result = return_result[qcel.util.parse_version(dftd3.__version__) < qcel.util.parse_version("1.3.0")]
    else:
        xtra_kw = {"params_tweaks": {"method": "r2scan", "atm": atm}}
        return_result = return_result[atm]

    spec = {
        "keywords": {"level_hint": "d3bj", **xtra_kw},
        "driver": "gradient",
        "model": {"method": "r2scan"},
    }

    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("mindless-03", return_dict=True)),
            specification=spec,
        )
    else:
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("mindless-03", return_dict=True)),
            **spec,
        )

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    atomic_result = qcng.compute(atomic_input, "s-dftd3", return_version=retver)
    atomic_result = checkver_and_convert(atomic_result, request.node.name, "post")

    assert atomic_result.success
    assert pytest.approx(return_result, abs=thr) == atomic_result.return_result
    assert pytest.approx(atomic_result.return_result, abs=thr) == return_result


@uusing("s-dftd3")
def test_dftd3_task_unknown_method(schema_versions, request):
    models, retver, models_out = schema_versions

    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            specification={
                "keywords": {"level_hint": "d3zero"},
                "model": {"method": "non-existent-method"},
                "driver": "energy",
            },
        )
    else:
        atomic_input = models.AtomicInput(
            molecule=models.Molecule(**qcng.get_molecule("water", return_dict=True)),
            keywords={"level_hint": "d3zero"},
            model={"method": "non-existent-method"},
            driver="energy",
        )
    error = models_out.ComputeError(
        error_type="input error", error_message="No entry for 'non-existent-method' present"
    )

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    atomic_result = qcng.compute(atomic_input, "s-dftd3", return_version=retver)
    atomic_result = checkver_and_convert(
        atomic_result, request.node.name, "post", vercheck=False, cast_dict_as="FailedOperation"
    )

    print(atomic_result.error)
    assert not atomic_result.success
    assert atomic_result.error == error


@uusing("s-dftd3")
def test_dftd3_task_cold_fusion(schema_versions, request):
    models, retver, models_out = schema_versions

    if from_v2(request.node.name):
        atomic_input = models.AtomicInput(
            molecule={
                "symbols": ["Li", "Li", "Li", "Li"],
                "geometry": [
                    [-1.58746019997201, +1.58746019997201, +1.58746019997201],
                    [-1.58746019997201, +1.58746019997201, +1.58746019997201],
                    [-1.58746019997201, -1.58746019997201, -1.58746019997201],
                    [+1.58746019997201, +1.58746019997201, -1.58746019997201],
                ],
                "validated": True,  # Force a nuclear fusion input, to make dftd3 fail
            },
            specification={
                "keywords": {"level_hint": "d3zero"},
                "model": {"method": "pbe"},
                "driver": "energy",
            },
        )
    else:
        atomic_input = models.AtomicInput(
            molecule={
                "symbols": ["Li", "Li", "Li", "Li"],
                "geometry": [
                    [-1.58746019997201, +1.58746019997201, +1.58746019997201],
                    [-1.58746019997201, +1.58746019997201, +1.58746019997201],
                    [-1.58746019997201, -1.58746019997201, -1.58746019997201],
                    [+1.58746019997201, +1.58746019997201, -1.58746019997201],
                ],
                "validated": True,  # Force a nuclear fusion input, to make dftd3 fail
            },
            keywords={"level_hint": "d3zero"},
            model={"method": "pbe"},
            driver="energy",
        )
    error = models_out.ComputeError(
        error_type="input error",
        error_message="Too close interatomic distances found",
    )

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    atomic_result = qcng.compute(atomic_input, "s-dftd3", return_version=retver)
    atomic_result = checkver_and_convert(atomic_result, request.node.name, "post", vercheck=False)

    pprint.pprint(atomic_result.error.model_dump())
    assert not atomic_result.success
    assert atomic_result.error == error
