"""Tests for adcc functionality"""
import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import checkver_and_convert, from_v2, schema_versions, using


@pytest.fixture
def h2o_data():
    return """
      O  0.0  0.000  -0.129
      H  0.0 -1.494  1.027
      H  0.0  1.494  1.027
      """


@using("adcc")
def test_run(h2o_data, schema_versions, request):
    models, retver, _ = schema_versions
    h2o = models.Molecule.from_data(h2o_data)

    if from_v2(request.node.name):
        inp = models.AtomicInput(
            molecule=h2o,
            specification={
                "driver": "properties",
                "model": {"method": "adc2", "basis": "sto-3g"},
                "keywords": {"n_singlets": 3},
            },
        )
    else:
        inp = models.AtomicInput(
            molecule=h2o, driver="properties", model={"method": "adc2", "basis": "sto-3g"}, keywords={"n_singlets": 3}
        )

    inp = checkver_and_convert(inp, request.node.name, "pre")
    ret = qcng.compute(
        inp, "adcc", raise_error=True, task_config={"ncores": 1}, return_dict=True, return_version=retver
    )
    ret = checkver_and_convert(ret, request.node.name, "post")
    # note dict-out

    ref_excitations = np.array([0.0693704245883876, 0.09773854881340478, 0.21481589246935925])
    ref_hf_energy = -74.45975898670224
    ref_mp2_energy = -74.67111187456267
    assert ret["success"] is True

    qcvars = ret["extras"]["qcvars"]

    assert qcvars["EXCITATION KIND"] == "SINGLET"
    assert compare_values(ref_excitations[0], ret["return_result"])
    assert compare_values(ref_hf_energy, ret["properties"]["scf_total_energy"])
    assert compare_values(ref_mp2_energy, ret["properties"]["mp2_total_energy"])
    assert compare_values(ref_excitations, qcvars["ADC2 EXCITATION ENERGIES"])
