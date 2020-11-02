"""Tests for adcc functionality"""
import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using


@pytest.fixture
def h2o():
    return qcel.models.Molecule.from_data(
        """
      O  0.0  0.000  -0.129
      H  0.0 -1.494  1.027
      H  0.0  1.494  1.027
      """
    )


# TODO: test against some ref values...
@using("adcc")
def test_run(h2o):
    inp = qcel.models.AtomicInput(
        molecule=h2o, driver="properties", model={"method": "adc2", "basis": "sto-3g"}, keywords={"n_singlets": 3}
    )
    ret = qcng.compute(inp, "adcc", raise_error=True, local_options={"ncores": 1})


# TODO: expand...
