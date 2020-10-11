import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive

import qcengine as qcng
from qcengine.testing import using


@using("qcore")
@pytest.mark.parametrize(
    "method, energy, gradient_norm",
    [
        # ({"method": "gfn1"}, 5, 5),
        ({"method": "wb97xd3", "basis": "6-31g"}, 5, 5),
        # ({"method": "hf", "basis": "6-31g"}, 5, 5),
    ],
)
def test_qcore_methods(method, energy, gradient_norm):

    atomic_input = qcel.models.AtomicInput(
        molecule=qcng.get_molecule("water"),
        model=method,
        driver="gradient",
    )

    atomic_result = qcng.compute(atomic_input, "qcore")

    assert atomic_result.success, atomic_result.error.error_message
    assert pytest.approx(atomic_result.return_result, thr) == return_result
