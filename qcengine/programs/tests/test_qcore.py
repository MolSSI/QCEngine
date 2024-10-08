import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive, compare_values

import qcengine as qcng
from qcengine.testing import checkver_and_convert, schema_versions, using


@using("qcore")
@pytest.mark.parametrize(
    "method, energy, gradient_norm",
    [
        ({"method": "gfn1"}, -5.765990520597323, 0.06294937),
        ({"method": "wb97xd3", "basis": "6-31g"}, -76.3624189841582, 0.04931507),
        ({"method": "hf", "basis": "6-31g"}, -75.98014477585764, 0.08866491),
    ],
)
def test_qcore_methods(method, energy, gradient_norm, schema_versions, request):
    models, retver, _ = schema_versions

    atomic_input = models.AtomicInput(
        molecule=models.Molecule(**qcng.get_molecule("water", return_dict=True)),
        model=method,
        driver="gradient",
    )

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    atomic_result = qcng.compute(atomic_input, "qcore", return_version=retver)
    atomic_result = checkver_and_convert(atomic_result, request.node.name, "post")

    assert atomic_result.success, atomic_result.error.error_message
    assert compare_values(atomic_result.properties.return_energy, energy)
    assert compare_values(np.linalg.norm(atomic_result.return_result), gradient_norm)
    assert atomic_result.wavefunction is None


@using("qcore")
def test_qcore_wavefunction(schema_versions, request):
    models, retver, _ = schema_versions

    atomic_input = models.AtomicInput(
        molecule=models.Molecule(**qcng.get_molecule("water", return_dict=True)),
        model={"method": "wb97xd3", "basis": "6-31g"},
        driver="gradient",
        protocols={"wavefunction": "all"},
    )

    atomic_input = checkver_and_convert(atomic_input, request.node.name, "pre")
    atomic_result = qcng.compute(atomic_input, "qcore", return_version=retver)
    atomic_result = checkver_and_convert(atomic_result, request.node.name, "post")

    assert atomic_result.success, atomic_result.error.error_message
    assert atomic_result.wavefunction is not None
    assert atomic_result.wavefunction.scf_orbitals_a is not None

    assert atomic_result.properties.calcinfo_nalpha == 5
    assert atomic_result.properties.calcinfo_nbasis == 13
