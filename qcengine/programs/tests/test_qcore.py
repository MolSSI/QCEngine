import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive

import qcengine as qcng
from qcengine.testing import qcengine_records, using

qcore_info = qcengine_records("qcore")


@pytest.mark.parametrize("test_case", qcore_info.list_test_cases())
def test_qcore_output_parser(test_case):

    # Get output file data
    data = qcore_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    output = qcng.get_program("qcore", check=False).parse_output(data, inp).dict()
    output.pop("provenance", None)

    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"]).dict()
    output_ref.pop("provenance", None)

    check = compare_recursive(output_ref, output)
    assert check, (output, output_ref)


@using("qcore")
@pytest.mark.parametrize("test_case", qcore_info.list_test_cases())
def test_qcore_input_formatter(test_case):

    # Get input file data
    data = qcore_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program("qcore", check=False).build_input(inp, qcng.get_config())
    assert input_file.keys() >= {"commands", "infiles"}


@pytest.mark.parametrize("test_case", qcore_info.list_test_cases())
def test_qcore_input_formatter_template(test_case):

    # Get input file data
    data = qcore_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program("qcore", check=False).build_input(inp, qcng.get_config(), template="Test template")
    assert input_file.keys() >= {"commands", "infiles"}


@using("qcore")
@pytest.mark.parametrize("test_case", qcore_info.list_test_cases())
def test_qcore_executor(test_case):
    # Get input file data
    data = qcore_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    # Run qcore
    result = qcng.compute(inp, "qcore")
    assert result.success is True

    # Get output file data
    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"])

    atol = 1e-6
    assert compare_recursive(output_ref.return_result, result.return_result, atol=atol)
