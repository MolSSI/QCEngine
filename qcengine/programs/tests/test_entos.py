import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive
from qcengine.testing import using_entos

import qcengine as qcng

# qcenginerecords not required, skips whole file
qcer = pytest.importorskip("qcenginerecords")

# Prep globals
entos_info = qcer.get_info('entos')


@pytest.mark.parametrize('test_case', entos_info.list_test_cases())
def test_entos_output_parser(test_case):

    # Get output file data
    data = entos_info.get_test_data(test_case)
    inp = qcel.models.ResultInput.parse_raw(data["input.json"])

    output = qcng.get_program('entos', check=False).parse_output(data, inp).dict()
    output.pop("provenance", None)

    output_ref = qcel.models.Result.parse_raw(data["output.json"]).dict()
    output_ref.pop("provenance", None)

    check = compare_recursive(output_ref, output)
    assert check, check


@pytest.mark.parametrize('test_case', entos_info.list_test_cases())
def test_entos_input_formatter(test_case):

    # Get input file data
    data = entos_info.get_test_data(test_case)
    inp = qcel.models.ResultInput.parse_raw(data["input.json"])

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program('entos', check=False).build_input(inp, qcng.get_config())
    assert input_file.keys() >= {"commands", "infiles"}


@pytest.mark.parametrize('test_case', entos_info.list_test_cases())
def test_entos_input_formatter_template(test_case):

    # Get input file data
    data = entos_info.get_test_data(test_case)
    inp = qcel.models.ResultInput.parse_raw(data["input.json"])

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program('entos',  check=False).build_input(inp, qcng.get_config(), template="Test template")
    assert input_file.keys() >= {"commands", "infiles"}


@using_entos
@pytest.mark.parametrize('test_case', entos_info.list_test_cases())
def test_entos_executor(test_case):
    # Get input file data
    data = entos_info.get_test_data(test_case)
    inp = qcel.models.ResultInput.parse_raw(data["input.json"])

    # Run entos
    result = qcng.compute(inp, 'entos')
    assert result.success is True

    # Get output file data
    output_ref = qcel.models.Result.parse_raw(data["output.json"])

    atol = 1e-6
    assert compare_recursive(output_ref.return_result, result.return_result, atol=atol)
