import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive

import qcengine as qcng
from qcengine.testing import checkver_and_convert, qcengine_records, schema_versions, using

molpro_info = qcengine_records("molpro")


@pytest.mark.parametrize("test_case", molpro_info.list_test_cases())
def test_molpro_output_parser(test_case):

    # Get output file data
    data = molpro_info.get_test_data(test_case)
    inp = qcel.models.v1.AtomicInput.parse_raw(data["input.json"])

    # only qcng.compute() handles schema versions. test_data returns v1 and parse_output returns v2, so need to convert
    inp = inp.convert_v(2)
    output = qcng.get_program("molpro", check=False).parse_output(data, inp)
    output = output.convert_v(1).dict()
    output.pop("provenance", None)
    output.pop("schema_version", None)

    output_ref = qcel.models.v1.AtomicResult.parse_raw(data["output.json"]).dict()
    output_ref.pop("provenance", None)
    output_ref.pop("schema_version", None)

    # TODO add `skip` to compare_recursive
    check = compare_recursive(output_ref, output, forgive={"stdout"})
    assert check, check


@pytest.mark.parametrize("test_case", molpro_info.list_test_cases())
def test_molpro_input_formatter(test_case):

    # Get input file data
    data = molpro_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    # only qcng.compute() handles schema versions. test_data returns v1 and build_input takes v2, so need to convert
    inp = inp.convert_v(2)

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program("molpro", check=False).build_input(inp, qcng.get_config())
    assert input_file.keys() >= {"commands", "infiles"}


@using("molpro")
@pytest.mark.parametrize("test_case", molpro_info.list_test_cases())
def test_molpro_executor(test_case, schema_versions, request):
    models, retver, _ = schema_versions

    # Get input file data
    data = molpro_info.get_test_data(test_case)
    inp = models.AtomicInput.parse_raw(data["input.json"])

    # Run Molpro
    inp = checkver_and_convert(inp, request.node.name, "pre")
    result = qcng.compute(inp, "molpro", task_config={"ncores": 4}, return_version=retver)
    result = checkver_and_convert(result, request.node.name, "post")
    assert result.success is True

    # Get output file data
    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"])

    atol = 1e-6
    assert compare_recursive(output_ref.return_result, result.return_result, atol=atol)
