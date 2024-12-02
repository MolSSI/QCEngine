import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive

import qcengine as qcng
from qcengine.testing import checkver_and_convert, qcengine_records, schema_versions, using

# Prep globals
terachem_info = qcengine_records("terachem")


@pytest.mark.parametrize("test_case", terachem_info.list_test_cases())
def test_terachem_output_parser(test_case):
    # Get output file data
    data = terachem_info.get_test_data(test_case)
    inp = qcel.models.v1.AtomicInput.parse_raw(data["input.json"])

    # only qcng.compute() handles schema versions. test_data returns v1 and parse_output returns v2, so need to convert
    inp = inp.convert_v(2)
    output = qcng.get_program("terachem", check=False).parse_output(data, inp)
    output = output.convert_v(1).dict()
    output_ref = qcel.models.v1.AtomicResult.parse_raw(data["output.json"]).dict()

    # Forgiving molecule since it is now sparse
    assert compare_recursive(output_ref, output, forgive={"stdout", "provenance", "molecule", "schema_version"})


@pytest.mark.parametrize("test_case", terachem_info.list_test_cases())
def test_terachem_input_formatter(test_case):
    # Get input file data
    data = terachem_info.get_test_data(test_case)
    inp = qcel.models.v1.AtomicInput.parse_raw(data["input.json"])

    # only qcng.compute() handles schema versions. test_data returns v1 and build_input takes v2, so need to convert
    inp = inp.convert_v(2)

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program("terachem", check=False).build_input(inp, qcng.get_config())
    assert input_file.keys() >= {"commands", "infiles"}


@using("terachem")
@pytest.mark.parametrize("test_case", terachem_info.list_test_cases())
def test_terachem_executor(test_case, schema_versions, request):
    models, retver, _ = schema_versions

    # Get input file data
    data = terachem_info.get_test_data(test_case)
    inp = models.AtomicInput.parse_raw(data["input.json"])

    # Run Terachem
    inp = checkver_and_convert(inp, request.node.name, "pre")
    result = qcng.compute(inp, "terachem", return_version=retver)
    result = checkver_and_convert(result, request.node.name, "post")

    # result = qcng.get_program('terachem').compute(inp, qcng.get_config())
    assert result.success is True

    # Get output file data
    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"])

    atol = 1e-6
    if result.driver == "gradient":
        atol = 1e-3
    assert compare_recursive(output_ref.return_result, result.return_result, atol=atol)
