import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive, compare_values

import qcengine as qcng
from qcengine.testing import checkver_and_convert, qcengine_records, schema_versions, using

qchem_info = qcengine_records("qchem")
qchem_logonly_info = qcengine_records("qchem_logonly")

qchem_forgive = [
    "root.molecule.provenance.version",
    "root.molecule.provenance.routine",
    "root.provenance.version",
    "root.provenance.routine",
    "root.provenance.creator",
    "root.extras",
    "root.molecule.extras",
    "root.schema_version",  # TODO temp until rearrangement befuddles
]


@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_output_parser(test_case):

    # Get output file data
    data = qchem_info.get_test_data(test_case)
    inp = qcel.models.v1.AtomicInput.parse_raw(data["input.json"])

    outfiles = qcel.util.deserialize(data["outfiles.msgpack"], "msgpack-ext")
    # only qcng.compute() handles schema versions. test_data returns v1 and parse_output returns v2, so need to convert
    inp = inp.convert_v(2)
    output = qcng.get_program("qchem", check=False).parse_output(outfiles, inp)
    output = output.convert_v(1).model_dump()
    output.pop("provenance", None)

    output_ref = qcel.models.v1.AtomicResult.parse_raw(data["output.json"]).model_dump()
    output_ref.pop("provenance", None)
    output_ref.pop("extras", None)
    output.pop("extras", None)

    # Modify ref to trim down total data as a molecule is now sparse
    output_ref["molecule"] = {k: v for k, v in output_ref["molecule"].items() if k in output["molecule"]}

    check, message = compare_recursive(output_ref, output, return_message=True, forgive=qchem_forgive)
    assert check, message


@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_input_formatter(test_case):

    # Get input file data
    data = qchem_info.get_test_data(test_case)
    inp = qcel.models.v1.AtomicInput.parse_raw(data["input.json"])

    # only qcng.compute() handles schema versions. test_data returns v1 and build_input takes v2, so need to convert
    inp = inp.convert_v(2)

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program("qchem", check=False).build_input(inp, qcng.get_config())
    assert input_file.keys() >= {"commands", "infiles"}


@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_input_formatter_template(test_case):

    # Get input file data
    data = qchem_info.get_test_data(test_case)
    inp = qcel.models.v1.AtomicInput.parse_raw(data["input.json"])

    # only qcng.compute() handles schema versions. test_data returns v1 and build_input takes v2, so need to convert
    inp = inp.convert_v(2)

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program("qchem", check=False).build_input(inp, qcng.get_config(), template="Test template")
    assert input_file.keys() >= {"commands", "infiles"}


@using("qchem")
@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_executor(test_case, schema_versions, request):
    models, retver, _ = schema_versions

    # Get input file data
    data = qchem_info.get_test_data(test_case)
    inp = models.AtomicInput.parse_raw(data["input.json"])

    # Run qchem
    inp = checkver_and_convert(inp, request.node.name, "pre")
    result = qcng.compute(inp, "qchem", return_version=retver)
    result = checkver_and_convert(result, request.node.name, "post")

    assert result.success is True

    # Get output file data
    output_ref = models.AtomicResult.parse_raw(data["output.json"])

    atol = 1e-6
    assert compare_recursive(output_ref.return_result, result.return_result, atol=atol)


@using("qchem")
def test_qchem_orientation(schema_versions, request):
    models, retver, _ = schema_versions

    mol = models.Molecule.from_data(
        """
        He 0.0  0.7  0.7
        He 0.0 -0.7 -0.7
        """
    )

    # Compare with rotation
    inp = {"molecule": mol, "driver": "gradient", "model": {"method": "HF", "basis": "6-31g"}}

    inp = checkver_and_convert(inp, request.node.name, "pre")
    ret = qcng.compute(inp, "qchem", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert compare_values(np.linalg.norm(ret.return_result, axis=0), [0, 0, 0.00791539])

    # Compare without rotations
    mol_noorient = mol.copy(update={"fix_orientation": True})
    inp = {"molecule": mol_noorient, "driver": "gradient", "model": {"method": "HF", "basis": "6-31g"}}

    inp = checkver_and_convert(inp, request.node.name, "pre")
    ret = qcng.compute(inp, "qchem", raise_error=True, return_version=retver)
    ret = checkver_and_convert(ret, request.node.name, "post")

    assert compare_values(np.linalg.norm(ret.return_result, axis=0), [0, 0.00559696541, 0.00559696541])


@pytest.mark.parametrize("test_case", qchem_logonly_info.list_test_cases())
def test_qchem_logfile_parser(test_case):

    # Get output file data
    data = qchem_logonly_info.get_test_data(test_case)
    outfiles = {"dispatch.out": data["qchem.out"]}
    with pytest.warns(Warning):
        output = qcng.get_program("qchem", check=False).parse_logfile(outfiles)
    # only qcng.compute() handles schema versions. above returns v2, so need to convert
    output = output.convert_v(1).model_dump()
    output["stdout"] = None

    output_ref = qcel.models.v1.AtomicResult.parse_raw(data["output.json"]).model_dump()
    for key in list(output["provenance"].keys()):
        if key not in output_ref["provenance"]:
            output["provenance"].pop(key)

    # Modify ref to trim down total data as a molecule is now sparse
    output_ref["molecule"] = {k: v for k, v in output_ref["molecule"].items() if k in output["molecule"]}

    check, message = compare_recursive(output_ref, output, return_message=True, forgive=qchem_forgive)
    assert check, message


@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_logfile_parser_qcscr(test_case):

    # Get output file data
    data = qchem_info.get_test_data(test_case)
    outfiles = qcel.util.deserialize(data["outfiles.msgpack"], "msgpack-ext")

    with pytest.warns(Warning):
        output = qcng.get_program("qchem", check=False).parse_logfile(outfiles)
    # only qcng.compute() handles schema versions. above returns v2, so need to convert
    output = output.convert_v(1).model_dump()
    output["stdout"] = None

    output_ref = qcel.models.v1.AtomicResult.parse_raw(data["output.json"]).model_dump()
    for key in list(output["provenance"].keys()):
        if key not in output_ref["provenance"]:
            output["provenance"].pop(key)

    output_ref["stdout"] = None

    # Modify ref to trim down total data as a molecule is now sparse
    output_ref["molecule"] = {k: v for k, v in output_ref["molecule"].items() if k in output["molecule"]}

    output_ref["model"]["method"] = output_ref["model"]["method"].lower()
    check, message = compare_recursive(output_ref, output, return_message=True, forgive=qchem_forgive)
    assert check, message
