import copy

import pytest
import numpy as np
import qcelemental as qcel
from qcelemental.testing import compare_recursive

import qcengine as qcng

# qcenginerecords not required, skips whole file
qcer = pytest.importorskip("qcenginerecords")

# Prep globals
terachem_info = qcer.get_info('terachem')

@pytest.mark.parametrize('test_case', terachem_info.list_test_cases())
def test_terachem_output_parser(test_case):

     # Get output file data
     data = terachem_info.get_test_data(test_case)
     inp = qcel.models.ResultInput.parse_raw(data["input.json"])

     output = qcng.get_program('terachem').parse_output(data, inp)

     output_ref = qcel.models.Result.parse_raw(data["output.json"])

     assert compare_recursive(output_ref.dict(), output.dict()) 

@pytest.mark.parametrize('test_case', terachem_info.list_test_cases())
def test_terachem_input_formatter(test_case):

    # Get output file data
    data = terachem_info.get_test_data(test_case)
    inp = qcel.models.ResultInput.parse_raw(data["input.json"])

    # Just test that it runs for now
    input_file = qcng.get_program('terachem').build_input(inp, qcng.get_config())
    assert input_file.keys() >= {"commands", "infiles"}
