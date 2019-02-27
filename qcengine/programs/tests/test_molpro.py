
import copy

import pytest
import numpy as np
import qcelemental as qcel
from qcelemental.testing import compare_recursive

from qcengine.programs import get_program

# qcenginerecords not required, skips whole file
qcer = pytest.importorskip("qcenginerecords")

# Prep globals
molpro_info = qcer.get_info('molpro')

@pytest.mark.parametrize('test_case', molpro_info.list_test_cases())
def test_molpro_output_parser(test_case):

    # Get output file data
    data = molpro_info.get_test_data(test_case)
    inp = qcel.models.ResultInput.parse_raw(data["input.json"])

    output = get_program('molpro').parse_output(data, inp)

    output_ref = qcel.models.Result.parse_raw(data["output.json"])

    assert compare_recursive(output_ref.dict(), output.dict())