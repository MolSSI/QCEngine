"""
Integrates the computes together
"""

import copy
import time

from . import util

# Single computes
from . import psi_compute
from . import rdkit_compute


def compute(input_data, program, raise_error=False, capture_output=True):
    """Executes a single quantum chemistry program given a QC Schema input.

    The full specification can be found at:
        http://molssi-qc-schema.readthedocs.io/en/latest/index.html#

    Parameters
    ----------
    input_data : dict
        A QC Schema input specification
    program : {"psi4", "rdkit"}
        The program to run the input under
    raise_error : bool, optional
        Determines if compute should raise an error or not.
    capture_output : bool, optional
        Determines if stdout/stderr should be captured.

    Returns
    -------
    ret : dict
        A QC Schema output

    """
    input_data = copy.deepcopy(input_data)

    # Run the program
    with util.compute_wrapper(capture_output=capture_output) as metadata:
        if program == "psi4":
            output_data = psi_compute.run_psi4(input_data)
        elif program == "rdkit":
            output_data = rdkit_compute.run_rdkit(input_data)
        else:
            output_data = input_data
            output_data["success"] = False
            output_data["error_message"] = "QCEngine Call Error:\nProgram {} not understood".format(program)

    return util.handle_output_metadata(output_data, metadata, raise_error=raise_error)


def compute_procedure(input_data, procedure, raise_error=False, capture_output=True):
    """Runs a procedure (a collection of the quantum chemistry executions)

    Parameters
    ----------
    input_data : dict
        A JSON input specific to the procedure executed
    procedure : {"geometric"}
        The name of the procedure to run
    raise_error : bool, option
        Determines if compute should raise an error or not.
    capture_output : bool, optional
        Determines if stdout/stderr should be captured.

    Returns
    ------
    dict
        A QC Schema representation of the requested output.
    """

    input_data = copy.deepcopy(input_data)

    # Run the procedure
    with util.compute_wrapper(capture_output=capture_output) as metadata:
        if procedure == "geometric":
            output_data = util.get_module_function("geometric", "run_json.geometric_run_json")(input_data)
        else:
            output_data = input_data
            output_data["success"] = False
            output_data["error_message"] = "QCEngine Call Error:\nProcedure {} not understood".format(program)

    return util.handle_output_metadata(output_data, metadata, raise_error=raise_error)