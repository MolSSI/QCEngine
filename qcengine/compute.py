"""
Integrates the computes together
"""

import copy

from .programs import get_program
from .util import compute_wrapper, handle_output_metadata, get_module_function
from .config import get_config


def compute(input_data, program, raise_error=False, capture_output=True, local_options=None):
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
    local_options : dict, optional
        A dictionary of local configuration options

    Returns
    -------
    ret : dict
        A QC Schema output

    """
    input_data = copy.deepcopy(input_data)
    config = get_config(local_options=local_options)

    # Run the program
    with compute_wrapper(capture_output=capture_output) as metadata:
        output_data = input_data
        try:
            output_data = get_program(program)(input_data, config)
        except KeyError as e:
            output_data["success"] = False
            output_data[
                "error_message"] = "QCEngine Call Error:\nProgram {} not understood.\nError Message: {}".format(
                    program, str(e))

    return handle_output_metadata(output_data, metadata, raise_error=raise_error)


def compute_procedure(input_data, procedure, raise_error=False, capture_output=True, local_options=None):
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
    local_options : dict, optional
        A dictionary of local configuration options

    Returns
    ------
    dict
        A QC Schema representation of the requested output.
    """

    input_data = copy.deepcopy(input_data)
    config = get_config(local_options=local_options)

    # Run the procedure
    with compute_wrapper(capture_output=capture_output) as metadata:
        output_data = input_data
        if procedure == "geometric":
            output_data = get_module_function("geometric", "run_json.geometric_run_json")(input_data)
        else:
            output_data["success"] = False
            output_data["error_message"] = "QCEngine Call Error:\nProcedure {} not understood".format(program)

    return handle_output_metadata(output_data, metadata, raise_error=raise_error)
