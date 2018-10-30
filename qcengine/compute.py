"""
Integrates the computes together
"""

import copy
import time
import traceback

from . import config
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
            output_data["success"] = False
            output_data["error_message"] = "QCEngine Call Error:\nProgram {} not understood".format(program)

    output_data["stdout"] = metadata["stdout"]
    output_data["stderr"] = metadata["stderr"]
    if metadata["success"] is not True:
        output_data["success"] = False
        output_data["error_message"] = metadata["error_message"]

    # Raise an error if one exists and a user requested a raise
    if raise_error and (output_data["success"] is not True):
        raise ValueError(output_data["error_message"])

    # Fill out provenance datadata
    if "provenance" in output_data:
        output_data["provenance"].update(config.get_provenance())
    else:
        output_data["provenance"] = config.get_provenance()

    output_data["provenance"]["wall_time"] = metadata["wall_time"]

    return output_data

def compute_procedure(input_data, procedure, raise_error=False):
    """Runs a procedure (a collection of the quantum chemistry executions)

    Parameters
    ----------
    input_data : dict
        A JSON input specific to the procedure executed
    procedure : {"geometric"}
        The name of the procedure to run
    raise_error : bool, option
        Determines if compute should raise an error or not.

    Returns
    ------
    dict
        A QC Schema representation of the requested output.
    """

    input_data = copy.deepcopy(input_data)

    # Run the procedure
    comp_time = time.time()
    try:
        if procedure == "geometric":
            output_data = util.get_module_function("geometric", "run_json.geometric_run_json")(input_data)
        else:
            output_data["error_message"] = "QCEngine: Procedure {} not understood".format(procedure)
    except Exception as e:
        output_data = input_data
        output_data["success"] = False
        output_data["error_message"] = "QCEngine compute_procedure error:\n" + traceback.format_exc()

    comp_time = time.time() - comp_time

    # Raise an error if one exists and a user requested a raise
    if raise_error and ("error_message" in output_data) and (output_data["error_message"] is not False):
        raise ValueError(output_data["error_message"])

    # Fill out provenance datadata
    if "provenance" in output_data:
        output_data["provenance"].update(config.get_provenance())
    else:
        output_data["provenance"] = config.get_provenance()

    output_data["provenance"]["wall_time"] = comp_time


    return output_data
