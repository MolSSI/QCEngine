"""
Integrates the computes together
"""

import copy
import time

from . import config
from . import util

# Single computes
from . import psi_compute
from . import rdkit_compute


def compute(input_data, program, raise_error=False):
    """Executes a single quantum chemistry program given a QC Schema input.

    The full specification can be found at:
        http://molssi-qc-schema.readthedocs.io/en/latest/index.html#

    Parameters
    ----------
    input_data : dict
        A QC Schema input specification
    program : {"psi4", "rdkit"}
        The program to run the input under
    raise_error : bool, option
        Determines if compute should raise an error or not.

    Returns
    -------
    ret : dict
        A QC Schema output

    """
    input_data = copy.deepcopy(input_data)

    # Run the program
    comp_time = time.time()
    if program == "psi4":
        output_data = psi_compute.run_psi4(input_data)
    elif program == "rdkit":
        output_data = rdkit_compute.run_rdkit(input_data)
    else:
        output_data["error"] = "QCEngine: Program {} not understood".format(program)
    comp_time = time.time() - comp_time

    # Raise an error if one exists and a user requested a raise
    if raise_error and ("error" in output_data) and (output_data["error"] is not False):
        raise ValueError(output_data["error"])

    # Fill out provenance datadata
    if "provenance" in output_data:
        output_data["provenance"].update(config.get_provenance())
    else:
        output_data["provenance"] = config.get_provenance()

    output_data["provenance"]["wall_time"] = comp_time

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
    if procedure == "geometric":
        output_data = util.get_module_function("geometric", "run_json.geometric_run_json")(input_data)
    else:
        output_data["error"] = "QCEngine: Procedure {} not understood".format(procedure)
    comp_time = time.time() - comp_time

    # Raise an error if one exists and a user requested a raise
    if raise_error and ("error" in output_data) and (output_data["error"] is not False):
        raise ValueError(output_data["error"])

    # Fill out provenance datadata
    if "provenance" in output_data:
        output_data["provenance"].update(config.get_provenance())
    else:
        output_data["provenance"] = config.get_provenance()

    output_data["provenance"]["wall_time"] = comp_time


    return output_data