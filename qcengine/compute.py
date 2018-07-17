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


def compute(input_data, program):
    """Executes a single quantum chemistry program given a QC Schema input.

    The full specification can be found at:
        http://molssi-qc-schema.readthedocs.io/en/latest/index.html#

    Parameters
    ----------
    input_data : dict
        A QC Schema input specification
    program : str, {"psi4", "rdkit"}
        The program to run the input under

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
        raise KeyError("Program %s not understood" % program)
    comp_time = time.time() - comp_time


    # Fill out provenance datadata
    if "provenance" in output_data:
        output_data["provenance"].update(config.get_provenance())
    else:
        output_data["provenance"] = config.get_provenance()

    output_data["provenance"]["wall_time"] = comp_time

    return output_data

def compute_procedure(input_data, procedure):
    """Runs a procedure

    Parameters
    ----------
    input_data : dict
        A JSON input specific to the procedure executed
    procedure : str, {"geometric"}
        The name of the procedure to run

    Raises
    ------
    KeyError
        Description
    """

    input_data = copy.deepcopy(input_data)

    # Run the procedure
    comp_time = time.time()
    if procedure == "geometric":
        output_data = util.get_module_function("geometric", "run_json.geometric_run_json")(input_data)
    else:
        raise KeyError("Procedure %s not understood" % procedure)
    comp_time = time.time() - comp_time

    # Fill out provenance datadata
    if "provenance" in output_data:
        output_data["provenance"].update(config.get_provenance())
    else:
        output_data["provenance"] = config.get_provenance()

    output_data["provenance"]["wall_time"] = comp_time


    return output_data