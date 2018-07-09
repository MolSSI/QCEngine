"""
Integrates the computes together
"""

import copy
import time

from . import config
from .psi_compute import run_psi4
from .rdkit_compute import run_rdkit


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
        output_data = run_psi4(input_data)
    elif program == "rdkit":
        output_data = run_rdkit(input_data)
    else:
        raise KeyError("Program %s not understood" % program)
    output_data["wall_time"] = time.time() - comp_time

    # Fill out data
    if "provenance" in output_data:
        output_data["provenance"].update(config.get_provenance())
    else:
        output_data["provenance"] = config.get_provenance()

    return output_data
