"""
Integrates the computes together
"""

import copy
import time

from . import config
from .psi_compute import run_psi4
from .rdkit_compute import run_rdkit


def compute(input_data, program):

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
