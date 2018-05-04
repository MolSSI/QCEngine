"""
Integrates the computes together
"""

import time

from .pass_compute import run_pass
from .psi_compute import run_psi4
from .rdkit_compute import run_rdkit
from . import dqm_config


def compute(json, program):

    # Run the program

    comp_time = time.time()
    if program == "pass":
        json = run_pass(json)
    elif program == "psi4":
        json = run_psi4(json)
    elif program == "rdkit":
        json = run_rdkit(json)
    else:
        raise KeyError("Program %s not understood" % program)
    json["wall_time"] = time.time() - comp_time

    # Fill out data
    if "provenance" in json:
        json["provenance"].update(dqm_config.get_provenance())
    else:
        json["provenance"] = dqm_config.get_provenance()

    return json
