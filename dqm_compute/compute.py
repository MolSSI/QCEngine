"""
Integrates the computes together
"""

from .pass_compute import run_pass
from .psi_compute import run_psi4


def compute(json, program):

    if program == "pass":
        return run_pass(json)
    elif program == "psi4":
        return run_psi4(json)
    else:
        raise KeyError("Program %s not understood" % program)
