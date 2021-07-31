from typing import Any, Dict

from ...exceptions import InputError


def muster_modelchem(method: str, derint: int) -> Dict[str, Any]:
    """Converts the QC method into GAMESS keywords."""

    method = method.lower()
    if method.endswith("-makefp"):
        method = method[:-7]
        run_makefp = True
    else:
        run_makefp = False

    opts = {}

    runtyp = {
        0: "energy",
        1: "gradient",
        2: "hessian",
        #'properties': 'prop',
    }[derint]

    if run_makefp and runtyp != "energy":
        raise InputError(f"GAMESS MAKEFP must be run with energy, not {runtyp}")

    if run_makefp:
        opts["contrl__runtyp"] = "makefp"
    else:
        opts["contrl__runtyp"] = runtyp

    if method == "gamess":
        pass

    elif method in ["scf", "hf"]:
        pass

        # opts['contrl__mplevl'] = 0
        # opts['contrl__cityp'] = 'none'
        # opts['contrl__cctyp'] = 'none'

    elif method == "mp2":
        opts["contrl__mplevl"] = 2

    elif method == "ccsd":
        opts["contrl__cctyp"] = "ccsd"

    elif method == "ccsd(t)":
        opts["contrl__cctyp"] = "ccsd(t)"

    elif method == "efp":
        opts["contrl__coord"] = "fragonly"

    else:
        raise InputError(f"Unrecognized method type '{method}'.")

    return opts
