from typing import Any, Dict, Tuple


def muster_modelchem(method: str, derint: int, use_tce: bool) -> Tuple[str, Dict[str, Any]]:
    """Converts the QC method into NWChem keywords."""

    method = method.lower()
    opts = {}

    runtyp = {
        0: "energy",
        1: "gradient",
        2: "hessian",
        #'properties': 'prop',
    }[derint]

    if method == "nwchem":
        pass

    elif method in ["scf", "hf"]:
        mdccmd = f"task scf {runtyp}\n\n"

    elif method == "mp2":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__mp2"] = True
        else:
            mdccmd = f"task mp2 {runtyp}\n\n"

    elif method == "mp3":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__mp3"] = True

    elif method == "mp4":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__mp4"] = True

    elif method == "ccd":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccd"] = True

    elif method == "ccsd":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccsd"] = True
        else:
            mdccmd = f"task ccsd {runtyp}\n\n"

    elif method == "ccsdt":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccsdt"] = True
        else:
            mdccmd = f"task ccsdt {runtyp}\n\n"

    elif method == "ccsd(t)":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccsd(t)"] = True
        else:
            mdccmd = f"task ccsd(t) {runtyp}\n\n"

    return mdccmd, opts
