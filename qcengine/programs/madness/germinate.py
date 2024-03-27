from typing import Any, Dict, Tuple

from qcengine.exceptions import InputError

# List of XC functionals known to NWChem
_xc_functionals = [
    "hf",
    "lda",
    "b3lyp",
]


def muster_modelchem(
    method: str,
    derint: int,
) -> Dict[str, Any]:
    """Converts the QC method into MADNESS keywords
    Options include energy calculation with moldft
    Geometry optimization with moldft
    propreties calculation with molresponse...runs moldft then molresponse

    Args:
       method (str): Name of the QC method to use
       derint (str): Index of the run type
    Returns:
       (dict): Any options for MADNESS
    """

    # Standardize the method name
    method = method.lower()

    # Initialize the options
    opts = {}
    moldft_opts = opts["moldft"] = {}
    molresponse_opts = opts["molresponse"] = {}

    runtyp = {"energy": "energy", "optimization": "gopt", "hessian": "hessian", "properties": "response"}[derint]

    # Write out the theory directive
    if runtyp == "energy":
        if method == "optimization":
            moldft_opts["dft__gopt"] = True

        elif method.split()[0] in _xc_functionals:
            moldft_opts["dft__xc"] = method

        else:
            raise InputError(f"Method not recognized: {method}")
    elif runtyp == "response":
        if method.split()[0] in _xc_functionals:

            moldft_opts["dft__xc"] = method
            molresponse_opts["response__xc"] = method
            molresponse_opts["response__archive"] = "restartdata"

        else:
            raise InputError(f"Method not recognized: {method}")

    return opts
