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
) -> Tuple[str, Dict[str, Any]]:
    """Converts the QC method into MADNESS keywords
    Options include energy calculation with moldft
    Geometry optimization with moldft
    propreties calculation with molresponse...runs moldft then molresponse

    Args:
       method (str): Name of the QC method to use
       derint (str): Index of the run type
    Returns:
       (str): Task command for MADNESS
       (dict): Any options for MADNESS
    """

    # Standardize the method name
    method = method.lower()

    opts = {}

    # Map the run type to
    # runtyp = {"energy": "energy", "gradient": "gradient", "hessian": "hessian", "properties": "property"}[derint]
    runtyp = {"energy": "energy", "optimization": "gopt", "hessian": "hessian", "properties": "molresponse"}[derint]

    # Write out the theory directive
    if runtyp == "energy":
        if method == "optimization":
            opts["dft__gopt"] = True
        elif method.split()[0] in _xc_functionals:
            opts["dft__xc"] = method
        else:
            raise InputError(f"Method not recognized: {method}")
        mdccmd = f""
    elif runtyp == "molresponse":
        if method.split()[0] in _xc_functionals:
            opts["dft__xc"] = method
            opts["response__xc"] = method
            opts["response__archive"] = "restartdata"
        else:
            raise InputError(f"Method not recognized: {method}")
        mdccmd = f"response"  ## we will split the options with the word response later

    ## all we have to do is add options to the dft block in order to change the run type
    ## default in energy
    # do nothing
    return mdccmd, opts
