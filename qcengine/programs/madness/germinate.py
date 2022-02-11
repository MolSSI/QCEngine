from typing import Any, Dict, Tuple

from qcengine.exceptions import InputError

# List of XC functionals known to NWChem
_xc_functionals = [
    "hf",
    "lda",
    "acm",
    "b3lyp",
    "beckehandh",
    "pbe0",
    "becke97",
    "becke97-1",
    "becke97-2",
    "becke97-3",
    "becke97-d",
    "becke98",
    "hcth",
    "hcth120",
    "hcth147",
    "hcth407",
    "becke97gga1",
    "hcth407p",
    "mpw91",
    "mpw1k",
    "xft97",
    "cft97",
    "ft97",
    "op",
    "bop",
    "pbeop",
    "xpkzb99",
    "cpkzb99",
    "xtpss03",
    "ctpss03",
    "xctpssh",
    "b1b95",
    "bb1k",
    "mpw1b95",
    "mpwb1k",
    "pw6b95",
    "pwb6k",
    "m05",
    "m05-2x",
    "vs98",
    "m06",
    "m06-hf",
    "m06-L",
    "m06-2x",
    "HFexch",
    "becke88",
    "xperdew91",
    "xpbe96",
    "gill96",
    "lyp",
    "perdew81",
    "perdew86",
    "perdew91",
    "cpbe96",
    "pw91lda",
    "slater",
    "vwn_1",
    "vwn_2",
    "vwn_3",
    "vwn_4",
    "vwn_5",
    "vwn_1_rpa",
    "xtpss03",
    "ctpss03",
    "bc95",
    "xpw6b95",
    "xpwb6k",
    "xm05",
    "xm05-2x",
    "cpw6b95",
    "cpwb6k",
    "cm05",
    "cm05-2x",
    "xvs98",
    "cvs98",
    "xm06-L",
    "xm06-hf",
    "xm06",
    "xm06-2x",
    "cm06-L",
    "cm06-hf",
    "cm06",
    "cm06-2x",
]


# def muster_modelchem(method: str, derint: int, use_tce: bool) -> Tuple[str, Dict[str, Any]]:
def muster_modelchem(method: str, derint: int,) -> Tuple[str, Dict[str, Any]]:
    """Converts the QC method into MADNESS keywords
    Options include energy calculation with moldft
    Geometry optimiazation with moldft
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
