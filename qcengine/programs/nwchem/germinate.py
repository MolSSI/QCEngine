from typing import Any, Dict, Tuple, Union

from qcengine.exceptions import InputError

# List of XC functionals known to NWChem
xc_functionals = [
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
    "hfexch",
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


def muster_modelchem(method: str, derint: int, qc_module: Union[bool, str]) -> Tuple[str, Dict[str, Any]]:
    """Converts the QC method into NWChem keywords

    Args:
        method (str): Name of the QC method to use
        derint (str): Index of the run type
        qc_module (bool or str): Module routing from the `qc_module` keyword. True or "tce" uses the Tensor
            Contraction Engine. "directmp2" uses the semidirect MP2 module for MP2. False uses the default module.
    Returns:
        (str): Task command for NWChem
        (dict): Any options for NWChem
    """

    # Standardize the method name
    method = method.lower()
    opts = {}
    qc_module = qc_module.lower() if isinstance(qc_module, str) else qc_module
    use_tce = qc_module is True or qc_module == "tce"
    mdccmd = None

    # Map the run type to
    runtyp = {
        "energy": "energy",
        "gradient": "gradient",
        "hessian": "hessian",
        "properties": "property",
        "optimize": "optimize",
    }[derint]

    # Write out the theory directive
    if method == "nwchem":
        mdccmd = ""

    elif method in ["scf", "hf"]:
        mdccmd = f"task scf {runtyp}\n\n"

    elif method == "mp2":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__mp2"] = True
        elif qc_module == "directmp2":
            mdccmd = f"task direct_mp2 {runtyp}\n\n"
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

    elif method == "ccsd+t(ccsd)":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccsd(t)"] = True
        else:
            mdccmd = f"task ccsd+t(ccsd) {runtyp}\n\n"

    elif method == "ccsd(t)":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccsd(t)"] = True
        else:
            mdccmd = f"task ccsd(t) {runtyp}\n\n"

    elif method == "ccsdt":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccsdt"] = True
        else:
            mdccmd = f"task ccsdt {runtyp}\n\n"

    elif method in ["lccd", "lccsd", "cisd", "qcisd", "cc2", "ccsdtq"]:
        # TCE-only methods
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts[f"tce__{method}"] = True

    elif method == "tddft":
        mdccmd = f"task tddft {runtyp}\n\n"

    elif method in ["sodft", "direct_mp2", "rimp2", "mcscf", "selci", "md", "pspw", "band"]:
        raise InputError(f'Method "{method}" not yet supported by QCEngine')

    elif method == "tce":
        raise InputError(
            f"Do not specify TCE as a method. Instead specify the desired method " f'as a keyword and "qc_module=True".'
        )

    elif method.split()[0] in xc_functionals:
        opts["dft__xc"] = method
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__"] = "dft"
        else:
            mdccmd = f"task dft {runtyp}\n\n"

    elif method == "pbe":
        opts["dft__xc"] = "xpbe96 cpbe96"
        mdccmd = f"task dft {runtyp}\n\n"

    elif method == "b2plyp":
        opts["dft__xc"] = "hfexch 0.53 becke88 0.47 lyp 0.73 mp2 0.27"
        opts["dft__dftmp2"] = True
        mdccmd = f"task dft {runtyp}\n\n"

    elif method == "b3lyp5":
        opts["dft__xc"] = "hfexch 0.2 slater 0.8 becke88 nonlocal 0.72 vwn_5 0.190 lyp 0.81"
        mdccmd = f"task dft {runtyp}\n\n"

    elif method == "dft":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__"] = "dft"
        else:
            mdccmd = f"task dft {runtyp}\n\n"

    else:
        raise InputError(f"Method not recognized: {method}")

    if mdccmd is None:
        raise InputError(f"Method {method} not available with qc_module={qc_module}. Try qc_module=tce.")

    return mdccmd, opts
