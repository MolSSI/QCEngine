from typing import Any, Dict

from qcengine.exceptions import InputError


def muster_modelchem(method: str, driver: "DriverEnum", cc_program: str = "") -> Dict[str, Any]:
    """Converts the QC method into CFOUR keywords.

    `cc_program` is the user's CFOUR CC_PROGRAM keyword value, if any, for methods whose CALC_LEVEL spelling
    depends on the module.
    """

    method = method.lower()
    opts = {}

    if driver == "properties":
        pass
    else:
        derint = driver.derivative_int()

        if derint == 0:
            if method == "cfour":
                pass  # permit clean operation of sandwich mode
            else:
                opts["deriv_level"] = "zero"

        elif derint == 1:
            opts["deriv_level"] = "first"

        elif derint == 2:
            opts["vibration"] = "exact"

    if method == "cfour":
        pass

    elif method in ["scf", "hf"]:
        opts["calc_level"] = "scf"

    elif method == "mp2":
        opts["calc_level"] = "mp2"

    elif method == "mp3":
        opts["calc_level"] = "mp3"

    elif method == "mp4(sdq)":
        opts["calc_level"] = "sdq-mp4"

    elif method == "mp4":
        opts["calc_level"] = "mp4"

    elif method == "cisd":
        opts["calc_level"] = "cisd"

    elif method == "qcisd":
        opts["calc_level"] = "qcisd"

    elif method == "qcisd(t)":
        opts["calc_level"] = "qcisd(t)"

    elif method == "lccd":
        opts["calc_level"] = "lccd"

    elif method == "lccsd":
        opts["calc_level"] = "lccsd"

    elif method == "ccd":
        opts["calc_level"] = "ccd"

    elif method == "bccd":
        opts["calc_level"] = "b-ccd"

    elif method == "cc2":
        opts["calc_level"] = "cc2"

    elif method == "ccsd":
        opts["calc_level"] = "ccsd"

    elif method == "cc3":
        opts["calc_level"] = "cc3"

    elif method == "ccsd+t(ccsd)":
        # ECC computes [T] alongside (T) under CCSD[T] and doesn't recognize CCSD+T(CCSD)
        opts["calc_level"] = "ccsd[t]" if cc_program.lower() == "ecc" else "ccsd+t(ccsd)"

    elif method == "ccsd(t)":
        # Can't use (T) b/c bug in xsymcor lops it off
        opts["calc_level"] = "ccsd[t]"

    elif method == "a-ccsd(t)":
        opts["calc_level"] = "ccsd(t)_l"

    elif method == "bccd(t)":
        # (T) on a Brueckner reference, as validated in qcdb. B-CCD(T) runs but output not harvested for module
        opts["calc_level"] = "ccsd[t]"
        opts["brueckner"] = "on"

    elif method == "ccsdt-1a":
        opts["calc_level"] = "ccsdt-1"

    elif method == "ccsdt-1b":
        opts["calc_level"] = "ccsdt-1b"

    elif method == "ccsdt-2":
        opts["calc_level"] = "ccsdt-2"

    elif method == "ccsdt-3":
        opts["calc_level"] = "ccsdt-3"

    elif method == "ccsdt":
        opts["calc_level"] = "ccsdt"

    elif method == "ccsdt(q)":
        opts["calc_level"] = "ccsdt(q)"

    elif method == "ccsdtq":
        opts["calc_level"] = "ccsdtq"

    else:
        raise InputError(f"Method not recognized: {method}")

    return opts
