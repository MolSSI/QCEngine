"""Map QCSchema method/driver/multiplicity to Gaussian route-section components."""

from typing import Any, Dict

from qcengine.exceptions import InputError

# rq-ef62979b rq-a88c706f rq-b49bb2a2 rq-4d763c6d
_SUPPORTED_DRIVERS = {"energy", "gradient", "hessian", "properties"}

# rq-ef62979b
_METHOD_MAP = {
    "hf": "HF",
    "scf": "HF",
    "mp2": "MP2",
    "ccsd": "CCSD",
    "ccsd(t)": "CCSD(T)",
}

# rq-ef62979b rq-c266b9bb rq-f1c116db rq-82942b90 rq-f512b60b
_DRIVER_JOB_TYPE = {
    "energy": "",
    "gradient": "Force=NoStep",
    "hessian": "Freq",
    "properties": "",
}


def muster_modelchem(method: str, driver: str, multiplicity: int) -> Dict[str, Any]:
    """Convert QCSchema method, driver, and multiplicity to Gaussian route keywords.

    Parameters
    ----------
    method:
        Method string from AtomicInput.model.method (case-insensitive).
    driver:
        One of "energy", "gradient", "hessian", "properties".
    multiplicity:
        Molecular spin multiplicity from AtomicInput.molecule.molecular_multiplicity.

    Returns
    -------
    dict with keys:
        "method_string": Gaussian method token (e.g. "HF", "UB3LYP").
        "job_type":      Gaussian job-type keyword (e.g. "Force=NoStep", "Freq").
        "extra_keywords": Dict of additional route tokens (e.g. {"Pop": "Full"}).

    Raises
    ------
    InputError
        If driver is not one of the four supported values.
    """
    # rq-b49bb2a2 rq-e008c6fe
    driver = driver.lower()
    if driver not in _SUPPORTED_DRIVERS:
        raise InputError(
            f"Driver '{driver}' is not supported by the Gaussian harness. "
            f"Supported drivers: {sorted(_SUPPORTED_DRIVERS)}"
        )

    # rq-a88c706f
    method_lower = method.lower()
    # rq-ef62979b rq-0199bb02 rq-7049aaa5 rq-546faafb rq-2052a2fa rq-92492893 rq-851d782c rq-079b1e22
    method_string = _METHOD_MAP.get(method_lower, method.upper())

    # rq-4d763c6d rq-972b6470 rq-6d73bb9e rq-c946c285 rq-25c4ad49 rq-fe13245b
    if multiplicity > 1:
        upper = method_string.upper()
        if not (upper.startswith("U") or upper.startswith("RO")):
            method_string = "U" + method_string

    # rq-c266b9bb rq-f1c116db rq-82942b90 rq-f512b60b
    job_type = _DRIVER_JOB_TYPE[driver]

    # rq-f512b60b
    extra_keywords: Dict[str, str] = {}
    if driver == "properties":
        extra_keywords["Pop"] = "Full"

    return {
        "method_string": method_string,
        "job_type": job_type,
        "extra_keywords": extra_keywords,
    }
