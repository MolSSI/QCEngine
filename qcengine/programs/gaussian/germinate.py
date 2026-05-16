"""Map QCSchema method/driver/multiplicity to Gaussian route-section components."""

from typing import Any, Dict, Optional, Tuple

from qcengine.exceptions import InputError
from qcengine.programs.empirical_dispersion_resources import get_dispersion_aliases

# rq-ef62979b rq-a88c706f rq-b49bb2a2 rq-4d763c6d
_SUPPORTED_DRIVERS = {"energy", "gradient", "hessian", "properties"}

# rq-ef62979b
_METHOD_MAP = {
    "hf": "HF",
    "scf": "HF",
    "mp2": "MP2",
    # rq-2cb2c5c0 — Gaussian's analytic MP2 gradient requires the Full
    # option be attached to the method spec (e.g. "MP2(Full)/basis"); the
    # alternative "MP2=Full" route keyword works for energy only.
    "mp2(full)": "MP2(Full)",
    "ccsd": "CCSD",
    "ccsd(t)": "CCSD(T)",
    "ccsd(t,full)": "CCSD(T,Full)",
}

# rq-ef62979b rq-c266b9bb rq-f1c116db rq-82942b90 rq-f512b60b
_DRIVER_JOB_TYPE = {
    "energy": "",
    "gradient": "Force=NoStep",
    "hessian": "Freq",
    "properties": "",
}

# rq-0c1b5066 — canonical-dashlevel-key → (Gaussian EmpiricalDispersion value, formal label).
# Only the levels Gaussian natively supports via its EmpiricalDispersion keyword
# are listed here; other dispersion levels in QCEngine's alias table (D3M, D3OP,
# D4, NL, etc.) trigger an InputError.
_GAUSSIAN_DISPERSION_MAP = {
    "d2":       ("GD2",   "D2"),
    "d3zero2b": ("GD3",   "D3"),
    "d3bj2b":   ("GD3BJ", "D3(BJ)"),
}


# rq-0c1b5066 rq-f357f8db
def recognize_dispersion(method: str) -> Tuple[str, Optional[str], Optional[str], Optional[str]]:
    """Strip a recognised dispersion suffix from a method string.

    Returns ``(functional, gaussian_keyword, formal_label, canonical_key)``:

    - ``functional``     — the method with any dispersion suffix removed
      (original casing preserved). Equal to ``method`` if no suffix found.
    - ``gaussian_keyword`` — Gaussian ``EmpiricalDispersion`` value (e.g.
      ``"GD3"``), or ``None`` if no suffix found.
    - ``formal_label``   — formal dispersion label (e.g. ``"D3(BJ)"``), or
      ``None`` if no suffix found.
    - ``canonical_key``  — QCEngine canonical dashlevel key (e.g. ``"d3bj2b"``),
      or ``None`` if no suffix found.

    Raises ``InputError`` when the method ends in a recognised dispersion suffix
    that Gaussian does not natively support (e.g. D3M, D3OP, D4, NL).
    """
    method_lower = method.lower()
    aliases = get_dispersion_aliases()

    # Longest alias wins so that "d3bj" matches before "d3" + stray "bj".
    for alias in sorted(aliases, key=len, reverse=True):
        suffix = f"-{alias}"
        if method_lower.endswith(suffix):
            canonical_key = aliases[alias]
            if canonical_key not in _GAUSSIAN_DISPERSION_MAP:
                supported = sorted({fmt for _, fmt in _GAUSSIAN_DISPERSION_MAP.values()})
                raise InputError(
                    f"Dispersion level '-{alias}' (canonical: '{canonical_key}') is not "
                    f"natively supported by the Gaussian harness. "
                    f"Supported dispersion levels: {supported}."
                )
            gaussian_keyword, formal_label = _GAUSSIAN_DISPERSION_MAP[canonical_key]
            functional = method[: -len(suffix)]
            return functional, gaussian_keyword, formal_label, canonical_key

    return method, None, None, None


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

    # rq-0c1b5066 — split off any empirical-dispersion suffix BEFORE mapping
    # the functional to a Gaussian route token. recognize_dispersion raises
    # InputError on unsupported dispersion levels.
    functional, emp_disp_value, dispersion_level, _ = recognize_dispersion(method)

    # rq-a88c706f
    functional_lower = functional.lower()
    # rq-ef62979b rq-0199bb02 rq-7049aaa5 rq-546faafb rq-2052a2fa rq-92492893 rq-851d782c rq-079b1e22
    method_string = _METHOD_MAP.get(functional_lower, functional.upper())

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

    # rq-0c1b5066 — empirical dispersion: emit the Gaussian route keyword and
    # surface dispersion_level/functional_name so downstream code (harvester,
    # conflict detection) can use them.
    result: Dict[str, Any] = {
        "method_string": method_string,
        "job_type": job_type,
        "extra_keywords": extra_keywords,
    }
    if emp_disp_value is not None:
        extra_keywords["EmpiricalDispersion"] = emp_disp_value
        result["dispersion_level"] = dispersion_level
        result["functional_name"] = functional.upper()

    return result
