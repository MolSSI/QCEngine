"""Assemble a Gaussian .com input file from its components."""

from typing import Any, Dict

# Reserved user-keyword names that are handled via Link 0, not the route section.
# rq-e7897f09
_RESERVED_KEYWORDS = {"memory", "nprocs"}


# rq-3c343d2a
def build_route_line(
    method_string: str,
    basis: str,
    job_type: str,
    extra_keywords: Dict[str, str],
    user_keywords: Dict[str, Any],
) -> str:
    """Format the Gaussian route (#P) section line.

    Parameters
    ----------
    method_string:
        Gaussian method token, e.g. "HF", "UB3LYP", "CCSD(T)".
    basis:
        Basis set name, e.g. "STO-3G", "6-31G*".
    job_type:
        Job-type keyword appended after method/basis, e.g. "Force=NoStep".
        Empty string produces no extra token.
    extra_keywords:
        Additional Key=Value pairs always appended (e.g. {"Pop": "Full"}).
        A value of "" produces a bare key.
    user_keywords:
        Pass-through keywords from AtomicInput.keywords.  Any key NOT in
        _RESERVED_KEYWORDS is appended verbatim.  A falsy value produces a
        bare key.

    Returns
    -------
    str
        The full route line, e.g. "#P HF/STO-3G Force=NoStep Pop=Full".
    """
    # rq-3c343d2a rq-af032242
    parts = [f"#P {method_string}/{basis}"]

    if job_type:
        parts.append(job_type)

    # rq-cd0d6637
    for key, val in extra_keywords.items():
        parts.append(f"{key}={val}" if val else key)

    # rq-7152455e rq-e7897f09
    for key, val in user_keywords.items():
        if key.lower() in _RESERVED_KEYWORDS:
            continue
        parts.append(f"{key}={val}" if val else key)

    return " ".join(parts)


# rq-0052a2c7
def build_com_file(
    link0: Dict[str, Any],
    route_line: str,
    title: str,
    charge: int,
    multiplicity: int,
    atom_block: str,
) -> str:
    """Assemble a complete Gaussian .com input file.

    Parameters
    ----------
    link0:
        Key/value pairs written as ``%Key=Value`` lines.
        Expected keys: ``"NProcShared"``, ``"Mem"``.
    route_line:
        Full route section string, already formatted (starts with ``#P``).
    title:
        Title card (single line, must be non-blank).
    charge:
        Molecular charge.
    multiplicity:
        Spin multiplicity.
    atom_block:
        Newline-separated ``"<symbol> <x> <y> <z>"`` lines (Ångström).

    Returns
    -------
    str
        Complete .com file content ending with a trailing blank line.
    """
    # rq-0052a2c7 rq-d7e8161d
    link0_lines = "\n".join(f"%{k}={v}" for k, v in link0.items())

    # rq-342a2ca3
    chg_mult = f"{charge} {multiplicity}"

    # rq-23d8b2c6
    sections = [
        link0_lines,
        route_line,
        "",
        title,
        "",
        chg_mult,
        atom_block,
        "",
    ]
    return "\n".join(sections) + "\n"
