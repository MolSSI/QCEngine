"""Extract results from MPQC stdout."""

import json
import re
from typing import Any, Dict, Tuple

from ...exceptions import UnknownError
from ..util import PreservingDict

_OUTPUT_MARKER = "Output KeyVal (format=JSON):"

#: Paths to the property block, in preference order. The harness always emits
#: the first; the second covers deprecated-form input and mpqc_input trees.
_PROPERTY_PATHS = (("mpqc", "property"), ("property",))


def extract_output_keyval(stdout: str) -> Dict[str, Any]:
    """Parse the last ``Output KeyVal (format=JSON):`` block in ``stdout``.

    The block is indented and followed by trailer lines, so its extent is found
    by brace matching. JSON is whitespace-insensitive, so no dedent is needed.

    Raises
    ------
    UnknownError
        If no block is present, or the block does not parse.
    """
    marker = stdout.rfind(_OUTPUT_MARKER)
    if marker == -1:
        raise UnknownError(f"No '{_OUTPUT_MARKER}' block found in MPQC stdout.")

    start = stdout.find("{", marker)
    if start == -1:
        raise UnknownError(f"'{_OUTPUT_MARKER}' block found but no opening brace followed it.")

    depth = 0
    in_string = False
    escaped = False
    for index in range(start, len(stdout)):
        char = stdout[index]
        if in_string:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                in_string = False
            continue
        if char == '"':
            in_string = True
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                block = stdout[start : index + 1]
                try:
                    return json.loads(block)
                except json.JSONDecodeError as exc:
                    raise UnknownError(f"Could not parse MPQC Output KeyVal block: {exc}")

    raise UnknownError("MPQC Output KeyVal block is truncated; no matching closing brace.")


def _as_real(value: Any) -> float:
    """Coerce an MPQC-serialized scalar to a float.

    MPQC writes every scalar as a string. ExcitationEnergy is declared over
    std::complex<double>, so a two-element [real, imag] pair is accepted and
    the real part taken.
    """
    if isinstance(value, (list, tuple)):
        if abs(float(value[1])) > 1e-6:
            raise TypeError(f"Could not convert to real: {value[1]}")
        return float(value[0])
    return float(value)


def harvest_property_value(keyval: Dict[str, Any]) -> Tuple[str, Any]:
    """Read the computed property from a parsed Output KeyVal tree.

    An ``Energy`` scalar is always written as a bare string, so a ``list``
    means "array property". Were a future MPQC to serialize a *complex scalar*
    as ``["0.3", "0.0"]`` this would misread it as a two-root array; that
    cannot happen for the real-only methods in scope.

    Returns
    -------
    property_type
        ``"Energy"`` or ``"ExcitationEnergy"`` as MPQC reported it.
    value
        A ``float`` for a scalar property, a ``list`` of floats for an array.

    Raises
    ------
    UnknownError
        If no property block carries a computed value.
    """
    for path in _PROPERTY_PATHS:
        block: Any = keyval
        for key in path:
            if not isinstance(block, dict) or key not in block:
                block = None
                break
            block = block[key]
        if not isinstance(block, dict):
            continue

        property_type = block.get("type", "Energy")
        raw = block.get("value", {}).get("value")
        if raw is None:
            raise UnknownError(
                f"MPQC property block of type '{property_type}' has no computed value; "
                f"the calculation may not have completed."
            )

        if isinstance(raw, list):
            # ExcitationEnergy: an array of roots. Each element is either a
            # scalar string or a [real, imag] pair; _as_real handles both.
            return property_type, [_as_real(item) for item in raw]

        # Energy: always a bare scalar string.
        return property_type, _as_real(raw)

    raise UnknownError("MPQC Output KeyVal contains no property block at 'mpqc:property' or 'property'.")


#: Best-effort stdout patterns. Label formats drift between MPQC versions,
#: e.g., older output uses `MP2 Energy   <val>`, current uses
# `MP2 energy = <val>`, so each pattern is optional and a miss only reduces
# qcvars.
_NUMBER = r"(-?\d+\.\d+(?:[eE][-+]?\d+)?)"

_SCALAR_PATTERNS = (
    ("NUCLEAR REPULSION ENERGY", re.compile(r"Nuclear repulsion energy\s*=\s*" + _NUMBER)),
    ("MP2 CORRELATION ENERGY", re.compile(r"MP2 [Ee]nergy\s*=?\s+" + _NUMBER)),
    ("CCSD CORRELATION ENERGY", re.compile(r"CCSD [Ee]nergy\s*=?\s+" + _NUMBER)),
    ("(T) CORRECTION ENERGY", re.compile(r"\(T\) [Ee]nergy:\s*" + _NUMBER)),
)

#: SCF iteration energies; the last one is the converged total. `\s+` rather
#: than a literal indent because older MPQC tab-indents these lines and current
#: MPQC uses four spaces. It deliberately does not match `(T) Energy:`, which
#: has text between the line start and the label.
_SCF_ITERATION = re.compile(r"^\s+Energy:\s*" + _NUMBER, re.MULTILINE)


def harvest_qcvars(stdout: str, method: str) -> PreservingDict:
    """Scrape optional QCVariables from MPQC stdout.

    Never raises on a missing pattern. ``return_result`` comes from the Output
    KeyVal block, so a scraping gap degrades ``qcvars`` only.
    """
    qcvars = PreservingDict()

    for key, pattern in _SCALAR_PATTERNS:
        match = pattern.search(stdout)
        if match:
            qcvars[key] = match.group(1)

    scf_energies = _SCF_ITERATION.findall(stdout)
    if scf_energies:
        qcvars["SCF TOTAL ENERGY"] = scf_energies[-1]
        qcvars["HF TOTAL ENERGY"] = scf_energies[-1]
        qcvars["CURRENT REFERENCE ENERGY"] = scf_energies[-1]

    return qcvars


def harvest(molecule, method: str, stdout: str) -> Tuple[PreservingDict, str, Any]:
    """Read an MPQC run into QCVariables plus the primary result.

    Parameters
    ----------
    molecule
        The input molecule; used only for ``N ATOMS``. MPQC prints no
        machine-readable output geometry, so no coordinates are harvested.
    method
        Lowercase QCSchema method name, used to label method-specific qcvars.
    stdout
        Complete MPQC stdout.

    Returns
    -------
    qcvars, property_type, value
    """
    keyval = extract_output_keyval(stdout)
    property_type, value = harvest_property_value(keyval)

    qcvars = harvest_qcvars(stdout, method)
    # "N ATOMS", not "CALCINFO_NATOM": qcvar_identities_resources.py:383 maps
    # "N ATOMS" -> properties.calcinfo_natom, and every other harness in the repo
    # uses that spelling.
    qcvars["N ATOMS"] = len(molecule.symbols)

    if property_type == "Energy":
        qcvars["CURRENT ENERGY"] = value
        qcvars[f"{method.upper()} TOTAL ENERGY"] = value
    else:
        # Excitation energies are not total energies; keep them out of the
        # CURRENT ENERGY slot and expose the raw array. sCI and EOM-* both strip
        # the ground state, so root 0 is a real first excitation in either case
        # and the two need no special-casing here.
        qcvars["MPQC EXCITATION ENERGIES"] = value

    return qcvars, property_type, value
