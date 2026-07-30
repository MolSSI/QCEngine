"""Extract results from MPQC stdout."""

import json
from typing import Any, Dict, Tuple

from ...exceptions import UnknownError

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
