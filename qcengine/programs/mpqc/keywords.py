"""Translate flat QCSchema keywords into MPQC's nested KeyVal tree."""

import copy
from typing import Any, Dict, Tuple

#: Keys consumed before ``__`` expansion because they do not obey the
#: generic nesting rule. See the design spec, "keywords.py".
MPQC_INPUT_KEY = "mpqc_input"
MPQC_ENV_KEY = "mpqc_env"
PROPERTY_PREFIX = "property__"


def deep_merge(base: Dict[str, Any], overlay: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge ``overlay`` onto ``base``; ``overlay`` wins.

    Neither argument is mutated.
    """
    merged = copy.deepcopy(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def extract_reserved(
    keywords: Dict[str, Any],
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Split reserved keys out of ``keywords``.

    Returns
    -------
    remaining : dict
        Every non-reserved keyword, still in flat ``block__key`` form.
    mpqc_input : dict
        A complete KeyVal tree replacing the generated skeleton, or ``{}``.
    mpqc_env : dict
        Environment-variable overrides, or ``{}``. Never reaches the input file.
    property_opts : dict
        ``property__*`` keywords with the prefix stripped, destined for
        ``mpqc:property``, or ``{}``.
    """
    remaining: Dict[str, Any] = {}
    mpqc_input: Dict[str, Any] = {}
    mpqc_env: Dict[str, Any] = {}
    property_opts: Dict[str, Any] = {}

    for key, value in keywords.items():
        if key == MPQC_INPUT_KEY:
            mpqc_input = copy.deepcopy(value)
        elif key == MPQC_ENV_KEY:
            mpqc_env = copy.deepcopy(value)
        elif key.startswith(PROPERTY_PREFIX):
            property_opts[key[len(PROPERTY_PREFIX) :]] = value
        else:
            remaining[key] = value

    return remaining, mpqc_input, mpqc_env, property_opts


def format_keywords(opts: Dict[str, Any]) -> Dict[str, Any]:
    """Expand flat ``block__key`` keywords into a nested dict.

    ``{"wfn__eom__manifold": "2h1p"}`` becomes
    ``{"wfn": {"eom": {"manifold": "2h1p"}}}``. Keys without ``__`` pass
    through at the top level.
    """
    tree: Dict[str, Any] = {}
    for key, value in opts.items():
        parts = key.split("__")
        cursor = tree
        for part in parts[:-1]:
            cursor = cursor.setdefault(part, {})
        cursor[parts[-1]] = value
    return tree
