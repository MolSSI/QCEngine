from typing import Any, Dict, Tuple

from qcengine.exceptions import InputError


def format_keywords(keywords: Dict[str, Any]) -> str:
    """Form keywords deck from dictionary `keywords` where keys are CFOUR keyword ("__" separating
    any nested-module keywords) strings and values are Python formatted.

    """
    text = []
    percents = {}

    keywords = {k.upper(): v for k, v in keywords.items()}
    for key, val in sorted(keywords.items()):
        if key.startswith("%"):
            percents[key] = keywords.pop(key)
        else:
            text.append("=".join(format_keyword(key, val)))

    text = "\n".join(text)
    text = "\n\n*CFOUR(" + text + ")\n\n"

    # not yet tested with multiple % sections
    percent_text = []
    for key, val in percents.items():
        percent_text.append("\n".join(format_percent(key, val)))
    percent_text = "\n".join(percent_text)

    return text, percent_text


def format_percent(keyword: str, val: Any) -> Tuple[str, str]:
    keyword = keyword.lower()

    if isinstance(val, list):
        if type(val[0]).__name__ == "list":
            if type(val[0][0]).__name__ == "list":
                raise InputError("Option has level of array nesting inconsistent with Cfour.")
            else:
                # for %extern_pot*
                text = "\n".join([" ".join([str(inner_v) for inner_v in outer_v]) for outer_v in val])
        else:
            # option is plain 1D array
            text = " ".join([str(v) for v in val])

    # No Transform
    else:
        text = str(val).upper()

    return keyword, text


def format_keyword(keyword: str, val: Any) -> Tuple[str, str]:
    """Reformat keyword's value from Python into CFOUR-speak. Arrays are the primary target."""
    keyword = keyword.upper()

    # Transform booleans into integers
    if val is True:
        text = "1"
    elif val is False:
        text = "0"

    # Transform list from [[3, 0, 1, 1], [2, 0, 1, 0]] --> 3-0-1-1/2-0-1-0
    elif isinstance(val, list):
        if type(val[0]).__name__ == "list":
            if type(val[0][0]).__name__ == "list":
                raise InputError("Option has level of array nesting inconsistent with CFOUR.")
            else:
                # option is 2D array
                text = "/".join("-".join(map(str, no)) for no in val)
        else:
            # option is plain 1D array
            if keyword in ["ESTATE_SYM", "CFOUR_ESTATE_SYM"]:
                # [3, 1, 0, 2] --> 3/1/0/2
                text = "/".join(map(str, val))
            else:
                # [3, 1, 0, 2] --> 3-1-0-2
                text = "-".join(map(str, val))

    # Transform the basis sets that *must* be lowercase
    elif keyword in ["CFOUR_BASIS", "BASIS"] and val.upper() in [
        "SVP",
        "DZP",
        "TZP",
        "TZP2P",
        "QZ2P",
        "PZ3D2F",
        "13S9P4D3F",
    ]:
        text = str(val.lower())

    # Transform the basis sets that *must* be mixedcase

    # * manual: (keyword is cc-pVXZ or PVXZ, where X = D, T, Q, 5, and 6, respectively)
    elif keyword in ["CFOUR_BASIS", "BASIS"] and val.upper() in [
        "CC-PVDZ",
        "CC-PVTZ",
        "CC-PVQZ",
        "CC-PV5Z",
        "CC-PV6Z",
    ]:
        text = val[:4].lower() + val[4:].upper()

    # manual: (keyword is PwCVXZ, where X = D, T, and Q, respectively)
    elif keyword in ["CFOUR_BASIS", "BASIS"] and val.upper() in ["PWCVDZ", "PWCVTZ", "PWCVQZ"]:
        text = "Pw" + val[3:].upper()

    # Transform the methods that *must* be mixed case
    elif keyword in ["CFOUR_CALC_LEVEL", "CALC_LEVEL"] and val.upper() == "CCSDT-1B":
        text = "CCSDT-1b"

    # No Transform
    else:
        text = str(val).upper()

    return keyword, text
