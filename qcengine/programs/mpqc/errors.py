"""Map MPQC's C++ exception classes onto QCEngine exceptions."""

import re
from typing import Optional, Tuple

from ...exceptions import ConvergenceError, InputError, ResourceError, UnknownError

#: MPQC prints `!! MPQC exception: exception:   <ClassName>` then `description:`.
_EXCEPTION_HEADER = re.compile(r"!! MPQC exception: exception:\s+(\w+)")
_DESCRIPTION = re.compile(r"^\s*description:\s*(.*)$", re.MULTILINE)

#: Exception classes are declared in src/mpqc/util/core/exception.h.
_EXCEPTION_MAP = {
    "InputError": InputError,
    "FeatureNotImplemented": InputError,
    "FeatureDisabled": InputError,
    "MaxIterExceeded": ConvergenceError,
    "ToleranceExceeded": ConvergenceError,
    "MemAllocFailed": ResourceError,
    "LimitExceeded": ResourceError,
    "AssertionFailed": UnknownError,
    "ProgrammingError": UnknownError,
    "SystemException": UnknownError,
    "FileOperationFailed": UnknownError,
    "SyscallFailed": UnknownError,
    "Uncomputable": UnknownError,
    "AlgorithmException": UnknownError,
}


def parse_mpqc_exception(stderr: str) -> Optional[Tuple[str, str]]:
    """Pull the exception class name and description out of MPQC stderr.

    Returns ``None`` when stderr carries no MPQC exception header.
    """
    header = _EXCEPTION_HEADER.search(stderr)
    if header is None:
        return None

    description_match = _DESCRIPTION.search(stderr, header.end())
    description = description_match.group(1).strip() if description_match else ""
    return header.group(1), description


def mpqc_exception_for(stderr: str) -> Tuple[type, str]:
    """Choose the QCEngine exception class for an MPQC failure.

    Returns
    -------
    exc_class
        A QCEngine exception class. ``UnknownError`` for anything unrecognized.
    message
        MPQC's own description, with the symbolized backtrace stripped. The
        full stderr is preserved separately on the result.
    """
    parsed = parse_mpqc_exception(stderr)
    if parsed is None:
        return UnknownError, "MPQC failed without a recognizable exception header."

    class_name, description = parsed
    exc_class = _EXCEPTION_MAP.get(class_name, UnknownError)
    return exc_class, f"MPQC {class_name}: {description}"
