"""Q-Chem extension callbacks for the shared cclib harness."""

import math
import os
import re
import stat
import tempfile
from typing import TYPE_CHECKING, Any, Dict, Mapping, Type

from qcelemental.util import safe_version

from ...config import TaskConfig
from ...exceptions import InputError, ResourceError, UnknownError
from ...util import execute
from .base import CCLibHarness, Job, ProgramDefinition, _validate_input_subset

if TYPE_CHECKING:
    from qcelemental.models.v2 import AtomicInput


QCHEM_RESERVED = {
    "JOBTYPE",
    "METHOD",
    "BASIS",
    "MEM_TOTAL",
    "INPUT_BOHR",
    "SCF_FINAL_PRINT",
    "PRINT_GENERAL_BASIS",
    "PRINT_ORBITALS",
    "MOLDEN_FORMAT",
}


def _render_qchem_scalar(key: str, value: Any) -> str:
    """Render one validated native Q-Chem keyword value."""

    if "\n" in key or "\r" in key:
        raise InputError(f"Q-Chem keyword contains a newline: {key!r}")
    if isinstance(value, str):
        if "\n" in value or "\r" in value:
            raise InputError(f"Q-Chem keyword {key!r} contains a newline")
        return value
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float) and math.isfinite(value):
        return str(value)
    raise InputError(f"Q-Chem keyword {key!r} must have a string, bool, int, or finite float value")


def build_input(input_model: "AtomicInput", config: TaskConfig, executable: str) -> Job:
    """Generate a Q-Chem job from the supported QCSchema subset."""

    driver, method, basis = _validate_input_subset(input_model)
    jobtype = {"energy": "sp", "gradient": "force", "hessian": "freq"}[driver.lower()]

    user_options: Dict[str, str] = {}
    for key, value in input_model.specification.keywords.items():
        if (
            not isinstance(key, str)
            or not key
            or any(character.isspace() or not character.isprintable() for character in key)
        ):
            raise InputError(f"Q-Chem keyword name must be exactly one non-empty native token: {key!r}")
        normalized_key = key.upper()
        if normalized_key in QCHEM_RESERVED:
            raise InputError(f"Q-Chem keyword {key!r} is reserved by CCLibHarness")
        if normalized_key in user_options:
            raise InputError(f"Q-Chem keyword collision after case normalization: {key!r}")
        user_options[normalized_key] = _render_qchem_scalar(key, value)

    molecule = input_model.molecule
    charge = int(molecule.molecular_charge)
    multiplicity = int(molecule.molecular_multiplicity)
    geometry_lines = [
        f"{symbol} {str(coordinates[0])} {str(coordinates[1])} {str(coordinates[2])}"
        for symbol, coordinates in zip(molecule.symbols, molecule.geometry)
    ]
    rem_lines = [
        f"JOBTYPE {jobtype}",
        f"METHOD {method}",
        f"BASIS {basis}",
        f"MEM_TOTAL {int(config.memory * 1024)}",
        "INPUT_BOHR TRUE",
        "SCF_FINAL_PRINT 2",
        "PRINT_GENERAL_BASIS TRUE",
        "PRINT_ORBITALS TRUE",
        "MOLDEN_FORMAT FALSE",
    ]
    rem_lines.extend(f"{key} {user_options[key]}" for key in sorted(user_options))
    input_text = (
        "$comment\n"
        "QCEngine CCLibHarness\n"
        "$end\n\n"
        "$molecule\n"
        f"{charge} {multiplicity}\n"
        + "\n".join(geometry_lines)
        + "\n$end\n\n"
        "$rem\n"
        + "\n".join(rem_lines)
        + "\n$end\n"
    )
    return Job(
        command=[executable, "-nt", str(config.ncores), "dispatch.in", "dispatch.out"],
        infiles={"dispatch.in": input_text},
        outfiles=["dispatch.out"],
        input_filename="dispatch.in",
        output_filename="dispatch.out",
        input_text=input_text,
        executable=executable,
    )


def select_output(outputs: Mapping[str, Any]) -> str:
    """Select Q-Chem's primary output file from execution results."""

    try:
        output = outputs["dispatch.out"]
    except KeyError as exc:
        raise UnknownError("Q-Chem did not produce dispatch.out") from exc
    if output is None:
        raise UnknownError("Q-Chem did not produce dispatch.out")
    return str(output)


def _path_has_mode(path: str, mode: int) -> bool:
    """Return whether a path has and permits the requested owner mode."""

    try:
        return bool(os.stat(path).st_mode & mode) and os.access(path, os.R_OK if mode == stat.S_IRUSR else os.X_OK)
    except OSError:
        return False


def preflight(executable: str, environment: Mapping[str, str]) -> Dict[str, str]:
    """Validate all Q-Chem resources and return an isolated child environment."""

    child_environment = dict(environment)
    invalid = []

    for variable in ("QC", "QCAUX"):
        value = child_environment.get(variable)
        if not value or not os.path.isdir(value) or not _path_has_mode(value, stat.S_IRUSR):
            invalid.append(f"{variable} is not set or does not identify a readable directory")

    qcprog = child_environment.get("QCPROG")
    if (
        not qcprog
        or not os.path.isfile(qcprog)
        or not _path_has_mode(qcprog, stat.S_IRUSR)
        or not _path_has_mode(qcprog, stat.S_IXUSR)
    ):
        invalid.append("QCPROG is not set or does not identify a readable executable program driver")

    if not os.path.isfile(executable) or not _path_has_mode(executable, stat.S_IXUSR):
        invalid.append(f"resolved qchem executable is not runnable: {executable}")

    if invalid:
        details = "; ".join(f"Q-Chem environment variable {item}" for item in invalid)
        raise ResourceError(f"{details}. Initialize the Q-Chem environment before running cclib-qchem.")

    child_environment.setdefault("QCSCRATCH", tempfile.gettempdir())
    return child_environment


def apply_scratch_environment(environment: Dict[str, str], scratch_directory: str) -> None:
    """Bind Q-Chem's native scratch variable to managed scratch."""

    environment["QCSCRATCH"] = scratch_directory


def probe(executable: str, environment: Mapping[str, str]) -> str:
    """Verify Q-Chem identity and return its normalized version."""

    success, outputs = execute(
        [executable, "version.in"],
        {"version.in": "$rem\n$end\n"},
        environment=dict(environment),
        timeout=15,
    )
    output = f"{outputs.get('stdout') or ''}\n{outputs.get('stderr') or ''}"
    if not success:
        raise ResourceError(f"Q-Chem identity/version probe failed for {executable}")
    if "A Quantum Leap Into The Future Of Chemistry" not in output or "Q-Chem" not in output:
        raise ResourceError(f"Executable {executable} failed Q-Chem identity verification")

    match = re.search(r"Q-Chem(?:\s+version:)?\s+([0-9]+(?:\.[0-9A-Za-z]+)+)", output, re.IGNORECASE)
    if match is None:
        raise ResourceError(f"Could not parse the Q-Chem version from executable {executable}")
    return safe_version(match.group(1))


def parser_type() -> Type[Any]:
    """Load the expected Q-Chem parser class lazily."""

    from cclib.parser.qchemparser import QChem

    return QChem


QCHEM_DEFINITION = ProgramDefinition(
    selector="cclib-qchem",
    executable="qchem",
    minimum_version="5.1",
    input_filename="dispatch.in",
    output_filename="dispatch.out",
    parser_name="QChem",
    parser_type=parser_type,
    normal_termination="Thank you very much for using Q-Chem",
    managed_scratch_suffix="_cclib_qchem_scratch",
    scratch_environment=apply_scratch_environment,
    generator=build_input,
    probe=probe,
    output_selector=select_output,
    preflight=preflight,
)


class QChemCCLibHarness(CCLibHarness):
    """Run Q-Chem and convert its output exclusively through cclib."""

    definition = QCHEM_DEFINITION
    _defaults = {**CCLibHarness._defaults, "name": "cclib-qchem", "node_parallel": False}
