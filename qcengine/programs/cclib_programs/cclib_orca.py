"""ORCA extension callbacks for the shared cclib harness."""

import re
from typing import TYPE_CHECKING, Any, Dict, Mapping, Optional, Type

from qcelemental import constants
from qcelemental.util import safe_version

from ...config import TaskConfig
from ...exceptions import InputError, ResourceError, UnknownError
from ...util import execute
from .base import CCLibHarness, Job, ProgramDefinition, _validate_input_subset

if TYPE_CHECKING:
    from qcelemental.models.v2 import AtomicInput


_ORCA_OUTPUT_DEFAULTS = [
    "PrintLevel Normal",
    "Print[P_Basis] 2",
    "Print[P_MOs] 1",
    "Print[P_Overlap] 1",
    "Print[P_Hirshfeld] 1",
]


def _validate_orca_block_body(name: str, body: str) -> None:
    """Allow block-local value lines only, never syntax that can escape the generated block."""

    for line_number, line in enumerate(body.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        first_token = stripped.split(None, 1)[0].casefold()
        if first_token in {"end", "$new_job"} or stripped.startswith(("%", "*")):
            raise InputError(
                f"ORCA block {name!r} body contains reserved outer syntax on line {line_number}: {line!r}"
            )


def build_input(input_model: "AtomicInput", config: TaskConfig, executable: str) -> Job:
    """Generate an ORCA job from the supported QCSchema subset."""

    driver, method, basis = _validate_input_subset(input_model)
    driver_keyword = {"energy": None, "gradient": "engrad", "hessian": "freq"}[driver.lower()]

    keywords = input_model.specification.keywords
    if not isinstance(keywords, Mapping):
        raise InputError("ORCA keywords must be a mapping")
    unknown = set(keywords) - {"simple", "blocks"}
    if unknown:
        raise InputError(f"ORCA keywords contain unknown top-level keys: {sorted(unknown)!r}")

    simple = keywords.get("simple", [])
    if not isinstance(simple, list):
        raise InputError("ORCA simple keywords must be a list")
    for value in simple:
        if not isinstance(value, str) or not value.strip() or "\n" in value or "\r" in value:
            raise InputError("ORCA simple keywords must be non-empty strings without newlines")

    blocks = keywords.get("blocks", {})
    if not isinstance(blocks, Mapping):
        raise InputError("ORCA blocks must be a mapping")
    output_body: Optional[str] = None
    user_blocks: Dict[str, str] = {}
    for name, body in blocks.items():
        if not isinstance(name, str) or re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) is None:
            raise InputError(f"ORCA block name is invalid: {name!r}")
        if not isinstance(body, str):
            raise InputError(f"ORCA block {name!r} body must be a string")
        _validate_orca_block_body(name, body)
        normalized_name = name.lower()
        if normalized_name in {"pal", "maxcore"}:
            raise InputError(f"ORCA block {name!r} is reserved for TaskConfig resources")
        if normalized_name == "coords":
            raise InputError("ORCA coordinate block 'coords' is reserved for QCSchema geometry")
        if normalized_name == "output":
            output_body = body
        else:
            user_blocks[name] = body

    simple_line = [method, basis]
    if driver_keyword is not None:
        simple_line.append(driver_keyword)
    simple_line.extend(simple)
    lines = ["! " + " ".join(simple_line), "%output", *_ORCA_OUTPUT_DEFAULTS]
    if output_body is not None:
        lines.extend(output_body.splitlines())
    lines.append("end")

    for name in sorted(user_blocks, key=lambda value: (value.lower(), value)):
        lines.append(f"%{name}")
        lines.extend(user_blocks[name].splitlines())
        lines.append("end")

    maxcore = max(1, int(config.memory * 1024 / config.ncores))
    lines.extend(["%pal", f"nprocs {config.ncores}", "end", f"%MaxCore {maxcore}"])

    molecule = input_model.molecule
    charge = int(molecule.molecular_charge)
    multiplicity = int(molecule.molecular_multiplicity)
    lines.append(f"* xyz {charge} {multiplicity}")
    for symbol, coordinates in zip(molecule.symbols, molecule.geometry):
        converted = [float(coordinate) * constants.bohr2angstroms for coordinate in coordinates]
        lines.append(f"{symbol} {str(converted[0])} {str(converted[1])} {str(converted[2])}")
    lines.append("*")
    input_text = "\n".join(lines) + "\n"

    return Job(
        command=[executable, "dispatch.inp"],
        infiles={"dispatch.inp": input_text},
        outfiles=[],
        input_filename="dispatch.inp",
        output_filename="dispatch.out",
        input_text=input_text,
        executable=executable,
    )


def select_output(outputs: Mapping[str, Any]) -> str:
    """Select ORCA's captured standard output for parsing."""

    try:
        output = outputs["stdout"]
    except KeyError as exc:
        raise UnknownError("ORCA did not produce captured stdout") from exc
    if output is None:
        raise UnknownError("ORCA did not produce captured stdout")
    return str(output)


def preflight(executable: str, environment: Mapping[str, str]) -> Dict[str, str]:
    """Return an isolated ORCA child environment."""

    return dict(environment)


def probe(executable: str, environment: Mapping[str, str]) -> str:
    """Verify ORCA identity and return its normalized version."""

    probe_input = "! HF STO-3G\n* xyz 0 2\nH 0.0 0.0 0.0\n*\n"
    success, outputs = execute(
        [executable, "version.inp"],
        {"version.inp": probe_input},
        environment=dict(environment),
        timeout=30,
    )
    output = f"{outputs.get('stdout') or ''}\n{outputs.get('stderr') or ''}"
    if not success:
        raise ResourceError(f"ORCA identity/version probe failed for {executable}")
    if re.search(r"O\s+R\s+C\s+A", output) is None:
        raise ResourceError(f"Executable {executable} failed ORCA identity verification")
    if "ORCA TERMINATED NORMALLY" not in output:
        raise ResourceError(f"ORCA version probe from {executable} did not reach normal termination")

    match = re.search(r"Program\s+Version\s+([0-9]+(?:\.[0-9]+)+)", output, re.IGNORECASE)
    if match is None:
        raise ResourceError(f"Could not parse the ORCA version from executable {executable}")
    return safe_version(match.group(1))


def parser_type() -> Type[Any]:
    """Load the expected ORCA parser class lazily."""

    from cclib.parser.orcaparser import ORCA

    return ORCA


ORCA_DEFINITION = ProgramDefinition(
    selector="cclib-orca",
    executable="orca",
    minimum_version="6.0",
    input_filename="dispatch.inp",
    output_filename="dispatch.out",
    parser_name="ORCA",
    parser_type=parser_type,
    normal_termination="ORCA TERMINATED NORMALLY",
    managed_scratch_suffix=None,
    scratch_environment=None,
    generator=build_input,
    probe=probe,
    output_selector=select_output,
    preflight=preflight,
)


class ORCACCLibHarness(CCLibHarness):
    """Run ORCA and convert its output exclusively through cclib."""

    definition = ORCA_DEFINITION
    _defaults = {**CCLibHarness._defaults, "name": "cclib-orca", "node_parallel": True}
