"""Shared cclib-backed program harnesses.

External cclib modules are intentionally imported only by the loader below so that
cclib remains an optional dependency of QCEngine.
"""

import math
import os
import re
import stat
import sys
import tempfile
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Dict, List, Literal, Mapping, Optional, Tuple, Type

from qcelemental import constants
from qcelemental.util import parse_version, safe_version, which

from ..config import TaskConfig
from ..exceptions import InputError, ResourceError, UnknownError
from ..util import execute, temporary_directory
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models.v2 import AtomicInput, AtomicResult


@dataclass(frozen=True)
class _CCLibAPI:
    version: str
    ccData: Type[Any]
    QCSchemaWriter: Type[Any]
    ccread: Any
    QChem: Type[Any]
    ORCA: Type[Any]

    @property
    def ccopen(self) -> Any:
        """Expose the auto-detecting opener while preserving the reviewed lazy API shape."""

        return self.ccread


_cclib_compatibility_cache: Optional[Tuple[bool, str]] = None


@dataclass(frozen=True)
class _Job:
    command: List[str]
    infiles: Dict[str, str]
    outfiles: List[str]
    input_filename: str
    output_filename: str
    input_text: str
    executable: str


@dataclass(frozen=True)
class _ExecutionResult:
    process_success: bool
    executable: str
    input_filename: str
    output_filename: str
    input_text: str
    output_text: str
    stdout: str
    stderr: str


def _validate_input_subset(input_model: "AtomicInput") -> Tuple[str, str, str]:
    """Validate and return the unmodified supported QCSchema request fields."""

    driver_value = input_model.specification.driver
    driver = driver_value.value if hasattr(driver_value, "value") else str(driver_value)
    method = input_model.specification.model.method
    basis = input_model.specification.model.basis

    if driver.lower() not in {"energy", "gradient", "hessian"}:
        raise InputError(f"Unsupported driver for CCLibHarness: {driver}")
    if method.lower() not in {"hf", "b3lyp", "bp86", "mp2", "ccsd"}:
        raise InputError(f"Unsupported method for CCLibHarness: {method}")
    if not isinstance(basis, str) or not basis.strip():
        raise InputError("CCLibHarness basis must be a non-empty string")
    if not all(bool(real) for real in input_model.molecule.real):
        raise InputError("CCLibHarness requires all atoms to be real; ghost atoms are unsupported")

    return driver, method, basis


def _load_cclib_api() -> _CCLibAPI:
    """Load the optional cclib interfaces only when a cclib harness is checked."""

    import cclib
    from cclib.io import ccopen
    from cclib.io.qcschemawriter import QCSchemaWriter
    from cclib.parser.data import ccData
    from cclib.parser.orcaparser import ORCA
    from cclib.parser.qchemparser import QChem

    return _CCLibAPI(
        version=cclib.__version__,
        ccData=ccData,
        QCSchemaWriter=QCSchemaWriter,
        ccread=ccopen,
        QChem=QChem,
        ORCA=ORCA,
    )


def _synthetic_ccdata(ccdata_type: Type[Any]) -> Any:
    """Construct the smallest writer input that exercises required cclib fields."""

    return ccdata_type(
        {
            "atomcoords": [[[1.0, 0.0, 0.0]]],
            "atomnos": [2],
            "atomcharges": {"mulliken": [0.0]},
            "charge": 0,
            "homos": [0],
            "metadata": {
                "package": "Synthetic",
                "package_version": "1.0",
                "methods": ["HF"],
                "basis_set": "sto-3g",
                "success": True,
            },
            "mult": 1,
            "natom": 1,
            "nbasis": 1,
            "nmo": 1,
            "scfenergies": [-24.6],
            "scftargets": [[[1.0e-6]]],
            "scfvalues": [[[1.0e-4], [1.0e-7]]],
        }
    )


def _validate_v1_atomic_result(output: Dict[str, Any]) -> Any:
    """Validate writer output against QCElemental's explicit v1 result model."""

    try:
        from qcelemental.models.v1 import AtomicResult

        return AtomicResult(**output)
    except RuntimeError:
        if sys.version_info < (3, 14):
            raise
        from qcelemental.models._v1v2 import AtomicResult

        return AtomicResult(**output)


def _check_cclib_compatibility() -> Tuple[bool, str]:
    """Probe and cache the cclib writer behavior required by this harness."""

    global _cclib_compatibility_cache
    if _cclib_compatibility_cache is not None:
        return _cclib_compatibility_cache

    try:
        api = _load_cclib_api()
        output = api.QCSchemaWriter(_synthetic_ccdata(api.ccData)).as_dict(validate=False)
        result = _validate_v1_atomic_result(output)

        extras = result.extras
        for attribute in ("atomcoords", "atomcharges"):
            if attribute not in extras:
                raise ValueError(f"missing flat {attribute} extra")

        expected_bohr = 1.8897261255
        molecule_coordinate = float(result.molecule.geometry[0][0])
        extra_coordinate = float(extras["atomcoords"][0][0][0])
        if abs(molecule_coordinate - expected_bohr) > 2.0e-9 or abs(extra_coordinate - expected_bohr) > 2.0e-9:
            raise ValueError("wrong geometry units; cclib must convert Angstrom to bohr")
        if extras["atomcharges"] != {"mulliken": [0.0]}:
            raise ValueError("missing flat atomcharges values")
    except (ModuleNotFoundError, ImportError) as exc:
        _cclib_compatibility_cache = (False, f"cclib is not importable: {exc}")
    except Exception as exc:
        _cclib_compatibility_cache = (False, f"incompatible cclib QCSchema writer: {exc}")
    else:
        _cclib_compatibility_cache = (True, str(api.version))

    return _cclib_compatibility_cache


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


def _build_qchem_input(input_model: "AtomicInput", config: TaskConfig, executable: str) -> _Job:
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
    return _Job(
        command=[executable, "-nt", str(config.ncores), "dispatch.in", "dispatch.out"],
        infiles={"dispatch.in": input_text},
        outfiles=["dispatch.out"],
        input_filename="dispatch.in",
        output_filename="dispatch.out",
        input_text=input_text,
        executable=executable,
    )


_ORCA_OUTPUT_DEFAULTS = [
    "PrintLevel Normal",
    "Print[P_Basis] 2",
    "Print[P_MOs] 1",
    "Print[P_Overlap] 1",
    "Print[P_Hirshfeld] 1",
]


def _build_orca_input(input_model: "AtomicInput", config: TaskConfig, executable: str) -> _Job:
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

    return _Job(
        command=[executable, "dispatch.inp"],
        infiles={"dispatch.inp": input_text},
        outfiles=[],
        input_filename="dispatch.inp",
        output_filename="dispatch.out",
        input_text=input_text,
        executable=executable,
    )


def _select_qchem_output(outputs: Mapping[str, Any]) -> str:
    try:
        output = outputs["dispatch.out"]
    except KeyError as exc:
        raise UnknownError("Q-Chem did not produce dispatch.out") from exc
    if output is None:
        raise UnknownError("Q-Chem did not produce dispatch.out")
    return str(output)


def _select_orca_output(outputs: Mapping[str, Any]) -> str:
    try:
        output = outputs["stdout"]
    except KeyError as exc:
        raise UnknownError("ORCA did not produce captured stdout") from exc
    if output is None:
        raise UnknownError("ORCA did not produce captured stdout")
    return str(output)


def _path_has_mode(path: str, mode: int) -> bool:
    try:
        return bool(os.stat(path).st_mode & mode) and os.access(path, os.R_OK if mode == stat.S_IRUSR else os.X_OK)
    except OSError:
        return False


def _preflight_none(executable: str, environment: Mapping[str, str]) -> Dict[str, str]:
    return dict(environment)


def _preflight_qchem(executable: str, environment: Mapping[str, str]) -> Dict[str, str]:
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


def _probe_qchem(executable: str, environment: Mapping[str, str]) -> str:
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


def _probe_orca(executable: str, environment: Mapping[str, str]) -> str:
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


@dataclass(frozen=True)
class _ProgramDefinition:
    selector: str
    executable: str
    minimum_version: str
    input_filename: str
    output_filename: str
    expected_parser: str
    normal_termination: str
    generator: Callable[..., Any]
    probe: Callable[[str, Mapping[str, str]], str]
    output_selector: Callable[[Mapping[str, Any]], str]
    preflight: Callable[[str, Mapping[str, str]], Dict[str, str]]

    @property
    def parser_name(self) -> str:
        """Return the reviewed parser identity under its result-metadata name."""

        return self.expected_parser


_PROGRAM_DEFINITIONS: Mapping[str, _ProgramDefinition] = {
    "qchem": _ProgramDefinition(
        selector="cclib-qchem",
        executable="qchem",
        minimum_version="5.1",
        input_filename="dispatch.in",
        output_filename="dispatch.out",
        expected_parser="QChem",
        normal_termination="Thank you very much for using Q-Chem",
        generator=_build_qchem_input,
        probe=_probe_qchem,
        output_selector=_select_qchem_output,
        preflight=_preflight_qchem,
    ),
    "orca": _ProgramDefinition(
        selector="cclib-orca",
        executable="orca",
        minimum_version="6.0",
        input_filename="dispatch.inp",
        output_filename="dispatch.out",
        expected_parser="ORCA",
        normal_termination="ORCA TERMINATED NORMALLY",
        generator=_build_orca_input,
        probe=_probe_orca,
        output_selector=_select_orca_output,
        preflight=_preflight_none,
    ),
}


def _diagnostic_tail(text: str, max_lines: int = 40, max_chars: int = 4000) -> str:
    """Return the bounded trailing portion used in execution diagnostics."""

    line_bounded = "\n".join(text.splitlines()[-max_lines:])
    return line_bounded[-max_chars:]


def _missing_qchem_environment_variable(diagnostic: str) -> Optional[str]:
    patterns = (
        r"undefined\s+environment\s+variable\s*[:=]?\s*['\"`]?\$?\{?([A-Za-z_][A-Za-z0-9_]*)",
        r"['\"`]?\b([A-Za-z_][A-Za-z0-9_]*)['\"`]?\s*:\s*(?:undefined|unbound)\s+variable",
        r"environment\s+variable\s+['\"`]?\$?\{?([A-Za-z_][A-Za-z0-9_]*)['\"`}]?\s+"
        r"(?:is\s+)?(?:undefined|not\s+(?:defined|set)|must\s+be\s+defined)",
    )
    for pattern in patterns:
        match = re.search(pattern, diagnostic, re.IGNORECASE)
        if match is not None:
            return match.group(1)
    return None


def _raise_execution_failure(
    definition: _ProgramDefinition,
    job: _Job,
    stage: str,
    diagnostic: str,
    cause: Optional[BaseException] = None,
) -> None:
    """Raise one consistently formatted, stage-aware execution failure."""

    missing_variable = (
        _missing_qchem_environment_variable(diagnostic) if definition.selector == "cclib-qchem" else None
    )
    is_environment_failure = definition.selector == "cclib-qchem" and (
        missing_variable is not None
        or re.search(
            r"undefined\s+(?:environment\s+)?variable|environment\s+variable.*"
            r"(?:undefined|not\s+(?:defined|set)|must\s+be\s+defined)",
            diagnostic,
            re.IGNORECASE,
        )
        is not None
    )
    is_license_failure = any(
        marker in diagnostic.lower()
        for marker in ("flexnet", "license checkout", "unable to validate license")
    )
    error_type = ResourceError if is_environment_failure or is_license_failure else UnknownError

    detail = ""
    if missing_variable is not None:
        detail = f"; undefined Q-Chem environment variable {missing_variable}"
    elif is_environment_failure:
        detail = "; undefined Q-Chem environment variable"
    elif is_license_failure:
        detail = "; license resource unavailable"
    message = (
        f"{definition.selector} failed for resolved executable {job.executable} "
        f"during {stage} stage{detail}.\nDiagnostic tail:\n{_diagnostic_tail(diagnostic)}"
    )
    if cause is None:
        raise error_type(message)
    raise error_type(message) from cause


def _execution_diagnostic(outputs: Mapping[str, Any], stdout: str, stderr: str) -> str:
    parts = [str(value) for value in outputs.values() if value not in (None, "")]
    for stream in (stdout, stderr):
        if stream and stream not in parts:
            parts.append(stream)
    return "\n".join(parts)


def _execute_job(definition: _ProgramDefinition, job: _Job, config: TaskConfig) -> _ExecutionResult:
    """Execute a generated job with QCEngine's managed scratch utilities."""

    environment = definition.preflight(job.executable, os.environ.copy())

    def run(scratch_directory: Optional[str]) -> _ExecutionResult:
        try:
            process_success, process = execute(
                job.command,
                infiles=job.infiles,
                outfiles=job.outfiles,
                scratch_directory=scratch_directory,
                scratch_messy=config.scratch_messy,
                environment=environment,
            )
        except Exception as exc:
            _raise_execution_failure(definition, job, "execution", str(exc), exc)

        stdout = str(process.get("stdout") or "")
        stderr = str(process.get("stderr") or "")
        outputs = dict(process.get("outfiles") or {})
        diagnostic = _execution_diagnostic(outputs, stdout, stderr)
        if not process_success:
            _raise_execution_failure(definition, job, "execution", diagnostic)

        selection_inputs = dict(outputs)
        if "stdout" in process:
            selection_inputs["stdout"] = process["stdout"]
        if "stderr" in process:
            selection_inputs["stderr"] = process["stderr"]
        try:
            output_text = definition.output_selector(selection_inputs)
        except Exception as exc:
            _raise_execution_failure(definition, job, "output selection", diagnostic, exc)

        if definition.normal_termination not in output_text:
            _raise_execution_failure(definition, job, "termination", output_text)

        return _ExecutionResult(
            process_success=process_success,
            executable=job.executable,
            input_filename=job.input_filename,
            output_filename=job.output_filename,
            input_text=job.input_text,
            output_text=output_text,
            stdout=stdout,
            stderr=stderr,
        )

    if definition.selector == "cclib-qchem":
        with temporary_directory(
            parent=config.scratch_directory,
            suffix="_cclib_qchem_scratch",
            messy=config.scratch_messy,
        ) as qcscratch:
            environment["QCSCRATCH"] = str(qcscratch)
            return run(str(qcscratch))
    return run(config.scratch_directory)


_METHOD_ALIASES = {
    "rhf": "hf",
    "uhf": "hf",
    "rmp2": "mp2",
    "ump2": "mp2",
    "rccsd": "ccsd",
    "uccsd": "ccsd",
}


def _raise_conversion_failure(
    definition: _ProgramDefinition,
    execution: _ExecutionResult,
    stage: str,
    cause: BaseException,
) -> None:
    """Raise one bounded, stage-aware parsing/conversion failure."""

    diagnostic = _diagnostic_tail(f"{execution.output_text}\n{cause}")
    raise UnknownError(
        f"{definition.selector} failed for resolved executable {execution.executable} "
        f"during {stage} stage.\nDiagnostic tail:\n{diagnostic}"
    ) from cause


def _native_files(input_model: "AtomicInput", execution: _ExecutionResult) -> Dict[str, str]:
    protocol = input_model.specification.protocols.native_files
    value = protocol.value if hasattr(protocol, "value") else str(protocol)
    if value == "none":
        return {}
    files = {"input": execution.input_text}
    if value == "all":
        files[execution.output_filename] = execution.output_text
    return files


def _parse_and_convert(
    definition: _ProgramDefinition,
    execution: _ExecutionResult,
    input_model: "AtomicInput",
) -> "AtomicResult":
    """Parse complete program output, validate QCSchema v1, and convert to v2."""

    try:
        api = _load_cclib_api()
    except Exception as exc:
        _raise_conversion_failure(definition, execution, "parser setup", exc)

    parser = None
    with tempfile.NamedTemporaryFile(mode="w", suffix=".out", encoding="utf-8") as temporary:
        temporary.write(execution.output_text)
        temporary.flush()
        # tempfile's wrapper is seekable but not itself an iterator, while
        # cclib's FileWrapper requires both. Reopen its named file as a normal
        # text stream so auto-detection sees the exact complete output.
        with open(temporary.name, encoding="utf-8") as source:
            try:
                parser = api.ccopen(source)
            except Exception as exc:
                _raise_conversion_failure(definition, execution, "parser auto-detection", exc)
            if parser is None:
                _raise_conversion_failure(
                    definition,
                    execution,
                    "parser auto-detection",
                    ValueError("cclib ccopen did not detect a parser"),
                )

            try:
                expected_parser = getattr(api, definition.parser_name)
                if type(parser) is not expected_parser:
                    _raise_conversion_failure(
                        definition,
                        execution,
                        "parser identity",
                        ValueError(
                            f"cclib selected {type(parser).__name__}; expected {definition.parser_name}"
                        ),
                    )

                try:
                    parsed = parser.parse()
                except Exception as exc:
                    _raise_conversion_failure(definition, execution, "parser parse", exc)
            finally:
                parser_input = getattr(parser, "inputfile", None)
                if parser_input is not None:
                    try:
                        parser_input.close()
                    except Exception:
                        pass

    metadata = getattr(parsed, "metadata", None)
    if not isinstance(metadata, Mapping) or metadata.get("success") is not True:
        _raise_conversion_failure(
            definition,
            execution,
            "parser result validation",
            ValueError("cclib parser result is incomplete: metadata['success'] is not true"),
        )

    try:
        writer_output = api.QCSchemaWriter(parsed).as_dict(validate=False)
    except Exception as exc:
        _raise_conversion_failure(definition, execution, "QCSchema writer", exc)

    required_fields = {
        "schema_name",
        "schema_version",
        "molecule",
        "provenance",
        "success",
        "extras",
        "driver",
        "model",
        "properties",
        "return_result",
    }
    missing_fields = required_fields - set(writer_output) if isinstance(writer_output, Mapping) else required_fields
    if missing_fields:
        _raise_conversion_failure(
            definition,
            execution,
            "QCSchema writer output",
            ValueError(f"writer output is missing required fields: {sorted(missing_fields)!r}"),
        )

    try:
        output = dict(writer_output)
        writer_extras = output.get("extras")
        if not isinstance(writer_extras, Mapping):
            raise TypeError("QCSchema writer extras must be a mapping")
        extras = dict(writer_extras)
        if "cclib_harness" in extras:
            raise ValueError("QCSchema writer extras already contain reserved key 'cclib_harness'")
        extras["cclib_harness"] = {
            "selector": definition.selector,
            "cclib_version": api.version,
            "parser": definition.parser_name,
            "executable": execution.executable,
        }
        output["extras"] = extras
        output["stdout"] = execution.output_text
        output["stderr"] = execution.stderr or None
        native_protocol = input_model.specification.protocols.native_files
        native_protocol_value = native_protocol.value if hasattr(native_protocol, "value") else str(native_protocol)
        output["protocols"] = {"native_files": native_protocol_value, "stdout": True}
        native_files = _native_files(input_model, execution)
        output["native_files"] = native_files
    except Exception as exc:
        _raise_conversion_failure(definition, execution, "QCSchema writer augmentation", exc)

    requested_driver, requested_method, requested_basis = _validate_input_subset(input_model)
    writer_model = output.get("model") if isinstance(output.get("model"), Mapping) else {}
    parsed_driver = output.get("driver")
    parsed_method = writer_model.get("method")
    parsed_basis = writer_model.get("basis")
    identities = (
        ("driver", requested_driver.casefold(), str(parsed_driver).casefold()),
        (
            "method",
            _METHOD_ALIASES.get(requested_method.casefold(), requested_method.casefold()),
            _METHOD_ALIASES.get(str(parsed_method).casefold(), str(parsed_method).casefold()),
        ),
        ("basis", requested_basis.casefold(), str(parsed_basis).casefold()),
    )
    for field, requested, parsed_value in identities:
        if requested != parsed_value:
            _raise_conversion_failure(
                definition,
                execution,
                f"parsed {field} mismatch",
                ValueError(f"requested {field} {requested!r}, parsed {parsed_value!r}"),
            )

    native_files_in_v1 = True
    try:
        result_v1 = _validate_v1_atomic_result(output)
    except Exception as exc:
        # Some older explicit v1 models predate native_files. Only use the
        # documented fallback when removing that field alone fixes validation.
        without_native = dict(output)
        without_native.pop("native_files")
        try:
            result_v1 = _validate_v1_atomic_result(without_native)
        except Exception:
            _raise_conversion_failure(definition, execution, "QCSchema v1 validation", exc)
        native_files_in_v1 = False
        if native_files:
            result_v1.extras["cclib_harness"]["native_input"] = execution.input_text

    try:
        result_v2 = result_v1.convert_v(2, external_input_data=input_model)
    except Exception as exc:
        _raise_conversion_failure(definition, execution, "QCSchema v2 conversion", exc)

    if not native_files_in_v1:
        # The fallback is represented only in cclib_harness metadata.
        assert result_v2.native_files is None or not result_v2.native_files
    return result_v2


def _probe_executable(
    harness: "CCLibHarness", executable: str, environment: Optional[Mapping[str, str]] = None
) -> str:
    if executable in harness.version_cache:
        return harness.version_cache[executable]

    definition = _PROGRAM_DEFINITIONS[harness.program]
    child_environment = dict(os.environ if environment is None else environment)
    try:
        version = definition.probe(executable, child_environment)
    except ResourceError:
        raise
    except Exception as exc:
        raise ResourceError(f"Failed to probe {definition.selector} executable {executable}: {exc}") from exc

    if parse_version(version) < parse_version(definition.minimum_version):
        raise ResourceError(
            f"{definition.selector} requires {harness.program.upper()} version {definition.minimum_version} or newer; "
            f"found {version} at {executable}"
        )

    harness.version_cache[executable] = version
    return version


class CCLibHarness(ProgramHarness):
    """A cclib-backed harness configured for one external QC program."""

    program: Literal["qchem", "orca"]

    _defaults: ClassVar[Dict[str, Any]] = {
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "managed_memory": True,
    }
    version_cache: ClassVar[Dict[str, str]] = {}

    def __init__(self, **kwargs: Any):
        kwargs["node_parallel"] = kwargs.get("program") == "orca"
        super().__init__(**kwargs)

    def found(self, raise_error: bool = False) -> bool:
        definition = _PROGRAM_DEFINITIONS[self.program]
        try:
            compatible, version_or_message = _check_cclib_compatibility()
            if not compatible:
                raise ResourceError(f"{definition.selector} requires compatible cclib: {version_or_message}")

            executable = which(definition.executable)
            if executable is None:
                raise ResourceError(
                    f"{definition.selector} executable '{definition.executable}' was not found on PATH"
                )

            environment = definition.preflight(executable, os.environ.copy())
            _probe_executable(self, executable, environment)
            return True
        except ResourceError:
            if raise_error:
                raise
            return False
        except Exception as exc:
            if raise_error:
                raise ResourceError(f"{definition.selector} availability probe failed: {exc}") from exc
            return False

    def get_version(self) -> str:
        definition = _PROGRAM_DEFINITIONS[self.program]
        executable = which(definition.executable)
        if executable is None:
            raise ResourceError(f"{definition.selector} executable '{definition.executable}' was not found on PATH")
        environment = definition.preflight(executable, os.environ.copy())
        return _probe_executable(self, executable, environment)

    def build_input(
        self, input_model: "AtomicInput", config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        definition = _PROGRAM_DEFINITIONS[self.program]
        executable = which(definition.executable)
        if executable is None:
            raise ResourceError(f"{definition.selector} executable '{definition.executable}' was not found on PATH")
        job = definition.generator(input_model, config, executable)
        return {
            "commands": job.command,
            "infiles": job.infiles,
            "outfiles": job.outfiles,
            "scratch_directory": config.scratch_directory,
        }

    def compute(self, input_data: "AtomicInput", config: TaskConfig) -> "AtomicResult":
        definition = _PROGRAM_DEFINITIONS[self.program]
        _validate_input_subset(input_data)
        executable = which(definition.executable)
        if executable is None:
            raise ResourceError(
                f"{definition.selector} executable '{definition.executable}' was not found on PATH"
            )
        job = definition.generator(input_data, config, executable)
        execution = _execute_job(definition, job, config)
        return _parse_and_convert(definition, execution, input_data)
