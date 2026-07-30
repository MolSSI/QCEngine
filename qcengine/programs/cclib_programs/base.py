"""Shared cclib-backed program harnesses.

External cclib modules are intentionally imported only by the loader below so that
cclib remains an optional dependency of QCEngine.
"""

import os
import re
import sys
import tempfile
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Dict, List, Mapping, Optional, Tuple, Type

from qcelemental.util import parse_version, which

from ...config import TaskConfig
from ...exceptions import InputError, ResourceError, UnknownError
from ...util import execute, temporary_directory
from ..model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models.v2 import AtomicInput, AtomicResult


@dataclass(frozen=True)
class CCLibAPI:
    """Lazily imported cclib interfaces used with real parser output."""

    version: str
    QCSchemaWriter: Type[Any]
    ccopen: Any


@dataclass(frozen=True)
class Job:
    """A complete native-program execution request."""

    command: List[str]
    infiles: Dict[str, str]
    outfiles: List[str]
    input_filename: str
    output_filename: str
    input_text: str
    executable: str


@dataclass(frozen=True)
class ExecutionResult:
    """Native-program output retained for cclib conversion."""

    process_success: bool
    executable: str
    input_filename: str
    output_filename: str
    input_text: str
    output_text: str
    stdout: str
    stderr: str


def _validate_native_token(field: str, value: Any) -> str:
    """Return one trimmed printable native token without structural whitespace."""

    if not isinstance(value, str):
        raise InputError(f"CCLibHarness {field} must be a non-empty string")
    if any(not character.isprintable() for character in value):
        raise InputError(f"CCLibHarness {field} must not contain control characters")
    normalized = value.strip()
    if not normalized:
        raise InputError(f"CCLibHarness {field} must be a non-empty string")
    if any(character.isspace() for character in normalized):
        raise InputError(f"CCLibHarness {field} must be exactly one native token")
    return normalized


def _input_fields(input_model: "AtomicInput") -> Tuple[str, str, str]:
    """Return native driver, method, and basis after shared structural validation."""

    driver_value = input_model.specification.driver
    driver = driver_value.value if hasattr(driver_value, "value") else str(driver_value)
    method = _validate_native_token("method", input_model.specification.model.method)
    basis = _validate_native_token("basis", input_model.specification.model.basis)

    if not all(bool(real) for real in input_model.molecule.real):
        raise InputError("CCLibHarness requires all atoms to be real; ghost atoms are unsupported")

    return driver, method, basis


def _load_cclib_api() -> CCLibAPI:
    """Import the optional cclib writer and auto-detecting opener on demand."""
    import cclib
    from cclib.io import ccopen
    from cclib.io.qcschemawriter import QCSchemaWriter

    return CCLibAPI(version=cclib.__version__, QCSchemaWriter=QCSchemaWriter, ccopen=ccopen)


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


@dataclass(frozen=True)
class ProgramDefinition:
    """Callbacks and identities required by one cclib-backed native program."""

    selector: str
    executable: str
    minimum_version: str
    input_filename: str
    output_filename: str
    parser_name: str
    parser_type: Callable[[], Type[Any]]
    normal_termination: str
    managed_scratch_suffix: Optional[str]
    scratch_environment: Optional[Callable[[Dict[str, str], str], None]]
    generator: Callable[["AtomicInput", TaskConfig, str], Job]
    probe: Callable[[str, Mapping[str, str]], str]
    output_selector: Callable[[Mapping[str, Any]], str]
    preflight: Callable[[str, Mapping[str, str]], Dict[str, str]]


def _diagnostic_tail(text: str, max_lines: int = 40, max_chars: int = 4000) -> str:
    """Return the bounded trailing portion used in execution diagnostics."""

    line_bounded = "\n".join(text.splitlines()[-max_lines:])
    return line_bounded[-max_chars:]


def _missing_environment_variable(diagnostic: str) -> Optional[str]:
    """Return a missing environment variable named by a native diagnostic."""
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
    definition: ProgramDefinition,
    job: Job,
    stage: str,
    diagnostic: str,
    cause: Optional[BaseException] = None,
) -> None:
    """Raise one consistently formatted, stage-aware execution failure."""

    missing_variable = _missing_environment_variable(diagnostic)
    is_environment_failure = missing_variable is not None or (
        re.search(
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
        detail = f"; undefined environment variable {missing_variable}"
    elif is_environment_failure:
        detail = "; undefined environment variable"
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
    """Combine distinct native output channels for failure reporting."""

    parts = [str(value) for value in outputs.values() if value not in (None, "")]
    for stream in (stdout, stderr):
        if stream and stream not in parts:
            parts.append(stream)
    return "\n".join(parts)


def _execute_job(definition: ProgramDefinition, job: Job, config: TaskConfig) -> ExecutionResult:
    """Execute a generated job with QCEngine's managed scratch utilities."""

    environment = definition.preflight(job.executable, os.environ.copy())

    def run(scratch_directory: Optional[str]) -> ExecutionResult:
        """Execute once within the selected scratch directory."""

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

        return ExecutionResult(
            process_success=process_success,
            executable=job.executable,
            input_filename=job.input_filename,
            output_filename=job.output_filename,
            input_text=job.input_text,
            output_text=output_text,
            stdout=stdout,
            stderr=stderr,
        )

    if definition.managed_scratch_suffix is not None:
        with temporary_directory(
            parent=config.scratch_directory,
            suffix=definition.managed_scratch_suffix,
            messy=config.scratch_messy,
        ) as managed_scratch:
            if definition.scratch_environment is not None:
                definition.scratch_environment(environment, str(managed_scratch))
            return run(str(managed_scratch))
    return run(config.scratch_directory)


def _raise_conversion_failure(
    definition: ProgramDefinition,
    execution: ExecutionResult,
    stage: str,
    cause: BaseException,
) -> None:
    """Raise one bounded, stage-aware parsing/conversion failure."""

    diagnostic = _diagnostic_tail(f"{execution.output_text}\n{cause}")
    raise UnknownError(
        f"{definition.selector} failed for resolved executable {execution.executable} "
        f"during {stage} stage.\nDiagnostic tail:\n{diagnostic}"
    ) from cause


def _native_files(input_model: "AtomicInput", execution: ExecutionResult) -> Dict[str, str]:
    """Select native files according to the requested protocol."""

    protocol = input_model.specification.protocols.native_files
    value = protocol.value if hasattr(protocol, "value") else str(protocol)
    if value == "none":
        return {}
    files = {"input": execution.input_text}
    if value == "all":
        files[execution.output_filename] = execution.output_text
    return files


def _parse_and_convert(
    definition: ProgramDefinition,
    execution: ExecutionResult,
    input_model: "AtomicInput",
) -> "AtomicResult":
    """Parse complete program output, validate QCSchema v1, and convert to v2."""

    try:
        api = _load_cclib_api()
    except Exception as exc:
        _raise_conversion_failure(definition, execution, "parser setup", exc)

    parser = None
    temporary_path: Optional[str] = None
    try:
        # cclib's FileWrapper requires an iterable stream. Create the named
        # output securely, finish writing and close it, then reopen read-only
        # so the lifecycle is portable to Windows as well as POSIX.
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".out", encoding="utf-8", delete=False
        ) as temporary:
            temporary_path = temporary.name
            temporary.write(execution.output_text)

        try:
            source = open(temporary_path, encoding="utf-8")
        except Exception as exc:
            _raise_conversion_failure(definition, execution, "parser auto-detection", exc)

        with source:
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
                try:
                    expected_parser = definition.parser_type()
                except Exception as exc:
                    _raise_conversion_failure(definition, execution, "parser type loading", exc)

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
    finally:
        if temporary_path is not None:
            try:
                os.unlink(temporary_path)
            except FileNotFoundError:
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

    requested_driver, _, requested_basis = _input_fields(input_model)
    writer_model = output.get("model") if isinstance(output.get("model"), Mapping) else {}
    parsed_driver = output.get("driver")
    parsed_basis = writer_model.get("basis")
    identities = (
        ("driver", requested_driver.casefold(), str(parsed_driver).casefold()),
        ("basis", requested_basis.strip().casefold(), str(parsed_basis).strip().casefold()),
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
    """Validate and cache the concrete harness executable version."""

    if executable in harness.version_cache:
        return harness.version_cache[executable]

    definition = harness.definition
    child_environment = dict(os.environ if environment is None else environment)
    try:
        version = definition.probe(executable, child_environment)
    except ResourceError:
        raise
    except Exception as exc:
        raise ResourceError(f"Failed to probe {definition.selector} executable {executable}: {exc}") from exc

    if parse_version(version) < parse_version(definition.minimum_version):
        raise ResourceError(
            f"{definition.selector} requires {definition.executable.upper()} version "
            f"{definition.minimum_version} or newer; found {version} at {executable}"
        )

    harness.version_cache[executable] = version
    return version


class CCLibHarness(ProgramHarness):
    """Shared execution and cclib conversion for a concrete native program."""

    definition: ClassVar[ProgramDefinition]
    _defaults: ClassVar[Dict[str, Any]] = {
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "managed_memory": True,
    }
    version_cache: ClassVar[Dict[str, str]] = {}

    def found(self, raise_error: bool = False) -> bool:
        """Return whether cclib and the configured native program are available."""

        definition = self.definition
        try:
            _load_cclib_api()
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
        """Return the validated native executable version."""

        definition = self.definition
        executable = which(definition.executable)
        if executable is None:
            raise ResourceError(f"{definition.selector} executable '{definition.executable}' was not found on PATH")
        environment = definition.preflight(executable, os.environ.copy())
        return _probe_executable(self, executable, environment)

    def build_input(
        self, input_model: "AtomicInput", config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        """Build the native command and file payload for an atomic input."""

        definition = self.definition
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
        """Execute the native program and convert its output through cclib."""

        definition = self.definition
        _input_fields(input_data)
        executable = which(definition.executable)
        if executable is None:
            raise ResourceError(
                f"{definition.selector} executable '{definition.executable}' was not found on PATH"
            )
        job = definition.generator(input_data, config, executable)
        execution = _execute_job(definition, job, config)
        return _parse_and_convert(definition, execution, input_data)
