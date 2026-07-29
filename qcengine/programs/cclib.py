"""Shared cclib-backed program harnesses.

External cclib modules are intentionally imported only by the loader below so that
cclib remains an optional dependency of QCEngine.
"""

import os
import re
import stat
import sys
import tempfile
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, ClassVar, Dict, Literal, Mapping, Optional, Tuple, Type

from qcelemental.util import parse_version, safe_version, which

from ..config import TaskConfig
from ..exceptions import InputError, ResourceError, UnknownError
from ..util import execute
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


_cclib_compatibility_cache: Optional[Tuple[bool, str]] = None


def _load_cclib_api() -> _CCLibAPI:
    """Load the optional cclib interfaces only when a cclib harness is checked."""

    import cclib
    from cclib.io import ccread
    from cclib.io.qcschemawriter import QCSchemaWriter
    from cclib.parser.data import ccData
    from cclib.parser.orcaparser import ORCA
    from cclib.parser.qchemparser import QChem

    return _CCLibAPI(
        version=cclib.__version__,
        ccData=ccData,
        QCSchemaWriter=QCSchemaWriter,
        ccread=ccread,
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


def _require_generation(*args: Any, **kwargs: Any) -> Any:
    """Mark the task boundary before native input generation is added."""

    raise InputError("Native input generation is not part of the availability probe")


def _select_qchem_output(outputs: Mapping[str, Any]) -> str:
    try:
        return str(outputs["dispatch.out"])
    except KeyError as exc:
        raise UnknownError("Q-Chem did not produce dispatch.out") from exc


def _select_orca_output(outputs: Mapping[str, Any]) -> str:
    try:
        return str(outputs["stdout"])
    except KeyError as exc:
        raise UnknownError("ORCA did not produce captured stdout") from exc


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


_PROGRAM_DEFINITIONS: Mapping[str, _ProgramDefinition] = {
    "qchem": _ProgramDefinition(
        selector="cclib-qchem",
        executable="qchem",
        minimum_version="5.1",
        input_filename="dispatch.in",
        output_filename="dispatch.out",
        expected_parser="QChem",
        normal_termination="Thank you very much for using Q-Chem",
        generator=_require_generation,
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
        generator=_require_generation,
        probe=_probe_orca,
        output_selector=_select_orca_output,
        preflight=_preflight_none,
    ),
}


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

    def compute(self, input_data: "AtomicInput", config: TaskConfig) -> "AtomicResult":
        raise NotImplementedError
