# CCLib Harness Review Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor the cclib-backed Q-Chem and ORCA harnesses into an extensible package, remove synthetic runtime compatibility machinery, pass methods through to native programs, align driver support with verified behavior, and consolidate tests.

**Architecture:** A shared `CCLibHarness` in `cclib_programs/base.py` owns lazy cclib loading, managed execution, parsing, and QCSchema conversion. `QChemCCLibHarness` and `ORCACCLibHarness` subclasses own immutable definitions whose callbacks live in `cclib_qchem.py` and `cclib_orca.py`; tests live beside the package under `cclib_programs/tests/`.

**Tech Stack:** Python 3.10+, QCEngine `ProgramHarness`, QCElemental QCSchema v1/v2, cclib, Pydantic, pytest, Sphinx.

## Global Constraints

- Keep cclib optional and lazily imported; importing QCEngine must not import external cclib modules.
- Do not modify native `qcengine/programs/qchem.py` or the cclib checkout.
- Do not add cclib to QCEngine core dependencies.
- Q-Chem supports energy, gradient, and Hessian; ORCA supports energy only.
- Native method strings pass through without a QCEngine allowlist or alias table.
- Retain writer-owned extras collision protection, parser identity checks, portable temporary-file cleanup, bounded stage-aware errors, and native-file protocols.
- Tests must live at `qcengine/programs/cclib_programs/tests/test_cclib.py` and be reduced rather than merely relocated.
- Version tests must use availability-gated installed software, not fabricated executable output.
- Reusable source and unit tests must not contain machine-local executable or setup paths.
- Keep `qchem_water_mp2.py` and `orca_water_ccsd.py` deleted as established by commit `71e6531`; remove their obsolete tests.

## File Structure

- Create `qcengine/programs/cclib_programs/__init__.py`: public package exports.
- Create `qcengine/programs/cclib_programs/base.py`: shared harness, execution, parser, and QCSchema conversion.
- Create `qcengine/programs/cclib_programs/cclib_qchem.py`: Q-Chem input, preflight, probe, and subclass.
- Create `qcengine/programs/cclib_programs/cclib_orca.py`: ORCA input, probe, and subclass.
- Create `qcengine/programs/cclib_programs/tests/__init__.py`: test package marker.
- Move and reduce `qcengine/programs/tests/test_cclib.py` to `qcengine/programs/cclib_programs/tests/test_cclib.py`.
- Delete `qcengine/programs/cclib.py` after imports and tests use the package.
- Modify `qcengine/programs/base.py`: register the two concrete subclasses.
- Modify `qcengine/tests/test_harness_canonical.py`: add canonical cclib selectors.
- Modify `docs/source/programs_cclib.rst`: correct support and explain extension structure.

---

### Task 1: Remove synthetic runtime compatibility machinery

**Files:**
- Modify: `qcengine/programs/cclib.py:28-184,907-955`
- Modify: `qcengine/programs/tests/test_cclib.py:1-367`

**Interfaces:**
- Consumes: existing `_load_cclib_api()`, `CCLibHarness.found()`, and `_parse_and_convert()`.
- Produces: `_CCLibAPI(version, QCSchemaWriter, ccopen)` with no `ccData`; `_ProgramDefinition.parser_type`; availability based on lazy imports plus the real executable probe.

- [ ] **Step 1: Rewrite the compatibility-focused tests before production changes**

Delete the entire `fake_cclib_api` fixture, `clear_cclib_caches`, and these tests:

```text
test_cclib_missing_is_reported_by_compatibility_probe
test_cclib_compatibility_validates_v1_geometry_and_flat_extras
test_cclib_compatibility_rejects_incomplete_writer
test_cclib_compatibility_probe_is_cached_once_per_process
test_real_cclib_compatibility_when_installed
test_executable_identity_and_supported_version
test_executable_identity_or_version_rejection
test_unrelated_orca_executable_is_rejected
test_executable_probe_cache_is_keyed_by_resolved_path
test_get_version_returns_external_program_version
all tests and helper fixtures that import `qchem_water_mp2` or `orca_water_ccsd`
```

Update availability tests to patch `_load_cclib_api` directly and verify that `found()` calls it before executable resolution:

```python
def test_found_lazily_loads_cclib_before_resolving_executable(monkeypatch):
    calls = []
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: calls.append("cclib") or object())
    monkeypatch.setattr(cclib_harness, "which", lambda command: calls.append(command) or None)
    harness = CCLibHarness(name="cclib-orca", program="orca")

    assert harness.found() is False
    assert calls == ["cclib", "orca"]
```

- [ ] **Step 2: Run the focused tests to establish the refactor baseline**

Run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest \
  qcengine/programs/tests/test_cclib.py \
  -k 'registration or import_qcengine or found or missing_executable' -q
```

Expected: the new availability-order test fails because `found()` still routes through the synthetic compatibility probe instead of calling the lazy loader as the availability check.

- [ ] **Step 3: Simplify the lazy API and availability path**

Change the production API shape to:

```python
@dataclass(frozen=True)
class _CCLibAPI:
    """Lazily imported cclib interfaces used with real parser output."""

    version: str
    QCSchemaWriter: Type[Any]
    ccopen: Any


def _load_cclib_api() -> _CCLibAPI:
    """Import the optional cclib writer and auto-detecting opener on demand."""
    import cclib
    from cclib.io import ccopen
    from cclib.io.qcschemawriter import QCSchemaWriter

    return _CCLibAPI(version=cclib.__version__, QCSchemaWriter=QCSchemaWriter, ccopen=ccopen)
```

Delete `_synthetic_ccdata`, `_check_cclib_compatibility`, and `_cclib_compatibility_cache`. Keep `_validate_v1_atomic_result` because real conversion uses it. Add `parser_type: Callable[[], Type[Any]]` to the existing `_ProgramDefinition`, define lazy `_qchem_parser_type()` and `_orca_parser_type()` callbacks, and change `_parse_and_convert()` from `getattr(api, definition.parser_name)` to `definition.parser_type()`. Change `CCLibHarness.found()` to call `_load_cclib_api()` and then resolve/preflight/probe the executable; convert import failures to the existing `ResourceError` behavior without running a writer against fabricated data.

- [ ] **Step 4: Update conversion test fakes and run all cclib tests**

In remaining conversion fakes, remove every `ccData=`, `QChem=`, and `ORCA=` field and rename `ccread=` to `ccopen=`. Set each existing program definition's `parser_type` callback to the corresponding fake parser type in tests that require parser identity.

Run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest qcengine/programs/tests/test_cclib.py -q
```

Expected: all retained tests pass.

- [ ] **Step 5: Commit the runtime cleanup**

```bash
git add qcengine/programs/cclib.py qcengine/programs/tests/test_cclib.py
git commit -m "refactor: remove synthetic cclib compatibility probe"
```

---

### Task 2: Split the monolithic harness into the cclib program package

**Files:**
- Create: `qcengine/programs/cclib_programs/__init__.py`
- Create: `qcengine/programs/cclib_programs/base.py`
- Create: `qcengine/programs/cclib_programs/cclib_qchem.py`
- Create: `qcengine/programs/cclib_programs/cclib_orca.py`
- Create: `qcengine/programs/cclib_programs/tests/__init__.py`
- Move: `qcengine/programs/tests/test_cclib.py` → `qcengine/programs/cclib_programs/tests/test_cclib.py`
- Delete: `qcengine/programs/cclib.py`
- Modify: `qcengine/programs/base.py:8,80-82`

**Interfaces:**
- Produces: `CCLibHarness`, `ProgramDefinition`, `Job`, `ExecutionResult`, `QChemCCLibHarness`, and `ORCACCLibHarness`.
- Concrete subclasses expose `definition: ClassVar[ProgramDefinition]`; the base has no program-name field or program-specific branch.

- [ ] **Step 1: Move the test and add failing package-import assertions**

Run:

```bash
mkdir -p qcengine/programs/cclib_programs/tests
git mv qcengine/programs/tests/test_cclib.py qcengine/programs/cclib_programs/tests/test_cclib.py
touch qcengine/programs/cclib_programs/tests/__init__.py
```

At the top of the moved test use:

```python
import qcengine.programs.cclib_programs.base as cclib_base
import qcengine.programs.cclib_programs.cclib_orca as cclib_orca
import qcengine.programs.cclib_programs.cclib_qchem as cclib_qchem
from qcengine.programs.cclib_programs import ORCACCLibHarness, QChemCCLibHarness
```

Update the registration assertion to require concrete types:

```python
assert type(qcng.get_program("cclib-qchem", check=False)) is QChemCCLibHarness
assert type(qcng.get_program("cclib-orca", check=False)) is ORCACCLibHarness
```

- [ ] **Step 2: Run the import test and verify failure**

Run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest \
  qcengine/programs/cclib_programs/tests/test_cclib.py::test_registered_instances_are_frozen_and_independent -q
```

Expected: collection fails because `qcengine.programs.cclib_programs` does not exist.

- [ ] **Step 3: Create the shared base module**

Move these responsibilities from the monolith into `cclib_programs/base.py`, preserving their tested bodies and adding docstrings to every helper:

```text
_CCLibAPI → CCLibAPI
_Job → Job
_ExecutionResult → ExecutionResult
_ProgramDefinition → ProgramDefinition
_load_cclib_api
_validate_v1_atomic_result
_diagnostic_tail
_raise_execution_failure
_execution_diagnostic
_execute_job
_raise_conversion_failure
_native_files
_parse_and_convert
_probe_executable
CCLibHarness
```

Define the extension interface as:

```python
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
    generator: Callable[["AtomicInput", TaskConfig, str], Job]
    probe: Callable[[str, Mapping[str, str]], str]
    output_selector: Callable[[Mapping[str, Any]], str]
    preflight: Callable[[str, Mapping[str, str]], Dict[str, str]]
```

Define the base class without a `program` model field:

```python
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
```

Inherited `found`, `get_version`, `build_input`, and `compute` must use `self.definition` only. Replace the Q-Chem scratch branch in shared execution with an optional definition callback or definition field rather than checking a selector string. Use `ProgramDefinition.managed_scratch_suffix: Optional[str]`, set only by Q-Chem, and create managed scratch when it is non-`None`.

- [ ] **Step 4: Create the program-specific modules**

Move Q-Chem symbols into `cclib_qchem.py`, rename extension callbacks consistently, and add docstrings:

```text
QCHEM_RESERVED
_render_qchem_scalar
build_input
select_output
_path_has_mode
preflight
probe
parser_type
QCHEM_DEFINITION
QChemCCLibHarness
```

Move ORCA symbols into `cclib_orca.py`, rename extension callbacks consistently, and add docstrings:

```text
_ORCA_OUTPUT_DEFAULTS
_validate_orca_block_body
build_input
select_output
preflight
probe
parser_type
ORCA_DEFINITION
ORCACCLibHarness
```

Concrete subclasses use complete defaults:

```python
class QChemCCLibHarness(CCLibHarness):
    """Run Q-Chem and convert its output exclusively through cclib."""

    definition = QCHEM_DEFINITION
    _defaults = {**CCLibHarness._defaults, "name": "cclib-qchem", "node_parallel": False}


class ORCACCLibHarness(CCLibHarness):
    """Run ORCA and convert its output exclusively through cclib."""

    definition = ORCA_DEFINITION
    _defaults = {**CCLibHarness._defaults, "name": "cclib-orca", "node_parallel": True}
```

`parser_type()` performs its cclib parser import inside the function so package import remains lazy.

- [ ] **Step 5: Export and register concrete harnesses**

Create `cclib_programs/__init__.py`:

```python
"""Extensible cclib-backed native program harnesses."""

from .base import CCLibHarness
from .cclib_orca import ORCACCLibHarness
from .cclib_qchem import QChemCCLibHarness

__all__ = ["CCLibHarness", "ORCACCLibHarness", "QChemCCLibHarness"]
```

Update `qcengine/programs/base.py` to import the concrete classes and register:

```python
register_program(QChemCCLibHarness())
register_program(ORCACCLibHarness())
```

Delete `qcengine/programs/cclib.py` only after all imports use the package.

- [ ] **Step 6: Update internal test references and run package tests**

Map former monolith references exactly:

```text
shared execution/conversion helpers → cclib_base
Q-Chem builder/preflight/probe → `cclib_qchem.build_input/preflight/probe`
ORCA builder/preflight/probe → `cclib_orca.build_input/preflight/probe`
_PROGRAM_DEFINITIONS["qchem"] → cclib_qchem.QCHEM_DEFINITION
_PROGRAM_DEFINITIONS["orca"] → cclib_orca.ORCA_DEFINITION
_Job/_ExecutionResult → Job/ExecutionResult
```

Run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest qcengine/programs/cclib_programs/tests/test_cclib.py -q
/tmp/qcengine-cclib-venv/bin/python -c 'import sys, qcengine; assert not any(n == "cclib" or n.startswith("cclib.") for n in sys.modules)'
```

Expected: package tests pass and importing QCEngine leaves cclib unloaded.

- [ ] **Step 7: Commit the package split**

```bash
git add qcengine/programs/base.py qcengine/programs/cclib_programs qcengine/programs/cclib.py qcengine/programs/tests/test_cclib.py
git commit -m "refactor: split cclib harnesses by native program"
```

---

### Task 3: Pass methods through and align program-specific drivers

**Files:**
- Modify: `qcengine/programs/cclib_programs/base.py`
- Modify: `qcengine/programs/cclib_programs/cclib_qchem.py`
- Modify: `qcengine/programs/cclib_programs/cclib_orca.py`
- Modify: `qcengine/programs/cclib_programs/tests/test_cclib.py`

**Interfaces:**
- Produces: common `_input_fields(input_model) -> Tuple[str, str, str]` that validates basis and real atoms but not method support.
- Q-Chem accepts `energy`, `gradient`, `hessian`; ORCA accepts `energy` only.

- [ ] **Step 1: Replace allowlist tests with pass-through and driver tests**

Add these focused cases before changing production:

```python
@pytest.mark.parametrize("program_module,method", [
    (cclib_qchem, "wb97x-d"),
    (cclib_orca, "dlpno-ccsd(t)"),
])
def test_native_methods_pass_through_without_qcengine_allowlist(program_module, method):
    job = program_module.build_input(
        _atomic_input(method=method), _task_config(), f"/{program_module.__name__.split('_')[-1]}"
    )
    assert method in job.input_text


@pytest.mark.parametrize("driver,jobtype", [
    ("energy", "sp"), ("gradient", "force"), ("hessian", "freq")
])
def test_qchem_supports_verified_drivers(driver, jobtype):
    assert f"JOBTYPE {jobtype}" in cclib_qchem.build_input(
        _atomic_input(driver=driver), _task_config(), "/qchem"
    ).input_text


@pytest.mark.parametrize("driver", ["gradient", "hessian"])
def test_orca_rejects_unverified_cclib_drivers(driver):
    with pytest.raises(InputError, match="cclib-orca.*energy"):
        cclib_orca.build_input(_atomic_input(driver=driver), _task_config(), "/orca")
```

Use safe placeholder executable names such as `/qchem` and `/orca`, not machine-local paths.

- [ ] **Step 2: Run the new behavior tests and verify failure**

Run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest \
  qcengine/programs/cclib_programs/tests/test_cclib.py \
  -k 'pass_through or verified_drivers or unverified_cclib_drivers' -q
```

Expected: arbitrary methods fail the old allowlist and ORCA gradient/Hessian are still accepted.

- [ ] **Step 3: Implement common-only input validation and program driver policies**

In the base, replace `_validate_input_subset` and delete `_METHOD_ALIASES`:

```python
def input_fields(input_model: "AtomicInput") -> Tuple[str, str, str]:
    """Return native driver, method, and basis after shared structural validation."""
    driver_value = input_model.specification.driver
    driver = driver_value.value if hasattr(driver_value, "value") else str(driver_value)
    method = input_model.specification.model.method
    basis = input_model.specification.model.basis
    if not isinstance(basis, str) or not basis.strip():
        raise InputError("CCLibHarness basis must be a non-empty string")
    if not all(bool(real) for real in input_model.molecule.real):
        raise InputError("CCLibHarness requires all atoms to be real; ghost atoms are unsupported")
    return driver, method, basis
```

Q-Chem maps exactly:

```python
jobtypes = {"energy": "sp", "gradient": "force", "hessian": "freq"}
```

ORCA checks `driver.casefold() == "energy"` and otherwise raises an `InputError` stating that `cclib-orca` currently supports only energy because gradient/Hessian cclib parsing is unsupported.

During conversion, retain parsed driver and trimmed/casefolded basis checks. Remove requested-versus-parsed method rejection and preserve the writer's parsed method metadata in v1 output.

- [ ] **Step 4: Run focused and full package tests**

Run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest \
  qcengine/programs/cclib_programs/tests/test_cclib.py -q
```

Expected: all retained package tests pass.

- [ ] **Step 5: Commit behavior alignment**

```bash
git add qcengine/programs/cclib_programs
git commit -m "fix: align cclib drivers and native methods"
```

---

### Task 4: Consolidate tests and add canonical/live coverage

**Files:**
- Modify: `qcengine/programs/cclib_programs/tests/test_cclib.py`
- Modify: `qcengine/tests/test_harness_canonical.py:19-46`

**Interfaces:**
- Consumes: concrete selectors and availability markers.
- Produces: reduced package suite, canonical energy/gradient coverage, real version checks, and Q-Chem Hessian coverage.

- [ ] **Step 1: Reduce the package test suite by behavior class**

Consolidate the retained suite to approximately 20–30 test functions using parametrization. Keep one focused function for each of these contracts:

```text
registration and independent concrete instances
lazy cclib import
missing executable behavior
Q-Chem environment preflight
basis/ghost validation
Q-Chem exact input and driver mapping
Q-Chem scalar keyword rendering and malformed/reserved rejection
ORCA exact input and energy-only policy
ORCA deterministic blocks and malformed/reserved rejection
Q-Chem/ORCA primary-output selection
managed scratch and QCSCRATCH isolation
bounded execution/resource failure classification
parser temporary-file cleanup on success and failure
parser detection/identity/parse failure classification
writer required fields/extras collision/validation failure
successful v1-to-v2 conversion and provenance preservation
parsed driver/basis mismatch rejection
native-file protocols and public schema conversion
live installed versions
live Q-Chem Hessian
live Q-Chem and ORCA energy calculations
real cclib fixture conversion
```

Delete separate tests whose assertions are subsumed by a parameterized contract, especially duplicate malformed-key variants and duplicate failure-chain cases. Do not weaken security boundaries for native keyword injection, extras ownership, or temporary-file cleanup.

- [ ] **Step 2: Add canonical selector entries**

Add to `_canonical_methods`:

```python
pytest.param("cclib-qchem", {"method": "hf", "basis": "sto-3g"}, {}, marks=using("cclib-qchem")),
pytest.param("cclib-orca", {"method": "hf", "basis": "sto-3g"}, {}, marks=using("cclib-orca")),
```

The existing canonical energy test will exercise both. The existing canonical gradient test parametrizes the same list, so skip ORCA there with an explicit branch analogous to ADCC:

```python
if program == "cclib-orca":
    with pytest.raises(qcng.exceptions.InputError, match="only energy"):
        qcng.compute(inp, program, raise_error=True, return_version=retver)
    return
```

This makes Q-Chem gradient canonical while verifying ORCA's documented rejection. Do not add these selectors to structured-basis coverage because both require native basis strings.

- [ ] **Step 3: Add real version and Hessian tests**

Use installed software through availability gates:

```python
@pytest.mark.parametrize("selector,minimum", [
    pytest.param("cclib-qchem", "5.1", marks=using("cclib-qchem")),
    pytest.param("cclib-orca", "6.0", marks=using("cclib-orca")),
])
def test_live_version_uses_available_software(selector, minimum):
    assert parse_version(qcng.get_program(selector).get_version()) >= parse_version(minimum)


@uusing("cclib-qchem")
def test_live_qchem_hessian():
    result = qcng.compute(
        _atomic_input(driver="hessian", method="hf", basis="sto-3g"),
        "cclib-qchem",
        raise_error=True,
        task_config={"ncores": 1, "memory": 1.0},
    )
    assert result.success is True
    assert result.return_result.shape == (6, 6)
```

- [ ] **Step 4: Run non-live and canonical tests**

Run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest \
  qcengine/programs/cclib_programs/tests/test_cclib.py \
  qcengine/tests/test_harness_canonical.py \
  -m 'not cclib_qchem and not cclib_orca' -q
```

Expected: non-live tests pass; availability-gated live cases are skipped or deselected.

- [ ] **Step 5: Run installed-software tests**

After sourcing the existing machine-local Q-Chem setup and adding ORCA to `PATH` only in the shell environment, run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest \
  qcengine/programs/cclib_programs/tests/test_cclib.py \
  qcengine/tests/test_harness_canonical.py \
  -m 'cclib_qchem or cclib_orca' -q
```

Expected: real version checks, Q-Chem Hessian, canonical Q-Chem energy/gradient, and canonical ORCA energy/rejection behavior pass.

- [ ] **Step 6: Commit the consolidated coverage**

```bash
git add qcengine/programs/cclib_programs/tests/test_cclib.py qcengine/tests/test_harness_canonical.py
git commit -m "test: consolidate cclib harness coverage"
```

---

### Task 5: Update documentation and complete integration verification

**Files:**
- Modify: `docs/source/programs_cclib.rst`

**Interfaces:**
- Documents the final package extension contract and verified support matrix.

- [ ] **Step 1: Update the support matrix and method contract**

Change the driver table to:

```rst
+--------------+-----------------+---------------------+
| Driver       | Q-Chem job type | ORCA support        |
+==============+=================+=====================+
| ``energy``   | ``sp``          | supported           |
+--------------+-----------------+---------------------+
| ``gradient`` | ``force``       | not supported       |
+--------------+-----------------+---------------------+
| ``hessian``  | ``freq``        | not supported       |
+--------------+-----------------+---------------------+
```

Replace the fixed method list with text stating that method strings pass through to the native program and unsupported methods fail downstream. Explain that ORCA gradient/Hessian are disabled because current cclib parsing does not support the verified output path.

- [ ] **Step 2: Document the extension layout**

Add an “Adding another cclib-backed program” section containing the exact package tree and these responsibilities:

```text
base.py: shared lazy loading, execution, parsing, and QCSchema conversion
cclib_<program>.py: native input, preflight, version probe, output selection, parser type, concrete subclass
ProgramDefinition: immutable callback contract consumed by CCLibHarness
```

State that program modules must not import cclib at module import time and every helper must have a responsibility-focused docstring.

- [ ] **Step 3: Build documentation with warnings as errors**

Run:

```bash
cd docs
/tmp/qcengine-cclib-venv/bin/python -m sphinx -W -b html source build/html
```

Expected: build succeeds with no warnings.

- [ ] **Step 4: Run fixture verification**

With `CCLIB_SOURCE_ROOT` set only in the shell environment, run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest \
  qcengine/programs/cclib_programs/tests/test_cclib.py \
  -k 'fixture' -q
```

Expected: seven supported fixtures pass; the deferred ORCA DFT frequency fixture is removed because ORCA is now energy-only.

- [ ] **Step 5: Run the complete non-addon suite**

Run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m pytest qcengine \
  -m 'not addon' -q
```

Expected: all non-addon tests pass.

- [ ] **Step 6: Run static and repository hygiene checks**

Run:

```bash
/tmp/qcengine-cclib-venv/bin/python -m compileall -q qcengine \
git diff --check
! grep -R -nE '/home/|/projects/|qchem_vars\.sh|orca_6_' \
  qcengine/programs/cclib_programs qcengine/tests/test_harness_canonical.py docs/source/programs_cclib.rst
```

Expected: all commands exit zero.

- [ ] **Step 7: Commit documentation**

```bash
git add docs/source/programs_cclib.rst
git commit -m "docs: describe extensible cclib harness package"
```

- [ ] **Step 8: Request independent whole-branch review**

Review the diff from `a97b46957b8e3d192cefda2855df47f5c8cd891a` through `HEAD`, with emphasis on all ten review comments, optional dependency behavior, module boundaries, support claims, test reduction, and machine-local path leakage. Resolve all Critical and Important findings, rerun affected tests, and record residual Minor findings.
