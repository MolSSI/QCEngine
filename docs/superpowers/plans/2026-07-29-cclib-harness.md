# Configurable cclib Program Harness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one configurable, optional cclib-backed QCEngine harness registered independently as `cclib-qchem` and `cclib-orca`, with generated native inputs, managed execution, cclib parsing/QCSchema conversion, live demonstrations, and documentation.

**Architecture:** A frozen `CCLibHarness` owns the shared lifecycle and selects behavior from a private immutable Q-Chem/ORCA definition table. External cclib imports remain lazy; executable identity/version probes are cached by resolved path; tests isolate proprietary execution behind monkeypatched discovery/execution and optional fixture/live gates. cclib writer output is validated as QCSchema v1, converted to QCSchema v2 with the original request as `input_data`, and returned to QCEngine's normal metadata/version-conversion layer.

**Tech Stack:** Python 3.10+, Pydantic 2, QCEngine `ProgramHarness`, QCElemental QCSchema v1/v2 models and discovery utilities, cclib development `QCSchemaWriter`, pytest, Sphinx/RST.

## Global Constraints

- Keep the native `qchem` selector registered and `qcengine/programs/qchem.py` unchanged.
- Register `CCLibHarness(name="cclib-qchem", program="qchem")` and `CCLibHarness(name="cclib-orca", program="orca")` as two frozen instances of the same class.
- Common harness settings are `scratch=True`, `thread_safe=False`, `thread_parallel=True`, and `managed_memory=True`; derive `node_parallel=False` for Q-Chem and `node_parallel=True` for ORCA from `program`.
- Treat ORCA `node_parallel=True` as scheduler capability metadata only; this implementation emits `%pal nprocs TaskConfig.ncores` and does not add a node launcher.
- Keep cclib optional. Importing QCEngine, listing registered programs, and using unrelated harnesses must not import cclib.
- Discover `qchem` and `orca` only through `PATH`; reusable code must not contain absolute executable/source paths or source shell scripts.
- Support only atomic `energy`, `gradient`, and `hessian` drivers; case-insensitive `hf`, `b3lyp`, `bp86`, `mp2`, and `ccsd`; non-empty string basis names; and real atoms.
- Generate native input from QCSchema; do not accept raw native input, ghost atoms, structured `BasisSet`, multi-job/restart files, or unsupported procedures/properties.
- Q-Chem derived/reserved `$rem` keys are authoritative and collisions raise `InputError`. ORCA accepts only the documented `simple` and `blocks` keyword structure; `pal` and `maxcore` are reserved.
- Use QCEngine `execute()` and managed scratch utilities. Parse complete output exclusively through cclib auto-detection and verify parser identity (`QChem` or `ORCA`).
- Call `QCSchemaWriter(ccdata).as_dict(validate=False)`, validate a QCSchema-v1 `AtomicResult`, then convert to internal QCSchema v2 with `external_input_data=input_model`.
- Preserve the original requested input, parsed output molecule, flat cclib extras, truthful QC-program provenance, complete primary output, and non-colliding `extras["cclib_harness"]` metadata.
- Never add ORCA dispersion twice or relabel driver/method/basis to conceal a parsed/requested mismatch.
- Raise `InputError`, `ResourceError`, and `UnknownError` according to the design. Every post-execution error includes selector, resolved executable, stage, and a tail bounded to 40 lines and 4000 characters, with exception chaining.
- `found(False)` returns `False` for every ordinary optional-resource/probe failure and never breaks global availability listing; `found(True)` raises a specific `ResourceError`.
- Optional fixture tests use only `CCLIB_SOURCE_ROOT`. Machine paths may appear only in clearly labeled local-verification documentation and shell commands, never reusable code or unit-test assumptions.
- Do not modify cclib. If an approved fixture exposes a cclib parser/writer defect, stop and report the fixture, parser, exception, and incorrect/missing fields before changing scope.
- Follow strict TDD: each production behavior first has a focused test that is observed failing for the expected missing behavior, then minimal production code, then green verification.

---

### Task 1: Register the configurable harness and implement optional-resource availability

**Files:**
- Create: `qcengine/programs/cclib.py`
- Modify: `qcengine/programs/base.py`
- Create: `qcengine/programs/tests/test_cclib.py`

**Interfaces:**
- Consumes: `ProgramHarness`, `TaskConfig`, `qcelemental.util.which`, `safe_version`, `parse_version`, QCEngine `execute`, `InputError`, `ResourceError`, and `UnknownError`.
- Produces:
  ```python
  class CCLibHarness(ProgramHarness):
      program: Literal["qchem", "orca"]
      version_cache: ClassVar[Dict[str, str]]
      def found(self, raise_error: bool = False) -> bool: ...
      def get_version(self) -> str: ...
      def compute(self, input_data: AtomicInput, config: TaskConfig) -> AtomicResult: ...
  ```
- Produces immutable `_ProgramDefinition` entries with selector, executable, minimum version, input/output names, expected parser, normal-termination marker, generator, probe, output selector, and preflight callable.
- Produces lazy `_load_cclib_api()` and cached `_check_cclib_compatibility()` helpers used by Task 4.

- [ ] **Step 1: Write registration/lazy-import tests**

  Add tests that case-insensitive lookup returns two distinct frozen `CCLibHarness` instances, native `qchem` remains separate, Q-Chem has `node_parallel=False`, ORCA has `node_parallel=True`, and a subprocess can import/list QCEngine while an import guard rejects `cclib` and `cclib.*`.

  ```python
  def test_registered_instances_are_frozen_and_independent():
      qchem = qcng.get_program("CCLIB-QCHEM", check=False)
      orca = qcng.get_program("cclib-orca", check=False)
      assert type(qchem) is type(orca) is CCLibHarness
      assert (qchem.program, qchem.node_parallel) == ("qchem", False)
      assert (orca.program, orca.node_parallel) == ("orca", True)
      assert qchem is not orca
      assert {"qchem", "cclib-qchem", "cclib-orca"} <= qcng.list_all_programs()
  ```

- [ ] **Step 2: Run the registration tests and verify RED**

  Run:
  ```bash
  python -m pytest \
    qcengine/programs/tests/test_cclib.py::test_registered_instances_are_frozen_and_independent \
    qcengine/programs/tests/test_cclib.py::test_import_qcengine_does_not_import_external_cclib -q
  ```
  Expected: collection/import failure because `qcengine.programs.cclib` and selectors do not exist.

- [ ] **Step 3: Add the minimal frozen class and registrations**

  Keep cclib imports out of module scope. Declare `program` because `ProgramHarness` forbids undeclared extras, derive `node_parallel` during construction, and register exactly:

  ```python
  register_program(CCLibHarness(name="cclib-qchem", program="qchem"))
  register_program(CCLibHarness(name="cclib-orca", program="orca"))
  ```

- [ ] **Step 4: Run registration tests and verify GREEN**

  Run the Step 2 command plus `python -m pytest qcengine/tests/test_program_utils.py -q`.

- [ ] **Step 5: Write failing cclib compatibility tests**

  Build a synthetic one-atom neutral singlet HF `ccData` with a 1 Å coordinate, SCF energy/history, basis metadata, and Mulliken charge. Test missing cclib, valid QCSchema-v1 writer output, Angstrom→bohr geometry near `1.8897261255`, preserved flat `atomcoords`/`atomcharges`, wrong-unit rejection, missing-flat-extra rejection, and one-process cache behavior.

- [ ] **Step 6: Verify compatibility RED**

  Run:
  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "compatibility or cclib_missing"
  ```
  Expected: helper/probe is absent.

- [ ] **Step 7: Implement lazy cclib loading and feature compatibility probe**

  Import cclib, `ccData`, `QCSchemaWriter`, and parser APIs inside helpers only. Validate writer output as explicit QCElemental v1 `AtomicResult` (use the repository's `_v1v2.AtomicResult` fallback only where Python 3.14 requires it), geometry conversion, and flat extras. Cache `(success, version_or_message)` and translate failures to `ResourceError` only at the availability boundary.

- [ ] **Step 8: Write failing executable identity/version/cache tests**

  Monkeypatch `which` and module-local `execute`. Test missing executables; Q-Chem banner plus version ≥5.1; ORCA banner plus `Program Version` ≥6.0 and normal termination; rejection of an unrelated executable named `orca`; `get_version()` returning the external-program version; and one probe per resolved path using fake `/opt/...` paths.

- [ ] **Step 9: Verify executable-probe RED**

  Run:
  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "executable or identity or version or probe_cache or unrelated_orca"
  ```

- [ ] **Step 10: Implement PATH discovery and cached probes**

  Use resolved paths as keys in public `version_cache`. Q-Chem uses a lightweight `$rem` input and requires both its banner/identity phrase and parseable version. ORCA uses a cached minimal H-doublet HF/STO-3G input when needed and requires banner, version, and normal termination. Compare with `parse_version`; normalize with `safe_version`; never cache a failed identity/version.

- [ ] **Step 11: Write failing Q-Chem environment and `found()` orchestration tests**

  Parametrize missing/invalid `QC`, `QCAUX`, and `QCPROG`, use temporary directories/files, assert one error names every invalid variable, require readable directories/executable driver/resolved executable, and confirm inherited `QCSCRATCH` is not required. Force each compatibility/discovery/probe failure and assert `found(False)` and `list_available_programs()` do not raise while `found(True)` raises a specific `ResourceError`.

- [ ] **Step 12: Verify environment/orchestration RED**

  Run:
  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "environment or preflight or found or available_programs"
  ```

- [ ] **Step 13: Implement environment preflight and complete `found()`**

  Build child environments from `os.environ.copy()`. Check compatibility → PATH → Q-Chem environment (when selected) → identity/version in order. Aggregate invalid environment resources. `found(False)` catches ordinary probe failures and returns `False`; `found(True)` preserves detailed `ResourceError` messages.

- [ ] **Step 14: Run Task 1 tests**

  ```bash
  python -m pytest \
    qcengine/programs/tests/test_cclib.py \
    qcengine/tests/test_program_utils.py -q \
    -k "registered or import_qcengine or compatibility or missing or identity or version or probe or environment or found or available"
  ```

- [ ] **Step 15: Commit Task 1**

  ```bash
  git add qcengine/programs/cclib.py qcengine/programs/base.py qcengine/programs/tests/test_cclib.py
  git commit -m "feat: register and probe cclib harnesses"
  ```

---

### Task 2: Validate the supported input subset and generate deterministic native inputs

**Files:**
- Modify: `qcengine/programs/cclib.py`
- Modify: `qcengine/programs/tests/test_cclib.py`

**Interfaces:**
- Consumes: Task 1 `_ProgramDefinition` and resolved executable.
- Produces:
  ```python
  @dataclass(frozen=True)
  class _Job:
      command: List[str]
      infiles: Dict[str, str]
      outfiles: List[str]
      input_filename: str
      output_filename: str
      input_text: str
      executable: str

  def _validate_input_subset(input_model: AtomicInput) -> Tuple[str, str, str]: ...
  def _build_qchem_input(input_model: AtomicInput, config: TaskConfig, executable: str) -> _Job: ...
  def _build_orca_input(input_model: AtomicInput, config: TaskConfig, executable: str) -> _Job: ...
  ```
- `CCLibHarness.build_input()` returns the existing harness-style dictionary containing `commands`, `infiles`, `outfiles`, and `scratch_directory` while delegating to `_Job` internally.

- [ ] **Step 1: Write shared input-subset failure tests**

  Use local v2 models, not QCEngineRecords. Parametrize accepted drivers/methods and reject unsupported driver/method, absent/empty/structured basis, and any `molecule.real is False` before execution. Include upper/mixed-case accepted methods and preserve requested method/basis spelling in generated text.

- [ ] **Step 2: Verify shared validation RED**

  Run:
  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "input_subset"
  ```

- [ ] **Step 3: Implement `_validate_input_subset()`**

  Normalize only comparisons. Return `(driver, method, basis)` strings; do not silently substitute aliases or alter the request.

- [ ] **Step 4: Write failing Q-Chem generator tests**

  Cover all driver mappings (`sp`, `force`, `freq`), all five methods, closed/open shell charge/multiplicity, the exact MP2 water bohr geometry, `MEM_TOTAL=int(memory*1024)`, `-nt ncores`, four harvesting defaults, deterministic ordinary keyword ordering, and scalar rendering (`str`, bool, int, finite float). Reject newline/`None`/container/non-finite values and case-insensitive collisions with:

  ```python
  QCHEM_RESERVED = {
      "JOBTYPE", "METHOD", "BASIS", "MEM_TOTAL", "INPUT_BOHR",
      "SCF_FINAL_PRINT", "PRINT_GENERAL_BASIS", "PRINT_ORBITALS", "MOLDEN_FORMAT",
  }
  ```

- [ ] **Step 5: Verify Q-Chem generator RED**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "qchem_input"
  ```

- [ ] **Step 6: Implement deterministic Q-Chem generation**

  Emit `$comment`, `$molecule`, and `$rem`; preserve atom order; emit geometry in bohr plus `INPUT_BOHR=TRUE`; make derived/default keys authoritative; render commands as:

  ```python
  [executable, "-nt", str(config.ncores), "dispatch.in", "dispatch.out"]
  ```

- [ ] **Step 7: Run Q-Chem generator tests and verify GREEN**

  Re-run Step 5.

- [ ] **Step 8: Write failing ORCA generator tests**

  Cover driver keywords (none/`engrad`/`freq`), all methods, charge/multiplicity, atom order, bohr→Å conversion, `%pal`, and `%MaxCore=max(1, int(memory*1024/ncores))`. Assert exact default `%output` lines; preserve `simple` list order; sort block names; append user output-body lines after defaults. Reject unknown top-level keys, malformed/non-list/newline simple values, malformed blocks, invalid block names, non-string bodies, and case-insensitive reserved `pal`/`maxcore` or coordinate injection.

  For `ncores=4`, `memory=2.734375`, assert:
  ```text
  %pal
  nprocs 4
  end
  %MaxCore 700
  ```

- [ ] **Step 9: Verify ORCA generator RED**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "orca_input"
  ```

- [ ] **Step 10: Implement deterministic ORCA generation**

  Emit `! method basis [driver] [simple...]`, `%output` and sorted user blocks, resource blocks, and `* xyz charge multiplicity`; command is `[executable, "dispatch.inp"]`.

- [ ] **Step 11: Run Task 2 tests**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "input_subset or qchem_input or orca_input"
  ```

- [ ] **Step 12: Commit Task 2**

  ```bash
  git add qcengine/programs/cclib.py qcengine/programs/tests/test_cclib.py
  git commit -m "feat: generate cclib program inputs"
  ```

---

### Task 3: Execute jobs in managed scratch and classify execution failures

**Files:**
- Modify: `qcengine/programs/cclib.py`
- Modify: `qcengine/programs/tests/test_cclib.py`

**Interfaces:**
- Consumes: Task 2 `_Job`, QCEngine `execute()` and `temporary_directory()`.
- Produces:
  ```python
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

  def _diagnostic_tail(text: str, max_lines: int = 40, max_chars: int = 4000) -> str: ...
  def _execute_job(definition: _ProgramDefinition, job: _Job, config: TaskConfig) -> _ExecutionResult: ...
  ```

- [ ] **Step 1: Write failing output-selection and scratch tests**

  Monkeypatch module-local `execute()`. Assert Q-Chem requests/selects `dispatch.out`; ORCA selects captured stdout; complete inherited environment is passed; Q-Chem child `QCSCRATCH` points to QCEngine-managed temporary scratch; and commands/infiles from Task 2 are unchanged.

- [ ] **Step 2: Verify execution RED**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "execution_selection or managed_scratch"
  ```

- [ ] **Step 3: Implement managed execution and output selection**

  Q-Chem runs with a managed `QCSCRATCH` and requested outfile. ORCA runs through the same QCEngine utility and consumes stdout. Do not add raw `subprocess` or machine paths.

- [ ] **Step 4: Write failing termination/error tests**

  Cover nonzero exit, missing primary output, zero exit without normal marker, undefined Q-Chem environment variable extraction, license strings (`FlexNet`, `license checkout`, `unable to validate license`), and execution exceptions. Assert environment/license cases are `ResourceError`; other failures are `UnknownError`; every message names selector, executable, stage, and bounded tail; original exceptions remain chained.

- [ ] **Step 5: Verify error-classification RED**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "termination or execution_error or license or diagnostic_tail"
  ```

- [ ] **Step 6: Implement termination checks and classifiers**

  Require successful process plus configured normal marker. Extract actual missing variable names when possible. Implement the exact 40-line/4000-character tail bound and use one stage-aware error formatter.

- [ ] **Step 7: Run Task 3 tests**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "execution or scratch or termination or license or diagnostic"
  ```

- [ ] **Step 8: Commit Task 3**

  ```bash
  git add qcengine/programs/cclib.py qcengine/programs/tests/test_cclib.py
  git commit -m "feat: execute cclib program jobs"
  ```

---

### Task 4: Parse with cclib and convert validated results to QCSchema v2

**Files:**
- Modify: `qcengine/programs/cclib.py`
- Modify: `qcengine/programs/tests/test_cclib.py`

**Interfaces:**
- Consumes: Task 1 lazy cclib API/definition, Task 3 `_ExecutionResult`, original v2 `AtomicInput`.
- Produces:
  ```python
  def _parse_and_convert(
      definition: _ProgramDefinition,
      execution: _ExecutionResult,
      input_model: AtomicInput,
  ) -> AtomicResult: ...
  ```
- Completes `CCLibHarness.compute()` as validate → build → execute → parse/convert.

- [ ] **Step 1: Write failing parser/writer failure tests**

  Use fake lazy cclib APIs (no installed cclib required). Cover auto-detection returning `None`, parser class mismatch, parser exception, incomplete result, missing/false `metadata["success"]`, writer exception/unsupported method, incomplete writer fields, and QCElemental v1 validation failure. Assert stage-aware `UnknownError`, bounded tail, and chaining.

- [ ] **Step 2: Verify parser failure RED**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "parser or writer_failure or validation_failure"
  ```

- [ ] **Step 3: Write failing successful-conversion tests**

  Fake Q-Chem and ORCA parser classes and writer dictionaries. Assert parser identity, v1 validation, internal v2 type, original input equality, parsed output molecule retention, complete output at top-level `stdout`, cclib flat extras unchanged, truthful QC-program provenance, and exactly one added object:

  ```python
  extras["cclib_harness"] = {
      "selector": definition.selector,
      "cclib_version": cclib_version,
      "parser": definition.parser_name,
      "executable": execution.executable,
  }
  ```

  Test driver/method/basis case-insensitive matching plus only these method aliases:

  ```python
  {"rhf": "hf", "uhf": "hf", "rmp2": "mp2", "ump2": "mp2",
   "rccsd": "ccsd", "uccsd": "ccsd"}
  ```

  Material mismatch must raise `UnknownError` rather than relabel data.

- [ ] **Step 4: Write failing native-file protocol tests**

  Assert protocol outcomes:
  ```python
  [("none", set()), ("input", {"input"}), ("all", {"input", "dispatch.out"})]
  ```
  Never emit keys named `stdout` or `stderr`, never add `extras["outfiles"]`, and keep ORCA captured output in top-level `stdout` (plus filename-like `dispatch.out` only for protocol `all`). If active v1 models reject `native_files`, retain generated input as `extras["cclib_harness"]["native_input"]` only for that specific validation incompatibility.

- [ ] **Step 5: Verify successful-conversion RED**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "conversion or mismatch or native_files or protocols"
  ```

- [ ] **Step 6: Implement cclib auto-detection and conversion**

  Feed complete text through a seekable named stream, call cclib auto-detection (`ccopen`) to retain parser identity, call `parse()`, and close the parser input in `finally`. Call `QCSchemaWriter(...).as_dict(validate=False)`, attach approved data without replacing extras, validate explicit v1, then:

  ```python
  result_v2 = result_v1.convert_v(2, external_input_data=input_model)
  ```

- [ ] **Step 7: Run mocked conversion tests**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "parser or writer or conversion or mismatch or native_files or protocols"
  ```

- [ ] **Step 8: Write and run public dispatch conversion test**

  Monkeypatch resources/execution/cclib, call `qcng.compute(..., "cclib-qchem", return_version=1, return_dict=False)`, and assert public schema v1 while direct harness output remains v2 and preserves the request.

- [ ] **Step 9: Add optional real-cclib fixture integration tests**

  Gate only on `CCLIB_SOURCE_ROOT`; skip if `<root>/data` or a named file is absent. Cover at least:
  - Q-Chem: `water_mp2.out`, a B3LYP/BP86 energy fixture, `water_ir.out`, `water_ccsd.out`.
  - ORCA: `water_mp2.out`, an HF/B3LYP energy fixture, `dvb_ir.out`, `water_ccsd.out`.

  Build requested molecule/identity from an initial real cclib writer pass, then send the same complete output through `_parse_and_convert()` and assert success, parser/provenance identity, driver/method/basis, non-empty flat extras, and original-input preservation.

- [ ] **Step 10: Run fixture integration**

  ```bash
  CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib \
    python -m pytest qcengine/programs/tests/test_cclib.py -q -k fixture
  ```

  If a failure is in cclib parsing/writer behavior rather than QCEngine integration, stop and report it without modifying cclib.

- [ ] **Step 11: Run Task 4 tests**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k "parser or writer or conversion or mismatch or native_files or protocols or fixture"
  ```

- [ ] **Step 12: Commit Task 4**

  ```bash
  git add qcengine/programs/cclib.py qcengine/programs/tests/test_cclib.py
  git commit -m "feat: convert cclib output to qcschema"
  ```

---

### Task 5: Add live-test integration and Q-Chem/ORCA demonstrations

**Files:**
- Modify: `qcengine/testing.py`
- Modify: `pyproject.toml`
- Modify: `qcengine/programs/tests/test_cclib.py`
- Create: `qchem_water_mp2.py`
- Create: `orca_water_ccsd.py`

**Interfaces:**
- Consumes: completed selectors and public `qcengine.compute()`.
- Produces pytest markers/availability keys and two executable demonstration scripts with `compare_result(result) -> List[str]` and `main(output_path=...) -> int`.

- [ ] **Step 1: Write failing testing-registration tests**

  Assert `has_program()` recognizes `cclib-qchem` and `cclib-orca`, and `uusing()` composes addon plus exact selector marks.

- [ ] **Step 2: Verify marker RED**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k testing_registration
  ```

- [ ] **Step 3: Add availability entries and markers**

  In `testing.py` add external-program version floors:
  ```python
  "cclib-qchem": is_program_new_enough("cclib-qchem", "5.1"),
  "cclib-orca": is_program_new_enough("cclib-orca", "6.0"),
  ```
  In `pyproject.toml` register explicit `cclib-qchem` and `cclib-orca` addon markers. Do not add cclib to core dependencies.

- [ ] **Step 4: Add separately marked live smoke tests**

  Add water MP2/STO-3G energy tests for both selectors and a water CCSD/STO-3G ORCA test. Mark only real execution tests with `@uusing(...)`; all unit/fixture tests remain offline.

- [ ] **Step 5: Write failing demonstration-unit tests**

  Import each script and test comparison logic with in-memory result dictionaries. Monkeypatch `qcengine.compute` in `main()` tests and assert selectors, `raise_error=True`, `return_version=1`, complete JSON output, concise PASS/FAIL reporting, and nonzero failure status. Scripts must not read historical `.out` or JSON files.

  Q-Chem comparison covers every structural/numerical/rich-extras criterion in the design, including seven historical MO energies and six SCF rows.

  ORCA uses the CCSD fixture geometry converted/stored in bohr, `task_config={"ncores": 4, "memory": 2.734375}`, and checks:
  - return/CCSD total `-75.013487814 ± 1e-6 Eh`;
  - SCF total `-74.96357424008319 ± 1e-6 Eh`;
  - CCSD correlation `-0.04991357391681104 ± 1e-6 Eh`;
  - `nbasis=nmo=7`, `nalpha=nbeta=5`, `natom=3`;
  - ORCA provenance, `cclib-orca` metadata, and representative flat CCSD/orbital/SCF extras.

  Preserve but do not numerically endorse the known anomalous MP2 properties in the historical ORCA CCSD writer output.

- [ ] **Step 6: Verify demonstration RED**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -k demonstration
  ```

- [ ] **Step 7: Implement `qchem_water_mp2.py` and `orca_water_ccsd.py`**

  Construct QCSchema inputs directly, call the exact selectors, print one line per comparison, write `qchem_water_mp2.result.json` or `orca_water_ccsd.result.json`, and exit nonzero if any comparison fails.

- [ ] **Step 8: Run offline Task 5 tests**

  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q \
    -m "not cclib-qchem and not cclib-orca" \
    -k "testing_registration or demonstration"
  ```

- [ ] **Step 9: Run configured live tests and demonstrations**

  After preparing PATH/Q-Chem environment exactly as documented:
  ```bash
  python -m pytest qcengine/programs/tests/test_cclib.py -q -m "cclib-qchem"
  python -m pytest qcengine/programs/tests/test_cclib.py -q -m "cclib-orca"
  python qchem_water_mp2.py
  python orca_water_ccsd.py
  ```

  If proprietary resources are unavailable, record skips/resource errors; do not claim live acceptance.

- [ ] **Step 10: Commit Task 5**

  ```bash
  git add pyproject.toml qcengine/testing.py qcengine/programs/tests/test_cclib.py qchem_water_mp2.py orca_water_ccsd.py
  git commit -m "test: add cclib live demonstrations"
  ```

---

### Task 6: Document the selectors and run the final validation gates

**Files:**
- Create: `docs/source/programs_cclib.rst`
- Modify: `docs/source/program_overview.rst`
- Modify: `docs/source/index.rst`
- Modify: `docs/source/changelog.rst`

**Interfaces:**
- Consumes: final behavior from Tasks 1–5.
- Produces: linked user documentation and unreleased changelog entry.

- [ ] **Step 1: Build docs before adding the page and record RED**

  Run:
  ```bash
  python -m sphinx -W -b html docs/source /tmp/qcengine-docs
  ```
  The behavioral RED condition is that no cclib selector documentation/page exists yet; record the current build status without changing code to force a failure.

- [ ] **Step 2: Write the dedicated cclib documentation**

  Document both selectors and continued native `qchem`; supported drivers/methods/string basis and molecule limits; Q-Chem flat keywords and ORCA `simple`/`blocks`; reserved resource keys; optional cclib development install; PATH-only discovery; environment initialization before Python; v1 writer validation/internal v2 conversion; flat extras and harness metadata; native-file protocols; fixture/live commands; and machine-local verification exports clearly labeled non-portable.

- [ ] **Step 3: Update overview, toctree, and changelog**

  Add separate non-production E/G/H rows for `cclib-qchem` and `cclib-orca` while retaining native Q-Chem. Link `programs_cclib` from the docs toctree. Under v0.51.0 New Features add:

  ```rst
  - cclib - add optional cclib-backed Q-Chem and ORCA program selectors.
  ```

  Do not invent a PR number or edit redirect-only `CHANGELOG.md`.

- [ ] **Step 4: Build docs and verify GREEN**

  ```bash
  python -m sphinx -W -b html docs/source /tmp/qcengine-docs
  ```

- [ ] **Step 5: Run the full offline and fixture suite**

  ```bash
  python -m pytest \
    qcengine/programs/tests/test_cclib.py \
    qcengine/tests/test_program_utils.py -q \
    -m "not cclib-qchem and not cclib-orca"

  CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib \
    python -m pytest qcengine/programs/tests/test_cclib.py -q -k fixture

  python -m pytest \
    qcengine/programs/tests/test_canonical_fields.py \
    qcengine/programs/tests/test_canonical_config.py -q

  python -m qcengine info
  ```

- [ ] **Step 6: Audit scope and portability**

  Verify from the final diff:
  - no module-scope external cclib import;
  - no reusable source/unit test contains `/projects/`, `/home/awallace43/`, or setup-script names;
  - all executable caches use resolved paths;
  - `found(False)` cannot leak ordinary probe failures;
  - child environments start from `os.environ.copy()`;
  - no native-file key is `stdout` or `stderr`;
  - parser/writer/validation errors are chained and bounded;
  - mismatches are rejected rather than relabeled;
  - `qcengine/programs/qchem.py` and the cclib checkout are unchanged.

- [ ] **Step 7: Commit Task 6**

  ```bash
  git add docs/source/programs_cclib.rst docs/source/program_overview.rst docs/source/index.rst docs/source/changelog.rst
  git commit -m "docs: describe cclib program harnesses"
  ```

- [ ] **Step 8: Request whole-branch review**

  Generate a review package from the branch merge-base through `HEAD`; give the final reviewer this plan, the approved design, the task ledger, and the complete diff. Fix one synthesized wave of Critical/Important findings, re-review that fix range once, and disposition residual findings according to the subagent-driven-development breaker rules.
