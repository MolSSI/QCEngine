# Feature: Gaussian Harness — Runner <!-- rq-436a54d5 -->

This document specifies the top-level harness class (`GaussianHarness`) that integrates
Gaussian 16 (or Gaussian 09) into QCEngine. It follows the multi-file subdirectory
pattern established by the GAMESS and NWChem harnesses. The harness lives under
`qcengine/programs/gaussian/` and is registered alongside other programs in
`qcengine/programs/base.py`.

Related requirements files:
- `input-builder.md` — `.com` file construction (germinate, keywords)
- `harvester.md`     — cclib-based output parsing

## Feature API <!-- rq-cf6fdf8a -->

### Class `GaussianHarness(ProgramHarness)` <!-- rq-998da273 -->

Located at `qcengine/programs/gaussian/runner.py`.

**`_defaults`**:

| Key              | Value       |
|------------------|-------------|
| `name`           | `"Gaussian"` |
| `scratch`        | `True`      |
| `thread_safe`    | `False`     |
| `thread_parallel`| `True`      |
| `node_parallel`  | `False`     |
| `managed_memory` | `True`      |

#### `found(raise_error: bool = False) -> bool`

- Returns `True` if and only if **both** of the following hold:
  1. `g16` or `g09` is present on `PATH` (checked in that priority order via `which`).
  2. `cclib` is importable (checked via `which_import`).
- If either condition fails and `raise_error=False`, returns `False`.
- If either condition fails and `raise_error=True`, raises `ModuleNotFoundError` with a
  descriptive install hint.
- Does not execute any binary.

#### `get_version() -> str`

- Calls `found(raise_error=True)`.
- Runs the located Gaussian executable with a blank `.com` input file; captures stdout.
- Searches for a line matching `Gaussian <NN>, Revision <TAG>,` (e.g. `Gaussian 16,
  Revision C.02,`) and returns a safe version string of the form `"<NN>.<TAG>"` (e.g.
  `"16.C.02"`), produced via `safe_version`.
- Caches the result in `version_cache: Dict[str, str]` keyed by the absolute executable
  path so the binary is not re-executed within the same process.
- Raises `UnknownError` if no matching version line is found in the output.

#### `compute(input_model: AtomicInput, config: TaskConfig) -> AtomicResult`

- Calls `found(raise_error=True)`.
- Raises `InputError` immediately if `input_model.specification.model.basis` is a `BasisSet` instance
  (only string basis names are supported).
- Calls `self.build_input(input_model, config)` to produce the job record.
- Calls `self.execute(job_inputs)` to run Gaussian.
- Inspects the collected log for known fatal error strings (see Error Dispatch below).
- On success, calls `self.parse_output(outfiles, input_model)` and returns the result.

**Error dispatch** (applied to log content before calling `parse_output`):

| Pattern in log | Exception raised |
|----------------|-----------------|
| `"galloc: could not allocate memory"` | `ResourceError` |
| `"Ernie: is not defined as a Gaussian input file"` or `"Ernie: unrecognized"` | `InputError` |
| `"Error in internal coordinate system"` | `InputError` |
| `"Normal termination"` absent and no matched pattern above | `UnknownError` |

#### `build_input(input_model: AtomicInput, config: TaskConfig, template: Optional[str] = None) -> Dict[str, Any]`

- Raises `InputError` if `input_model.specification.model.basis` is a `BasisSet` instance or `None`.
- Delegates method/driver/multiplicity mapping to `muster_modelchem()` (germinate.py).
- Delegates `.com` file assembly to `build_com_file()` (keywords.py).
- Returns a dict:
  ```python
  {
      "infiles": {"gaussian.com": <formatted .com string>},
      "command": [<absolute path to g16 or g09>],
      "scratch_directory": config.scratch_directory,
      "scratch_messy": config.scratch_messy,
  }
  ```

#### `execute(inputs: Dict[str, Any], ...) -> Tuple[bool, Dict]`

- Calls `qcengine.util.execute()` with:
  - `inputs["command"]` as the command.
  - `inputs["infiles"]` (containing `"gaussian.com"`) as input files.
  - `["gaussian.log"]` as the list of output files to collect.
  - `scratch_directory` and `scratch_messy` forwarded from `inputs`.
- Returns `(success, dexe)` where `dexe["outfiles"]["gaussian.log"]` holds the full log
  text and `dexe["stdout"]` / `dexe["stderr"]` hold the process streams.

#### `parse_output(outfiles: Dict[str, str], input_model: AtomicInput) -> AtomicResult`

- Pops `"stdout"`, `"stderr"`, and `"input"` from `outfiles` before passing the rest to
  `harvest()`.
- Calls `harvest(input_model.molecule, method, log_text)` from `harvester.py` to obtain
  `(qcvars, grad, hess, mol)`.
- Raises `UnknownError` (wrapping the raw exception) if `harvest()` raises.
- When `grad is not None`, sets both `"CURRENT GRADIENT"` and
  `"<METHOD_UPPER> TOTAL GRADIENT"` in `qcvars`.
- When `hess is not None`, sets both `"CURRENT HESSIAN"` and
  `"<METHOD_UPPER> TOTAL HESSIAN"` in `qcvars`.
- Resolves `return_result` from `qcvars`:
  - Driver `"properties"` → `qcvars["CURRENT ENERGY"]`
  - All other drivers → `qcvars[f"CURRENT {driver.upper()}"]`
- Raises `UnknownError` if the required `qcvars` key is absent.
- Calls `build_out(qcvars)` and `build_atomicproperties(qcvars)`.
- Returns:
  ```python
  AtomicResult(**{
      "schema_version": 2,
      "input_data": input_model,
      "molecule": mol,
      "extras": {**input_model.specification.extras, "qcvars": qcvars_serialised},
      "native_files": {"gaussian.com": input_text, "gaussian.log": log_text},
      "properties": atprop,
      "provenance": Provenance(creator="Gaussian", version=self.get_version(),
                               routine="gaussian").model_dump(),
      "return_result": retres,
      "stderr": stderr,
      "stdout": stdout,
      "success": True,
  })
  ```

### Registration <!-- rq-17196f6b -->

`GaussianHarness` must be imported and registered in `qcengine/programs/base.py`:

```python
from .gaussian import GaussianHarness
# ...
register_program(GaussianHarness())
```

It must also be exported from `qcengine/programs/gaussian/__init__.py`.

## Gaussian Executable Details <!-- rq-cb0f61f2 -->

- Preferred executable: `g16`; fall back to `g09` if `g16` is not found.
- Gaussian reads its input from a `.com` file passed as its sole command-line argument.
- Gaussian writes its log to a file whose name derives from the `.com` file
  (e.g. `gaussian.com` → `gaussian.log` in the scratch directory).
- A non-zero exit code does **not** reliably signal failure; the log content is the
  authoritative source of success/failure status.
- `"Normal termination of Gaussian"` at the end of the log indicates success.

## Gherkin Scenarios <!-- rq-2bd5902e -->

```gherkin
Feature: GaussianHarness runner

  # --- found() ---

  @rq-4cf7c460
  Scenario: found() returns True when g16 and cclib are both available
    Given "g16" is present on the PATH
    And "cclib" is importable
    When GaussianHarness.found() is called
    Then it returns True

  @rq-d3b8d8c3
  Scenario: found() falls back to g09 when g16 is absent
    Given "g16" is not on the PATH
    And "g09" is present on the PATH
    And "cclib" is importable
    When GaussianHarness.found() is called
    Then it returns True

  @rq-8e21b7d7
  Scenario: found() returns False when neither g16 nor g09 is on PATH
    Given "g16" is not on the PATH
    And "g09" is not on the PATH
    When GaussianHarness.found() is called
    Then it returns False

  @rq-b5e72c6d
  Scenario: found() returns False when cclib is not importable
    Given "g16" is present on the PATH
    But "cclib" is not importable
    When GaussianHarness.found() is called
    Then it returns False

  @rq-bd98d417
  Scenario: found(raise_error=True) raises ModuleNotFoundError when g16 is missing
    Given "g16" is not on the PATH
    And "g09" is not on the PATH
    When GaussianHarness.found(raise_error=True) is called
    Then ModuleNotFoundError is raised

  @rq-242ef1d8
  Scenario: found(raise_error=True) raises ModuleNotFoundError when cclib is missing
    Given "g16" is present on the PATH
    But "cclib" is not importable
    When GaussianHarness.found(raise_error=True) is called
    Then ModuleNotFoundError is raised

  # --- get_version() ---

  @rq-186c5060
  Scenario: get_version() returns a normalised version string for Gaussian 16
    Given "g16" is present on the PATH
    And running g16 with blank input produces a line "Gaussian 16, Revision C.02,"
    When GaussianHarness.get_version() is called
    Then it returns "16.C.02"

  @rq-850ede4f
  Scenario: get_version() returns a normalised version string for Gaussian 09
    Given only "g09" is present on the PATH
    And running g09 with blank input produces a line "Gaussian 09, Revision D.01,"
    When GaussianHarness.get_version() is called
    Then it returns "09.D.01"

  @rq-3e134642
  Scenario: get_version() caches the result after the first call
    Given "g16" is present on the PATH and has been queried once
    When GaussianHarness.get_version() is called a second time
    Then the Gaussian binary is not re-executed
    And the same version string is returned

  @rq-20d9a784
  Scenario: get_version() raises UnknownError if version cannot be parsed
    Given "g16" is present on the PATH
    And running g16 produces output with no "Gaussian NN, Revision" line
    When GaussianHarness.get_version() is called
    Then UnknownError is raised

  # --- compute() error paths ---

  @rq-c77a2c13
  Scenario: compute() raises InputError when model.basis is a BasisSet object
    Given an AtomicInput whose model.basis is a BasisSet object (not a string)
    When GaussianHarness.compute() is called
    Then InputError is raised without executing Gaussian

  @rq-165224fc
  Scenario: compute() raises InputError when the log contains an unrecognised-keyword error
    Given Gaussian runs and its log contains "Ernie: is not defined as a Gaussian input file"
    When GaussianHarness.compute() is called
    Then InputError is raised

  @rq-edb73b83
  Scenario: compute() raises InputError on internal coordinate error
    Given Gaussian runs and its log contains "Error in internal coordinate system"
    When GaussianHarness.compute() is called
    Then InputError is raised

  @rq-a4c449da
  Scenario: compute() raises ResourceError on memory allocation failure
    Given Gaussian runs and its log contains "galloc: could not allocate memory"
    When GaussianHarness.compute() is called
    Then ResourceError is raised

  @rq-71aceb93
  Scenario: compute() raises UnknownError on abnormal termination with no known pattern
    Given Gaussian runs and its log does not contain "Normal termination"
    And the log does not match any known error pattern
    When GaussianHarness.compute() is called
    Then UnknownError is raised

  # --- compute() success paths ---

  @rq-99807237
  Scenario: compute() returns AtomicResult for a successful energy calculation
    Given a valid AtomicInput for HF/STO-3G energy on water
    And Gaussian terminates normally
    When GaussianHarness.compute() is called
    Then an AtomicResult is returned with success=True
    And return_result is a float representing the total energy in Hartree
    And native_files contains keys "gaussian.com" and "gaussian.log"
    And provenance.creator equals "Gaussian"

  @rq-e62578e2
  Scenario: compute() returns AtomicResult with gradient for gradient driver
    Given a valid AtomicInput for HF/STO-3G gradient on water (3 atoms)
    And Gaussian terminates normally
    When GaussianHarness.compute() is called
    Then an AtomicResult is returned with success=True
    And return_result is a flat list of 9 floats (3 atoms × 3 Cartesian components)

  @rq-4d683dac
  Scenario: compute() returns AtomicResult with Hessian for hessian driver
    Given a valid AtomicInput for HF/STO-3G hessian on water (3 atoms)
    And Gaussian terminates normally
    When GaussianHarness.compute() is called
    Then an AtomicResult is returned with success=True
    And return_result is a flat list of 81 floats (9 × 9 Hessian)

  @rq-1a051bc2
  Scenario: compute() returns AtomicResult for properties driver
    Given a valid AtomicInput for HF/STO-3G properties on water
    And Gaussian terminates normally
    When GaussianHarness.compute() is called
    Then an AtomicResult is returned with success=True
    And return_result is a float (the total energy)
    And extras["qcvars"] contains "CURRENT ENERGY"
```
