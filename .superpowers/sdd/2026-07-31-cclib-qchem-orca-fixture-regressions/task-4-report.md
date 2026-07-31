# Task 4 report — ORCA cclib fixture regressions

## Status

Complete. Added the availability-guarded ORCA HF/STO-3G helium end-to-end smoke test while preserving the existing Task 2 ORCA 5.0 and 6.0 fixture expectation tables. ORCA 6.1 has no collector-eligible rows: its only fixture is recorded as excluded in the Task 1 inventory.

## Changed files

- `qcengine/programs/cclib_programs/tests/test_orca.py`
  - Imports `qcengine` and the `uusing` availability guard.
  - Adds `test_live_hf_single_point`, a `cclib-orca` HF/STO-3G He energy calculation guarded by `@uusing("cclib-orca")` and marked `cclib_orca`.
  - Asserts successful execution and a scalar `float` result.

## Validation

```bash
git diff --check
uv run --with black black --check qcengine/programs/cclib_programs/tests/test_orca.py
uv run python -m py_compile qcengine/programs/cclib_programs/tests/test_orca.py
```

Passed. Black found no formatting changes; the test module compiled successfully.

```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib uv run --with /home/awallace43/gits/cclib --with pytest pytest qcengine/programs/cclib_programs/tests/test_orca.py -q
```

Passed: `20 passed, 1 skipped, 10 warnings in 2.67s`. The 20 fixture/input rows passed. The live test skipped because `cclib-orca` is unavailable in this environment, as intended. The warnings are the pre-existing expected Pydantic deprecation warnings from `result.dict(exclude={"extras"})` in the fixture tests.

## Commit

- `3dded64c test: add orca cclib fixture regressions`

## Concerns

- The live ORCA path was not executable in this environment; its availability guard was exercised through the expected skip.
- The pre-existing untracked `uv.lock` remains unmodified and unstaged.

## Review fix: defer live-smoke availability probes

The original `uusing` decorators were eager: importing either new fixture module initialized `qcengine.testing`, which ran native availability probes during pytest collection. Removed `uusing` from both modules and retained their program markers. Each live smoke test now performs `qcng.get_program(<selector>, check=False).found()` inside its test body and calls `pytest.skip` when unavailable. Therefore fixture collection and fixture-only selections do not execute `qchem` or `orca`, while the live tests still skip if their program is unavailable.

### Changed files

- `qcengine/programs/cclib_programs/tests/test_qchem.py`
  - Removes the eager `uusing("cclib-qchem")` decorator/import.
  - Keeps `@pytest.mark.cclib_qchem`; adds the test-body-local `cclib-qchem` availability check and skip.
- `qcengine/programs/cclib_programs/tests/test_orca.py`
  - Removes the eager `uusing("cclib-orca")` decorator/import.
  - Keeps `@pytest.mark.cclib_orca`; adds the test-body-local `cclib-orca` availability check and skip.

### Validation evidence

```bash
uv run --with black black --check qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py
uv run python -m py_compile qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py
git diff --check
```

Passed: both modules are Black-formatted, compile, and have no whitespace errors.

A temporary `PATH` shim supplied failing `qchem` and `orca` executables that append every invocation to a probe log. With that shim first on `PATH`:

```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib uv run --with /home/awallace43/gits/cclib --with pytest pytest qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py --collect-only -q
```

Passed: `118 tests collected`; the probe log was empty.

```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib uv run --with /home/awallace43/gits/cclib --with pytest pytest -q qcengine/programs/cclib_programs/tests/test_qchem.py::test_generated_inputs qcengine/programs/cclib_programs/tests/test_qchem.py::test_parsed_outputs qcengine/programs/cclib_programs/tests/test_orca.py::test_generated_inputs qcengine/programs/cclib_programs/tests/test_orca.py::test_parsed_outputs
```

Passed: `116 passed, 58 warnings in 1.63s`; the probe log was empty. The warnings are the existing Pydantic deprecations from fixture-output serialization.

```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib uv run --with /home/awallace43/gits/cclib --with pytest pytest -q -rs qcengine/programs/cclib_programs/tests/test_qchem.py::test_live_hf_single_point qcengine/programs/cclib_programs/tests/test_orca.py::test_live_hf_single_point
```

Passed: `2 skipped`. Under the same shim, the probe log contained exactly `qchem version.in` and `orca version.inp`, proving probes occur only when the selected live test body runs and that unavailable programs retain skip behavior.
