# Task 3 report — Q-Chem fixture regressions

## Status

Complete. The versioned Q-Chem 5.1, 5.4, and 6.0 fixture input/output tables populated by Task 2 were preserved. Task 3 adds the remaining availability-guarded, full end-to-end Q-Chem HF/STO-3G helium single-point smoke test.

## Changed files

- `qcengine/programs/cclib_programs/tests/test_qchem.py`
  - Imports the public `qcengine.compute` API and `uusing` availability guard.
  - Adds `test_live_hf_single_point`, marked for `cclib_qchem` and skipped when `cclib-qchem` is unavailable.
  - Runs an HF/STO-3G helium energy calculation through `cclib-qchem` with one core and 1 GiB memory, then asserts `success is True` and that `return_result` is a scalar `float`.

## Validation

```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_qchem.py -q
```

Could not start because `pytest` is not on `PATH` (`/bin/bash: pytest: command not found`).

```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib uv run --with /home/awallace43/gits/cclib --with pytest pytest qcengine/programs/cclib_programs/tests/test_qchem.py -q
```

Passed: `96 passed, 1 skipped, 48 warnings in 3.28s`. All 48 generated-input and 48 parsed-output fixture rows passed. The live smoke test skipped because Q-Chem/cclib-qchem is unavailable in this environment. The 48 existing Pydantic deprecation warnings come from the required fixture serialization using `result.dict(exclude={"extras"})`.

```bash
uv run --with pytest pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_live_hf_single_point -q -rs
```

Passed: `1 skipped`; pytest reported `Not detecting module cclib-qchem`, confirming the availability guard prevents an unavailable executable from failing the suite.

```bash
uv run --with black black --check qcengine/programs/cclib_programs/tests/test_qchem.py
uv run python -m py_compile qcengine/programs/cclib_programs/tests/test_qchem.py
git diff --check
```

Passed. The test file is Black-formatted, compiles, and has no whitespace errors.

## Commit

- `6130355c test: add qchem cclib fixture regressions`

## Concerns

- Q-Chem was unavailable locally, so the new live path was correctly skipped rather than executed. It must execute in an environment with a configured `cclib-qchem` installation.
- The pre-existing untracked `uv.lock` remains unmodified and unstaged.
