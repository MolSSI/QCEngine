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
