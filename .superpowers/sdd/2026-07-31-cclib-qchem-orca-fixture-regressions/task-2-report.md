# Task 2 report — inline expectation scaffolding

## Status

Complete within Task 2 scope. Added one exact-input regression row per program for a minimal HF/STO-3G helium energy request. The rows are self-contained and invoke only `build_input`; they do not discover or invoke `qchem` or `orca`.

The Task 1 inventory is the authoritative fixture-status source. It records 58 conversion-eligible paired fixtures (48 Q-Chem and 10 ORCA) and their exclusions. Per the Task 2 boundary, those program/version fixture rows, parsed-output assertions, and live executable tests were not added; Tasks 3 and 4 own them.

## Changed files

- `qcengine/programs/cclib_programs/tests/test_qchem.py`
  - New Q-Chem exact-input module with one `"minimum"` HF/STO-3G He row.
  - Pins the generated Q-Chem input literal, including the one-core 1024 MB resource value.
- `qcengine/programs/cclib_programs/tests/test_orca.py`
  - New ORCA exact-input module with one `"minimum"` HF/STO-3G He row.
  - Pins the generated ORCA input literal, including `%pal`, `%MaxCore`, and default output settings.

## TDD and collection method

1. Each module was first created with the requested `("minimum", HE_INPUT, "")` parameter row.
2. After supplying all required `TaskConfig` fields used by this checkout (`ncores=1`, `nnodes=1`, `memory=1.0`, `scratch_directory=None`, `retries=0`, and `mpiexec_command=None`), each focused test failed because generated native text differed from the intentionally empty expected literal.
3. A local one-shot Python collector constructed the same `AtomicInput`, constructed that fixed `TaskConfig`, and called each program's `build_input(input_model, config, "/resolved/program")`. Its `repr(job.input_text)` output was copied exactly into the corresponding test row as adjacent Python string literals.
4. The focused tests then passed together without executable invocation.

The brief's complete-fixture `ExecutionResult`/`_parse_and_convert` JSON collection is intentionally not materialized in these modules: Task 3 begins and owns Q-Chem parsed-output rows, and Task 4 begins and owns ORCA parsed-output rows. Task 1 already ran that collector against every authoritative included fixture to establish the statuses; see `task-1-report.md`.

## Commands and results

1. Required commands as written:

   ```bash
   pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_generated_inputs -q
   pytest qcengine/programs/cclib_programs/tests/test_orca.py::test_generated_inputs -q
   ```

   Result: neither could start because `pytest` is not on `PATH` (`/bin/bash: pytest: command not found`).

2. Repository-managed equivalent while establishing the failing rows:

   ```bash
   uv run --with pytest pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_generated_inputs -q
   uv run --with pytest pytest qcengine/programs/cclib_programs/tests/test_orca.py::test_generated_inputs -q
   ```

   Result: both initially exposed that this checkout's `TaskConfig` requires explicit `nnodes`, `scratch_directory`, `retries`, and `mpiexec_command`; after those fixed configuration values were supplied, each test failed as intended on `expected_input == ""` (one failure per module).

3. Exact native-input collector:

   ```bash
   uv run python - <<'PY'
   # construct the module's HF/STO-3G He AtomicInput and fixed TaskConfig,
   # then print repr(builder(input_model, config, "/resolved/program").input_text)
   # for cclib_qchem.build_input and cclib_orca.build_input
   PY
   ```

   Result: passed. The returned Q-Chem and ORCA strings were copied verbatim into their `"minimum"` rows.

4. Focused regression validation:

   ```bash
   uv run --with pytest pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_generated_inputs qcengine/programs/cclib_programs/tests/test_orca.py::test_generated_inputs -q
   ```

   Result: passed, `2 passed in 0.36s`.

5. Whitespace and worktree check before the test commit:

   ```bash
   git diff --check
   git status --short
   ```

   Result: `git diff --check` passed. Only the two new test modules and pre-existing untracked `uv.lock` were shown.

## Commits

- `bfa56a39 test: scaffold cclib qchem and orca fixtures`
- `docs: report cclib fixture scaffolding` (this report)

## Concerns

- `uv.lock` was already untracked before Task 2 work and remains deliberately unmodified and unstaged.
- Fixture-specific input/output rows and complete parsed-result dictionaries remain deferred to Tasks 3 and 4 to avoid widening this scaffolding task.

---

## Task 2 gap fix

### Status

Implemented the reviewer-requested complete fixture expectation coverage. The two modules now preserve exact generated native input text and serialized parsed-result dictionaries for every Task 1 included fixture: 48 Q-Chem rows (5.1, 5.4, 6.0) and 10 ORCA rows (5.0, 6.0). The rows are statically grouped by native fixture version and do not invoke Q-Chem or ORCA.

`test_parsed_outputs` reads each fixture only from `$CCLIB_SOURCE_ROOT/data`, creates an `ExecutionResult` containing the complete fixture output, replays it through `_parse_and_convert`, and compares the sorted JSON serialization of `result.dict(exclude={"extras"})`. The serialization helper converts only NumPy arrays to their JSON representation; `extras` is the sole excluded result field.

### Files

- `qcengine/programs/cclib_programs/tests/test_qchem.py`
  - Adds `QCHEM_5_1`, `QCHEM_5_4`, and `QCHEM_6_0` input/output expectation groups.
  - Adds 48 exact-input rows and 48 full `ExecutionResult` replay/output rows.
- `qcengine/programs/cclib_programs/tests/test_orca.py`
  - Adds `ORCA_5_0` and `ORCA_6_0` input/output expectation groups.
  - Adds 10 exact-input rows and 10 full `ExecutionResult` replay/output rows.

### Commands and results

```bash
uv run --with /home/awallace43/gits/cclib python /tmp/generate_cclib_fixture_tests.py
```

Passed. Replayed every Task 1 included fixture and generated inline `AtomicInput`, native input, and sorted serialized expected result values.

```bash
uv run --with black black qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py
uv run python -m py_compile qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py
```

Passed. Both generated modules are formatted and compile successfully.

```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib uv run --with /home/awallace43/gits/cclib --with pytest pytest qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py -q
```

Passed: `116 passed in 1.59s`. The run emitted 58 expected Pydantic deprecation warnings for `result.dict`; that method is intentionally retained because the task requires serialization of `result.dict(exclude={"extras"})`.

```bash
git diff --check
```

Passed.

### Commit

- `038cc06e test: add cclib fixture expectations`

### Concerns

- The fixture tests require the local cclib checkout through `CCLIB_SOURCE_ROOT` and the same cclib writer/parser revision used by Task 1. They skip clearly when that root is unavailable.
- The pre-existing untracked `uv.lock` remains unmodified and unstaged.
