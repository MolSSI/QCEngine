# CCLib Q-Chem and ORCA Fixture Regressions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add version-grouped, inline Q-Chem and ORCA input/output regression tables and one live HF smoke calculation per harness.

**Architecture:** Each dedicated module owns an explicit fixture inventory as parametrization rows. Input rows compare the generator's complete native text without native execution. Output rows replay cclib fixture output through the existing conversion path and compare saved dictionaries excluding `extras`; live availability-marked rows alone invoke native software.

**Tech Stack:** Python, pytest, QCEngine, QCElemental, cclib fixture corpus.

## Global Constraints

- Fixture rows are grouped and identified by the native version encoded in their cclib data directory: Q-Chem 5.1/5.4/6.0 and ORCA 5.0/6.0/6.1.
- Input and output expected values are inline `pytest.mark.parametrize` arguments.
- Input/output parametrized tests must never execute Q-Chem or ORCA.
- Output comparisons remove only `extras`.
- Output fixtures are loaded from `$CCLIB_SOURCE_ROOT/data`; unavailable fixture roots skip clearly.
- Live HF/STO-3G single-point tests use the existing program availability marks.

---

### Task 1: Build and record the fixture inventory

**Files:**
- Create: `docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md`

**Interfaces:**
- Consumes: `/home/awallace43/gits/cclib/data/{QChem,ORCA}`.
- Produces: explicit per-version inclusion/exclusion inventory used to form test rows.

- [ ] **Step 1: Enumerate paired fixture candidates and parser outcomes**

Run:
```bash
python - <<'PY'
from pathlib import Path
for program in ("QChem", "ORCA"):
    root = Path("/home/awallace43/gits/cclib/data") / program
    for output in sorted(root.glob("**/*.out")):
        source = output.with_suffix(".in" if program == "QChem" else ".inp")
        print(program, output.parent.name, output.stem, "paired" if source.exists() else "missing-input")
PY
```
Record every row under its version directory. Mark `.log`-only and missing-pair candidates excluded.

- [ ] **Step 2: Record collector failures as exclusions**

For every paired `.out`, run `ccopen`, the expected parser class check, parse, `QCSchemaWriter(...).as_dict(validate=False)`, and `_parse_and_convert`. Record `included-input`, `included-output`, or the exact failed stage.

- [ ] **Step 3: Commit the inventory baseline**

```bash
git add -f docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md
git commit -m "docs: track cclib fixture regressions"
```

### Task 2: Collect inline expectations

**Files:**
- Create: `qcengine/programs/cclib_programs/tests/test_qchem.py`
- Create: `qcengine/programs/cclib_programs/tests/test_orca.py`

**Interfaces:**
- Consumes: fixture path, `AtomicInput`, `TaskConfig`, and program `build_input`/`_parse_and_convert` APIs.
- Produces: `(version, input_model, expected_input)` and `(version, relative_output, expected_output)` rows.

- [ ] **Step 1: Write a failing exact-input row**

Start each module with one minimal HF/STO-3G He row:
```python
@pytest.mark.parametrize("version,input_model,expected_input", [("minimum", HE_INPUT, "")])
def test_generated_inputs(version, input_model, expected_input):
    assert build_input(input_model, TASK_CONFIG, "/resolved/program").input_text == expected_input
```

- [ ] **Step 2: Verify the test fails because expected text is empty**

Run:
```bash
pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_generated_inputs -q
pytest qcengine/programs/cclib_programs/tests/test_orca.py::test_generated_inputs -q
```
Expected: each fails with generated native text differing from `""`.

- [ ] **Step 3: Collect exact native input and parsed dictionaries**

Use a local one-shot Python collector to create `AtomicInput` from the cclib writer result, call the program generator with fixed `TaskConfig(ncores=1, memory=1.0)`, replay the complete `.out` through `ExecutionResult` and `_parse_and_convert`, then serialize `result.dict(exclude={"extras"})` using `json.dumps(..., sort_keys=True)`. Preserve the input text and dictionary as Python literals in their corresponding version group.

- [ ] **Step 4: Verify the exact input tests pass without executables**

Run:
```bash
pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_generated_inputs qcengine/programs/cclib_programs/tests/test_orca.py::test_generated_inputs -q
```
Expected: PASS without `qchem` or `orca` invocation.

### Task 3: Implement Q-Chem regressions

**Files:**
- Modify: `qcengine/programs/cclib_programs/tests/test_qchem.py`

**Interfaces:**
- Consumes: Q-Chem 5.1, 5.4, and 6.0 fixture paths and saved expectation rows.
- Produces: `test_generated_inputs`, `test_parsed_outputs`, and `test_live_hf_single_point`.

- [ ] **Step 1: Add a failing output assertion for `basicQChem5.1/water_mp2.out`**

```python
@pytest.mark.parametrize("version,relative_output,expected_output", [
    ("5.1", "QChem/basicQChem5.1/water_mp2.out", {}),
])
def test_parsed_outputs(version, relative_output, expected_output):
    assert normalize(convert_fixture(relative_output)) == expected_output
```

- [ ] **Step 2: Run it with fixture data and verify it fails on `{}`**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_parsed_outputs -q
```
Expected: FAIL because the converted result contains fields absent from `{}`.

- [ ] **Step 3: Add all collected Q-Chem version groups and live smoke test**

Populate `5.1`, `5.4`, and `6.0` rows using collector output. Add an availability-guarded `cclib-qchem` HF/STO-3G He (or H2) energy calculation asserting `success` and scalar `return_result`.

- [ ] **Step 4: Run Q-Chem regression tests**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_qchem.py -q
```
Expected: fixture rows PASS; live row PASS or SKIP based on program availability.

- [ ] **Step 5: Commit Q-Chem coverage**

```bash
git add qcengine/programs/cclib_programs/tests/test_qchem.py
git commit -m "test: add qchem cclib fixture regressions"
```

### Task 4: Implement ORCA regressions

**Files:**
- Modify: `qcengine/programs/cclib_programs/tests/test_orca.py`

**Interfaces:**
- Consumes: ORCA 5.0, 6.0, and 6.1 fixture paths and saved expectation rows.
- Produces: `test_generated_inputs`, `test_parsed_outputs`, and `test_live_hf_single_point`.

- [ ] **Step 1: Add a failing output assertion for `basicORCA6.0/dvb_sp_hf.out`**

```python
@pytest.mark.parametrize("version,relative_output,expected_output", [
    ("6.0", "ORCA/basicORCA6.0/dvb_sp_hf.out", {}),
])
def test_parsed_outputs(version, relative_output, expected_output):
    assert normalize(convert_fixture(relative_output)) == expected_output
```

- [ ] **Step 2: Verify it fails on the empty expected dictionary**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_orca.py::test_parsed_outputs -q
```
Expected: FAIL because the converted result contains fields absent from `{}`.

- [ ] **Step 3: Add all collected ORCA version groups and live smoke test**

Populate `5.0`, `6.0`, and `6.1` rows. Include generator-supported energy requests in input rows; record other-driver input exclusions in the inventory. Add an availability-guarded `cclib-orca` HF/STO-3G He energy calculation asserting `success` and scalar `return_result`.

- [ ] **Step 4: Run ORCA regression tests**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_orca.py -q
```
Expected: fixture rows PASS; live row PASS or SKIP based on program availability.

- [ ] **Step 5: Commit ORCA coverage**

```bash
git add qcengine/programs/cclib_programs/tests/test_orca.py
git commit -m "test: add orca cclib fixture regressions"
```

### Task 5: Verify suite boundaries and complete tracking

**Files:**
- Modify: `docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md`
- Test: `qcengine/programs/cclib_programs/tests/test_qchem.py`
- Test: `qcengine/programs/cclib_programs/tests/test_orca.py`
- Test: `qcengine/programs/cclib_programs/tests/test_cclib.py`

- [ ] **Step 1: Mark each completed or excluded inventory row**

For every candidate, retain its version, source/output path, final status, and concrete exclusion reason if omitted.

- [ ] **Step 2: Run the targeted suite with fixture data**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py qcengine/programs/cclib_programs/tests/test_cclib.py -q
```
Expected: PASS, with only availability-guarded live tests skipped when native software is unavailable.

- [ ] **Step 3: Commit tracking and verification updates**

```bash
git add -f docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md
git commit -m "docs: complete cclib fixture regression plan"
```
