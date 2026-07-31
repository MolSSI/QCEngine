# Task 1 fixture inventory report

## Status

Complete. Task 1 only was implemented: the plan now contains the complete approved-version fixture inventory and collector results. No Q-Chem or ORCA test module was created or modified.

## Files changed

- `docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md`
  - Marked all three Task 1 steps complete.
  - Recorded every paired `.out` candidate under Q-Chem 5.1/5.4/6.0 and ORCA 5.0/6.0/6.1 with independent `included-input`/`included-output` statuses or exact failure stages.
  - Recorded all 18 `.log`-only candidates as excluded and recorded that no `.out` candidate lacks its paired native input.

No test files were added or updated, per task scope.

## Commands and results

1. Prescribed inventory command:

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

   Result: the checkout has no `python` executable on `PATH` (`/bin/bash: python: command not found`).

2. Same prescribed inventory body executed with the repository tool environment:

   ```bash
   uv run --with /home/awallace43/gits/cclib python - <<'PY'
   # prescribed inventory body above
   PY
   ```

   Result: passed; 91 `.out` candidates enumerated, all `paired` (55 Q-Chem, 36 ORCA).

3. Collector check:

   ```bash
   uv run --with /home/awallace43/gits/cclib python - <<'PY'
   # one-shot collector: ccopen, exact parser class, parse,
   # QCSchemaWriter(...).as_dict(validate=False), AtomicInput construction,
   # native input generation, and _parse_and_convert for each paired .out
   PY
   ```

   Result: passed for all 91 paired candidates. It used cclib `1.8.1.post1235+5c3639de` from `/home/awallace43/gits/cclib` and `TaskConfig(ncores=1, nnodes=1, memory=1.0, scratch_directory=None, retries=0, mpiexec_command=None)`. Exact per-fixture collector outcomes are in the plan.

4. Inventory-to-plan audit:

   ```bash
   python3 - <<'PY'
   # assert every JSONL collector row and every .log-only exclusion is present
   # in docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md
   PY
   ```

   Result: passed; validated 91 paired rows and 18 `.log`-only exclusions against the collector output.

5. Whitespace validation:

   ```bash
   git diff --check
   ```

   Result: passed before the baseline commit.

6. Baseline commit:

   ```bash
   git add -f docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md
   git commit -m "docs: track cclib fixture regressions"
   ```

   Result: passed; commit `e1eecb9b docs: track cclib fixture regressions`.

## Fixture counts and statuses

| Program/version | Paired `.out` | Included inputs | Included outputs | Excluded `.out` | Excluded `.log`-only | Missing-pair |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Q-Chem 5.1 | 21 | 17 | 17 | 4 | 0 | 0 |
| Q-Chem 5.4 | 21 | 18 | 18 | 3 | 0 | 0 |
| Q-Chem 6.0 | 13 | 13 | 13 | 0 | 0 | 0 |
| ORCA 5.0 | 17 | 4 | 4 | 13 | 9 | 0 |
| ORCA 6.0 | 18 | 6 | 6 | 12 | 9 | 0 |
| ORCA 6.1 | 1 | 0 | 0 | 1 | 0 | 0 |
| **Total** | **91** | **58** | **58** | **33** | **18** | **0** |

The 33 paired `.out` exclusions are exactly recorded in the plan. Their collector-stage causes are: one Q-Chem parser result with `metadata.success` not true; six Q-Chem writer failures for CCD, CCSD(T), or MP4; 22 ORCA writer `KeyError: 'functional'` failures; and four ORCA writer failures for CCSD(T) or MP3.

## Commits

- `e1eecb9b docs: track cclib fixture regressions` — Task 1 inventory baseline.
- `docs: report cclib fixture inventory` — this detailed Task 1 report follow-up commit.

## Concerns

- The released `cclib` package resolved by `uv run --with cclib` lacks `cclib.io.qcschemawriter`, which the current harness imports. The collector therefore correctly used the required local cclib checkout, whose source provides that module.
- Collector exclusions are intentionally baseline facts for the current local cclib checkout. A later cclib update may change the writer outcomes and should prompt inventory recollection.
