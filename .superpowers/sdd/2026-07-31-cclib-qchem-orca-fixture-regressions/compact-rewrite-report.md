# Compact cclib fixture rewrite report

## Scope and result

Rewrote only the two dedicated cclib fixture suites and their existing design/plan records. The suites now use compact representative `AtomicInput` rows with regex fragment assertions and compact `(version, relative_output, dotted_path, expected)` output rows. Output values are read from converted QCSchema results and checked with `numpy.allclose(atol=1e-6)`.

- Q-Chem output coverage: 5.1 and 5.4 HF/DFT/MP2/CCSD energies; 6.0 eligible HF solvent energy.
- ORCA output coverage: 5.0 and 6.0 HF/MP2/CCSD energies. ORCA DFT conversion and the sole 6.1 fixture remain excluded because cclib's writer raises `KeyError: 'functional'`.
- Live smoke tests still defer availability checking until their test bodies. Both use H2 HF/STO-3G; ORCA no longer uses He.
- No cclib checkout files were changed (`git -C /home/awallace43/gits/cclib status --short` was empty).

## Changed files

- `qcengine/programs/cclib_programs/tests/test_qchem.py`
- `qcengine/programs/cclib_programs/tests/test_orca.py`
- `docs/superpowers/specs/2026-07-31-cclib-qchem-orca-fixture-regressions-design.md`
- `docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md`

## Line counts

| File | Before | After | Reduction |
| --- | ---: | ---: | ---: |
| `test_qchem.py` | 35,937 | 164 | 35,773 |
| `test_orca.py` | 11,765 | 161 | 11,604 |
| **Total** | **47,702** | **325** | **47,377** |

## Commands and results

Activated environment for all Python validation:

```bash
source /home/awallace43/miniconda3/etc/profile.d/conda.sh
conda activate cclib_qcng
```

```bash
python -m py_compile qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py
black --check qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py
```

Passed: both modules compiled and Black reported no changes.

```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest -q qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py
```

Passed: `21 passed, 2 skipped in 1.10s`. The two deferred live smoke tests skipped because native executables were unavailable.

```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest -q qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py qcengine/programs/cclib_programs/tests/test_cclib.py
```

Passed: `90 passed, 7 skipped in 3.94s`.

```bash
git diff --check
git -C /home/awallace43/gits/cclib status --short
```

Passed: no whitespace errors; cclib worktree remained clean.

## Diff and commit

Implementation diff before commit: 365 insertions and 47,860 deletions across the two rewritten suites and their specification/plan updates. The compact rewrite is committed as:

- `e255f03e test: compact cclib fixture regressions`

## Residual risk

Live Q-Chem and ORCA execution paths were skipped due to unavailable executables. Fixture conversion, generation, and shared cclib harness tests passed against the activated `cclib_qcng` environment and fixed cclib source root.
