# CCLib Q-Chem and ORCA Fixture Regression Tests

## Goal

Add dedicated Q-Chem and ORCA regression modules that lock down native-input generation, cclib-backed output conversion, and one live end-to-end HF single-point calculation for each harness.

## Scope

- Create `qcengine/programs/cclib_programs/tests/test_qchem.py`.
- Create `qcengine/programs/cclib_programs/tests/test_orca.py`.
- Cover the remaining compatible paired input/output examples under `/home/awallace43/gits/cclib/data/QChem` and `/home/awallace43/gits/cclib/data/ORCA`.
- Keep `test_cclib.py`'s shared harness/unit tests intact; do not refactor unrelated code.

## Test Design

Each module contains exactly three public regression-test shapes:

1. **Native input parametrization**: one `pytest.mark.parametrize` table of small, representative `AtomicInput` requests. Each row targets one unique generator setting for a supported native-program version and stores only the required regular expressions (for driver mapping, method/basis, resources, coordinates/units, or program-specific keyword blocks). The test calls the program-specific input generator and requires every expression to match. It never invokes the executable.
2. **Parsed output parametrization**: one `pytest.mark.parametrize` table whose rows contain the cclib fixture-relative output path, a dotted QCSchema result path, and a scalar or small numeric-array expectation. The test reads the fixture output, passes it through the existing cclib harness conversion path, retrieves the named field, and compares numerics with `numpy.allclose(atol=1e-6)`. It never invokes the executable.
3. **Live smoke calculation**: one small HF/STO-3G single-point calculation, guarded by the existing program-availability marks. It exercises native input generation, executable invocation, normal-termination handling, cclib parsing, and QCSchema conversion. It is skipped when the corresponding program is not configured.

The input regexes and output `(path, expected)` pairs are inline arguments in their respective parametrization rows. Every row carries the native-program version derived from its cclib fixture directory (for example `5.1`, `5.4`, or `6.0` for Q-Chem; `5.0`, `6.0`, or `6.1` for ORCA), and its pytest ID includes that version. Tables are grouped by version so parser/generator behavior can be diagnosed and extended per software release. Each version includes representative unique settings rather than redundant full-corpus snapshots. Output checks deliberately omit `extras` and compare only fields meaningful to the calculation class. Both deferred live smoke calculations use H2.

## Representative Coverage and Exclusions

Q-Chem 5.1 and 5.4 each retain HF, B3LYP SCF, MP2, and CCSD `properties.return_energy` fixtures; 6.0 retains an HF solvent energy because that is its only eligible fixture class. ORCA 5.0 and 6.0 each retain HF, MP2, and CCSD energies. ORCA DFT outputs and the only 6.1 output are excluded because the current cclib writer fails their QCSchema conversion with `KeyError: 'functional'`. Other eligible files are intentionally excluded as redundant corpus coverage. The implementation plan is the complete eligibility and exclusion inventory.

## Fixture Eligibility

A fixture is included when it has a paired native input and normal-termination output in a versioned cclib data directory, is detected as the expected cclib parser type, parses successfully, produces a QCSchema result accepted by the current harness, and can form a supported QCSchema request for the generator. Versioned fixture coverage is retained even when the live executable smoke test uses only the supported minimum version.

Fixtures with `.log` rather than `.out`, absent pairs, unsupported drivers, unsupported native features, parser failures, writer/QCSchema validation failures, or redundant eligible representatives are recorded in the implementation plan as excluded with their concrete reason. The test suite must not silently discover files at runtime: the parametrization tables are the explicit compact coverage inventory.

## Snapshot Collection

A collection helper is run once against the installed cclib and the fixed cclib checkout to produce candidate inline dictionaries. It serializes only stable, JSON-compatible fields and omits the top-level `extras` field. The generated data is reviewed and copied into the two test modules; normal tests consume the fixture text and committed expected values, not a collector-generated artifact.

The output tests are skipped with a clear message if `CCLIB_SOURCE_ROOT` is absent or does not contain the named fixture. The expected dictionaries remain committed so parser-behavior changes are visible as a test diff whenever the fixture source is available.

## Acceptance Criteria

- The two requested test files exist and each has the three test shapes above.
- Every parametrized fixture row identifies and is grouped by its Q-Chem or ORCA software version.
- All covered input rows assert their required native fragments with regular expressions.
- All covered output rows assert a calculation-specific scalar or small array at its saved QCSchema path with an absolute tolerance of `1e-6`.
- No input/output parametrized test executes Q-Chem or ORCA.
- Each program has one availability-guarded live HF single-point smoke test.
- The tracking plan lists every candidate fixture's eligibility and the compact representative inclusion/exclusion decision.
- Targeted tests pass with fixture data available; live tests either pass or are explicitly skipped due to missing executable/configuration.
