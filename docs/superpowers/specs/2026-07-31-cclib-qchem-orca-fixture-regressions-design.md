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

1. **Native input parametrization**: one `pytest.mark.parametrize` table whose rows contain an `AtomicInput`-compatible request and the complete expected generated native input text. The test calls the program-specific input generator and asserts exact text equality. It never invokes the executable.
2. **Parsed output parametrization**: one `pytest.mark.parametrize` table whose rows contain the cclib fixture-relative output path and the saved expected parsed-result dictionary. The test reads the fixture output, passes it through the existing cclib harness conversion path, removes `extras`, and compares the remaining JSON-compatible result data to the embedded expected dictionary. It never invokes the executable.
3. **Live smoke calculation**: one small HF/STO-3G single-point calculation, guarded by the existing program-availability marks. It exercises native input generation, executable invocation, normal-termination handling, cclib parsing, and QCSchema conversion. It is skipped when the corresponding program is not configured.

The expected native text and expected output dictionary are inline arguments in their respective parametrization rows. Every row carries the native-program version derived from its cclib fixture directory (for example `5.1`, `5.4`, or `6.0` for Q-Chem; `5.0`, `6.0`, or `6.1` for ORCA), and its pytest ID includes that version. Tables are grouped by version so parser/generator behavior can be diagnosed and extended per software release. Adding coverage means adding exactly one row to the applicable version group. Output assertions deliberately exclude `extras` for this first suite because harness and parser metadata are not the behavior under test.

## Fixture Eligibility

A fixture is included when it has a paired native input and normal-termination output in a versioned cclib data directory, is detected as the expected cclib parser type, parses successfully, produces a QCSchema result accepted by the current harness, and can form a supported QCSchema request for the generator. Versioned fixture coverage is retained even when the live executable smoke test uses only the supported minimum version.

Fixtures with `.log` rather than `.out`, absent pairs, unsupported drivers, unsupported native features, parser failures, or writer/QCSchema validation failures are recorded in the implementation plan as excluded with their concrete reason. The test suite must not silently discover files at runtime: the parametrization tables are the explicit coverage inventory.

## Snapshot Collection

A collection helper is run once against the installed cclib and the fixed cclib checkout to produce candidate inline dictionaries. It serializes only stable, JSON-compatible fields and omits the top-level `extras` field. The generated data is reviewed and copied into the two test modules; normal tests consume the fixture text and committed expected values, not a collector-generated artifact.

The output tests are skipped with a clear message if `CCLIB_SOURCE_ROOT` is absent or does not contain the named fixture. The expected dictionaries remain committed so parser-behavior changes are visible as a test diff whenever the fixture source is available.

## Acceptance Criteria

- The two requested test files exist and each has the three test shapes above.
- Every parametrized fixture row identifies and is grouped by its Q-Chem or ORCA software version.
- All covered input rows assert byte-for-byte generated input text.
- All covered output rows assert the saved parsed dictionary after excluding `extras`.
- No input/output parametrized test executes Q-Chem or ORCA.
- Each program has one availability-guarded live HF single-point smoke test.
- The tracking plan lists every candidate fixture and its included/excluded status.
- Targeted tests pass with fixture data available; live tests either pass or are explicitly skipped due to missing executable/configuration.
