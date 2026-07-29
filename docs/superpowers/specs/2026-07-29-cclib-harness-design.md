# Configurable cclib Program Harness Design

## Summary

Add a configurable `CCLibHarness` to QCEngine so multiple external quantum-chemistry programs can share one execution-and-parsing implementation. Register two independent public selectors:

- `cclib-qchem` for Q-Chem 5.1 and later
- `cclib-orca` for ORCA 6.0 and later

The existing native `qchem` harness remains registered and unchanged. QCEngine generates native input from QCSchema `AtomicInput`, runs the selected executable, passes its text output to cclib, converts cclib data through `QCSchemaWriter`, and returns a QCEngine-compatible `AtomicResult`.

The first release covers atomic energy, gradient, and Hessian calculations for HF, B3LYP, BP86, MP2, and CCSD with string basis names. It uses representative input and output behavior from:

- `/home/awallace43/gits/cclib/data/QChem/basicQChem5.1/`
- `/home/awallace43/gits/cclib/data/ORCA/basicORCA6.0/`

The required demonstration is a live Q-Chem MP2/STO-3G water calculation through `qcengine.compute(..., "cclib-qchem")` that reproduces the structure, core numerical results, and representative cclib extras of `/home/awallace43/projects/cclib_qcng_example/qcschema_mp2.json` within stated tolerances.

## Goals

1. Provide cclib-backed QCEngine access to Q-Chem without displacing the existing native Q-Chem harness.
2. Add ORCA support without creating a standalone ORCA parser in QCEngine.
3. Prove that one `CCLibHarness` implementation can be instantiated and registered for multiple programs.
4. Generate native inputs from QCSchema rather than accepting fixture output or requiring raw native input.
5. Preserve the rich flat cclib attributes emitted by `QCSchemaWriter` in `AtomicResult.extras`.
6. Produce explicit, actionable errors for program discovery, Q-Chem environment configuration, execution, parsing, and schema conversion failures.
7. Keep cclib optional and avoid machine-specific paths or setup commands in reusable harness code.

## Non-goals

The initial implementation does not:

- Replace or modify the native `qchem` selector.
- Redesign QCEngine's one-name-per-harness registry.
- Modify cclib source code.
- Accept raw native input through `AtomicInput.extras` or another bypass.
- Support geometry optimization, relaxed or unrelaxed scans, molecular dynamics, excited-state procedures, NMR, Raman, polarizability, or other specialized property drivers.
- Support MP3, MP4, CCSD(T), ADC, EOM-CCSD, STEOM, ROCIS, or methods the current cclib `QCSchemaWriter` cannot represent.
- Support ghost atoms, structured QCSchema `BasisSet` objects, multi-job files, restart files, or node-parallel execution.
- Guarantee byte-for-byte equality with historical fixture JSON across program builds.
- Encode local executable paths or source shell initialization scripts automatically.

If implementation reveals a cclib defect that prevents an approved case, work stops and records the concrete gap before any cclib change or scope expansion.

## Existing constraints

QCEngine's registry reads one `ProgramHarness.name` value from each instance and stores the instance under the lowercase name. Dispatch does not pass the requested registry name to `compute()`. Therefore, a list-valued alias API would require changes across registration, dispatch, availability, CLI behavior, and tests.

The design avoids that change by creating two configured instances of one class:

```python
CCLibHarness(name="cclib-qchem", program="qchem")
CCLibHarness(name="cclib-orca", program="orca")
```

Both instances are registered statically in `qcengine/programs/base.py` like other built-in program harnesses.

cclib auto-detects Q-Chem from `A Quantum Leap Into The Future Of Chemistry` and ORCA from `O   R   C   A`. Its current development branch converts parsed data to QCSchema v1 through `cclib.io.qcschemawriter.QCSchemaWriter`. QCEngine harnesses operate internally on QCSchema v2, so conversion must cross a validated v1 `AtomicResult` before returning a v2 model to QCEngine.

## Architecture

### `CCLibHarness`

Create `qcengine/programs/cclib.py` with a frozen `ProgramHarness` subclass named exactly `CCLibHarness`. Each instance has a `program` discriminator whose supported values are `qchem` and `orca`. Common harness defaults are:

- `scratch=True`
- `thread_safe=False`
- `thread_parallel=True`
- `node_parallel=False`
- `managed_memory=True`

The class owns the shared lifecycle:

1. Check cclib and executable availability.
2. Probe and cache executable identity/version.
3. Validate the approved `AtomicInput` subset.
4. Dispatch to the selected native input generator.
5. Execute in QCEngine-managed scratch space.
6. Select the program output text.
7. Parse through cclib auto-detection.
8. Verify the detected parser matches the configured program.
9. Convert through `QCSchemaWriter` and QCElemental.
10. Validate requested-versus-parsed identity and attach QCEngine/cclib metadata.

### Program definitions

A private immutable program-definition table maps each `program` discriminator to:

- Registry selector
- Executable basename
- Minimum version
- Input and output filenames
- cclib parser class expected from auto-detection
- Version/identity probe
- Native input generator
- Command builder
- Normal-termination marker
- Output-text selector
- Program-specific environment preflight

The table contains data and narrow callables; it does not create another public registry. Adding a future cclib-supported program requires one definition, one input generator, one configured `CCLibHarness` instance, and focused tests.

### Registration

`qcengine/programs/base.py` imports `CCLibHarness` and registers:

```python
register_program(CCLibHarness(name="cclib-qchem", program="qchem"))
register_program(CCLibHarness(name="cclib-orca", program="orca"))
```

`list_all_programs()` always includes both selectors. `list_available_programs()` includes each selector only when that configured instance's cclib and executable checks succeed.

## Supported input contract

### Drivers

The first implementation accepts only:

| QCSchema driver | Q-Chem job type | ORCA simple keyword |
| --- | --- | --- |
| `energy` | `sp` | none |
| `gradient` | `force` | `engrad` |
| `hessian` | `freq` | `freq` |

Every other driver raises `InputError` before execution.

### Methods

The initial accepted method names, compared case-insensitively, are:

- `hf`
- `b3lyp`
- `bp86`
- `mp2`
- `ccsd`

Native restricted/unrestricted behavior follows charge, multiplicity, program defaults, and explicit non-reserved program keywords. Parsed aliases are normalized only for result consistency checks:

- `rhf` and `uhf` compare as `hf`
- `rmp2` and `ump2` compare as `mp2`
- `rccsd` and `uccsd` compare as `ccsd`

Other methods raise `InputError`. Extending this allowlist requires a cclib-writer-compatible fixture or live test.

### Basis and molecule

- `AtomicInput.specification.model.basis` must be a non-empty string.
- Structured `BasisSet` values are rejected.
- All atoms must be real; ghost atoms are rejected.
- Charge and multiplicity come from `AtomicInput.molecule`.
- Atom order is preserved.
- The harness accepts the QCSchema Cartesian geometry in bohr and emits the units required by the target input format.

## Native input generation

### Q-Chem

The Q-Chem generator emits:

1. `$comment` identifying QCEngine and `CCLibHarness`.
2. `$molecule` with charge, multiplicity, and Cartesian coordinates.
3. `$rem` containing driver, method, basis, resource, harvesting, and user options.

Geometry is emitted in bohr and `$rem` includes `INPUT_BOHR=TRUE`.

Derived options are authoritative:

- `JOBTYPE` from the driver table
- `METHOD` from the QCSchema model
- `BASIS` from the QCSchema model
- `MEM_TOTAL` as `int(TaskConfig.memory * 1024)` MiB
- `INPUT_BOHR=TRUE`

Default non-energy-changing harvesting options are:

- `SCF_FINAL_PRINT=2`
- `PRINT_GENERAL_BASIS=TRUE`
- `PRINT_ORBITALS=TRUE`
- `MOLDEN_FORMAT=FALSE`

Ordinary `AtomicInput.specification.keywords` keys become uppercase `$rem` keys. Keys are emitted deterministically. The derived keys above are reserved; attempting to provide any reserved key raises `InputError` rather than silently overriding either QCSchema or user intent.

The execution command is equivalent to:

```text
qchem -nt <ncores> dispatch.in dispatch.out
```

The resolved executable path replaces `qchem` in the actual command.

### ORCA

The ORCA generator emits:

1. A `!` line containing method, basis, mapped driver, and additional simple keywords.
2. Optional `%...` blocks.
3. `%pal` with `nprocs TaskConfig.ncores`.
4. `%MaxCore` with per-process memory in MiB, calculated as `max(1, int(TaskConfig.memory * 1024 / TaskConfig.ncores))`.
5. `* xyz <charge> <multiplicity>` Cartesian coordinates and closing `*`.

Geometry is converted from bohr to Angstrom.

ORCA-specific `AtomicInput.specification.keywords` has this exact shape:

```python
{
    "simple": ["rks", "usesym"],
    "blocks": {
        "output": "PrintLevel Normal\nPrint[P_Basis] 2"
    },
}
```

- `simple` must be a list of non-empty strings without newlines.
- `blocks` must map valid block names to string bodies.
- User block bodies are emitted after harness defaults for the same block, so later user lines have ORCA's normal last-value precedence.
- The `pal` and `maxcore` block names are reserved because resources come from `TaskConfig`.
- Coordinates are not accepted through keywords.
- Unknown top-level ORCA keyword keys raise `InputError`.

The default `%output` body requests data useful to cclib:

```text
PrintLevel Normal
Print[P_Basis] 2
Print[P_MOs] 1
Print[P_Overlap] 1
Print[P_Hirshfeld] 1
```

The execution command is equivalent to:

```text
orca dispatch.inp
```

ORCA's captured standard output is the parseable output text.

## Availability and executable identity

### General rules

The harness discovers programs only through `PATH` using QCEngine/QCElemental discovery utilities. It does not accept harness-specific absolute executable paths. Agents, users, schedulers, and test environments must configure `PATH` and program environment variables before importing or running availability-dependent code.

`found(raise_error=False)` returns `False` for:

- cclib not importable
- the cached cclib/QCElemental compatibility probe fails
- executable absent from `PATH`
- executable identity mismatch
- unparseable or unsupported program version
- failed program-specific environment preflight

With `raise_error=True`, the harness raises a specific `ResourceError` describing the failed resource check.

The identity probe is required because a PATH-resolved executable named `orca` may be an unrelated program. A basename match alone is insufficient.

A cached cclib compatibility probe constructs a minimal synthetic HF `ccData`, calls `QCSchemaWriter(...).as_dict(validate=False)`, and validates the dictionary as a QCElemental QCSchema v1 `AtomicResult`. It also checks the expected Angstrom-to-bohr geometry conversion and presence of flat cclib extras. This rejects an installed cclib whose writer predates the required development behavior without relying on a local path or unreleased version number.

### Version probes

- Q-Chem uses the existing lightweight version-input pattern and requires a Q-Chem banner plus a parseable version at least 5.1.
- ORCA uses a minimal cached probe that produces the ORCA banner and `Program Version` line, and requires a parseable version at least 6.0.
- Probe results are cached by resolved executable path, not only by registry name or class.

The ORCA probe may execute a minimal one-atom HF/STO-3G input when a no-calculation version flag cannot prove identity. This probe runs at most once per resolved executable path in a process.

### Local verification setup

Reusable code contains no paths from this section. On the target machine, live verification prepares the shell with:

```bash
export PATH=/projects/cos-lab-cs207/common/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
source ~/qchem_vars.sh
```

After sourcing, Q-Chem should resolve to:

```text
/projects/cos-lab-cs207/common/software/qchem5.1/bin/qchem
```

ORCA should resolve to:

```text
/projects/cos-lab-cs207/common/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca
```

The Python environment used for verification installs both source trees editably and satisfies QCEngine's QCElemental version range. The cclib validation baseline is the `qcng` branch at commit `467ae3f7` or a descendant containing the same QCSchema-writer behavior:

```bash
python -m pip install -e /home/awallace43/gits/cclib
python -m pip install -e '/home/awallace43/gits/qcengine[test]'
```

## Q-Chem environment handling

The harness does not source a shell script. It validates the environment inherited by the Python process.

Before Q-Chem execution or version probing, it checks the Q-Chem 5.1 resources needed by the wrapper, including:

- `QC` identifies a readable Q-Chem installation directory.
- `QCAUX` identifies a readable auxiliary-data directory.
- `QCPROG` identifies an executable/readable program driver.
- The resolved `qchem` executable is runnable.

QCEngine supplies a managed scratch directory and sets `QCSCRATCH` for the child process. An inherited `QCSCRATCH` value is not required to equal the QCEngine scratch path.

Preflight errors identify every invalid variable without naming a machine-specific setup script. Example:

```text
Q-Chem environment variable QCAUX is not set or does not identify a readable directory. Initialize the Q-Chem environment before running cclib-qchem.
```

If Q-Chem output reports an undefined environment variable not caught by preflight, the harness extracts and names the actual variable when possible. Environment and license failures are classified as `ResourceError`, not parser or unknown failures.

## Execution and output selection

QCEngine's execution utilities own scratch creation, input writing, timeout handling, environment passing, stdout/stderr capture, and cleanup.

### Q-Chem

- Input file: `dispatch.in`
- Primary output: `dispatch.out`
- Success requires a successful process and Q-Chem's normal termination marker.
- `dispatch.out` is passed to cclib.

### ORCA

- Input file: `dispatch.inp`
- Primary output: captured stdout
- Success requires a successful process and `ORCA TERMINATED NORMALLY`.
- Captured stdout is passed to cclib through a disposable seekable text stream or scratch output file.

The harness retains a bounded diagnostic tail before schema conversion so failure messages remain useful even when parsing cannot begin.

## cclib parsing and QCSchema conversion

The shared conversion pipeline is:

1. Call cclib's auto-detecting parse API on the complete text output.
2. Reject `None` or an incomplete parser result.
3. Verify the selected parser class is `QChem` for `cclib-qchem` or `ORCA` for `cclib-orca`.
4. Require `ccData.metadata["success"] is True`.
5. Build a dictionary with `QCSchemaWriter(ccdata).as_dict(validate=False)`.
6. Add captured stdout, native input/output, and harness metadata without replacing cclib's flat extras.
7. Construct a QCElemental QCSchema v1 `AtomicResult` to validate the writer output.
8. Convert that result to QCSchema v2 for the internal QCEngine harness contract.
9. Set v2 `input_data` to the original requested `AtomicInput`; retain the parsed program-native result molecule as the output molecule.
10. Compare normalized parsed driver, method, and basis to the requested calculation.

Requested-versus-parsed comparison is case-insensitive and uses only the documented method aliases. Basis comparison is case-insensitive after trimming whitespace. A material mismatch raises `UnknownError`; the harness never relabels parsed data to make it appear consistent.

ORCA's cclib parser already incorporates its documented dispersion-energy adjustment into parsed SCF energies. The harness must not add dispersion a second time.

## Result data and provenance

The result preserves `QCSchemaWriter`'s flat cclib extras, including available atom charges, coordinates, atomic numbers, orbital energies and symmetries, SCF histories, moments, MP energies, counts, and cclib unit annotations.

The harness adds one non-colliding metadata object:

```python
extras["cclib_harness"] = {
    "selector": "cclib-qchem",  # or cclib-orca
    "cclib_version": "...",
    "parser": "QChem",          # or ORCA
    "executable": "...resolved path...",
}
```

Provenance semantics are:

- `creator`: actual parsed QC program (`QChem` or `ORCA`)
- `version`: program version parsed from the calculation output
- `routine`: cclib QCSchema writer routine
- QCEngine wrapper fields: added by QCEngine's normal metadata handling

`stdout` contains the complete primary program output. `native_files` contains the generated native input and primary native output when accepted by the active QCElemental model. If model-version constraints prevent a native file field, the generated input is retained under `extras["cclib_harness"]` metadata rather than discarded.

Runtime-dependent stdout, wrapper timing, hostname, executable path, program patch version, signed zero, and orientation details are not expected to match historical JSON byte-for-byte.

## Error model

### `InputError`

Raise before execution for:

- unsupported driver or method
- absent or structured basis
- ghost atoms
- malformed Q-Chem keyword values
- reserved Q-Chem key overrides
- malformed ORCA `simple` or `blocks` values
- reserved ORCA resource blocks
- unknown top-level ORCA keyword keys

### `ResourceError`

Raise for:

- missing or incompatible cclib
- executable absent from `PATH`
- wrong executable identity
- unsupported executable version
- missing or invalid program environment resource
- Q-Chem license/environment startup failure

### `UnknownError`

Raise for:

- nonzero program exit not classified as a resource/input failure
- absent normal-termination marker
- cclib auto-detection failure
- parser-class mismatch
- cclib unsuccessful metadata
- cclib parser exception
- `QCSchemaWriter` unsupported or incomplete result
- QCElemental result validation failure
- requested-versus-parsed driver, method, or basis mismatch

Every post-execution error includes selector, executable path, failure stage, and a bounded output tail. Exception chaining preserves the original parser or validation exception.

## Testing strategy

### Registration and availability

Tests verify:

- both selectors are registered and case-insensitive
- native `qchem` remains registered separately
- each selector has independent availability
- cclib import and writer feature checks
- executable absence
- wrong executable identity, including an unrelated `orca`
- version-floor rejection
- version cache separation by resolved path

These tests monkeypatch discovery and probes and require no proprietary executable.

### Input generation

Tests compare generated content and commands for Q-Chem and ORCA across:

- energy, gradient, and Hessian drivers
- HF, B3LYP, BP86, MP2, and CCSD
- closed- and open-shell molecules where supported
- core and memory mapping
- bohr/Angstrom coordinate handling
- harvesting defaults
- deterministic keyword ordering
- program-specific user keywords
- every pre-execution validation failure

Representative expectations derive from the approved cclib Q-Chem 5.1 and ORCA 6.0 input directories. Tests compare semantic sections and required directives rather than comments, whitespace, or original Z-matrix syntax.

### Execution and conversion

Mocked execution tests verify output selection, termination checks, parser-class enforcement, original input preservation, parsed output molecule retention, flat extras preservation, added harness metadata, stdout/native files, provenance, and schema v1/v2 conversion.

Real cclib integration tests parse representative Q-Chem 5.1 and ORCA 6 fixture outputs when `CCLIB_SOURCE_ROOT` names a cclib source checkout containing `data/`. For local validation, set `CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib`. These tests are optional integration tests rather than a duplicated cclib parser suite. The QCEngine unit suite does not contain absolute fixture paths.

At minimum, fixture integration covers:

- Q-Chem water MP2 energy
- Q-Chem B3LYP or BP86 energy
- Q-Chem B3LYP Hessian
- Q-Chem water CCSD energy
- ORCA water MP2 energy
- ORCA HF or B3LYP energy
- ORCA B3LYP Hessian
- ORCA water CCSD energy

### Failure tests

Tests cover:

- missing executable and cclib
- wrong executable identity
- each invalid Q-Chem environment variable
- generic extraction of a missing variable named by Q-Chem output
- license failure classification
- process failure and abnormal termination
- parser detection failure and parser mismatch
- unsuccessful cclib metadata
- unsupported writer method
- incomplete cclib attributes
- driver, method, and basis mismatch

### Live tests

Add `cclib-qchem` and `cclib-orca` pytest markers and availability entries in `qcengine/testing.py`. Live tests skip unless the fully configured executable and cclib are available.

On the target machine, live tests use the exact installations selected by the prepared `PATH`:

- Q-Chem 5.1 from `/projects/cos-lab-cs207/common/software/qchem5.1/bin/qchem`
- ORCA 6.1.1 from `/projects/cos-lab-cs207/common/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca`

A live water MP2/STO-3G smoke calculation is required for each selector.

## Required Q-Chem demonstration

Create `/home/awallace43/gits/qcengine/qchem_water_mp2.py`.

The script constructs the QCSchema equivalent of:

```text
/home/awallace43/gits/cclib/data/QChem/basicQChem5.1/water_mp2.in
```

It uses neutral singlet water, driver `energy`, method `mp2`, basis `sto-3g`, and the historical Cartesian geometry in bohr:

```python
symbols = ["O", "H", "H"]
geometry = [
    -0.0, 0.0, 0.22517858316070177,
    -1.4941103633283772, -0.0, -0.9007143324538345,
    1.4941103633283772, -0.0, -0.9007143324538345,
]
```

It calls:

```python
qcengine.compute(
    atomic_input,
    "cclib-qchem",
    raise_error=True,
    return_version=1,
)
```

The script does not read the historical `.out` file. It prints a concise pass/fail comparison and writes the complete result to `qchem_water_mp2.result.json`.

### Numerical acceptance

Using absolute tolerances unless stated otherwise:

| Field | Expected | Tolerance |
| --- | ---: | ---: |
| `return_result` | `-75.00228214` | `1e-6` Eh |
| `properties.return_energy` | `-75.00228214` | `1e-6` Eh |
| `properties.scf_total_energy` | `-74.9643287618` | `1e-6` Eh |
| `properties.mp2_total_energy` | `-75.00228214` | `1e-6` Eh |
| `properties.mp2_correlation_energy` | `-0.0379533782` | `1e-6` Eh |
| `properties.scf_dipole_moment` | `[0, 0, -0.6584056190]` | `1e-5` e·bohr per component |
| `calcinfo_nbasis` | `7` | exact |
| `calcinfo_nmo` | `7` | exact |
| `calcinfo_nalpha` | `5` | exact |
| `calcinfo_nbeta` | `5` | exact |
| `calcinfo_natom` | `3` | exact |
| `scf_iterations` | `6` | exact |

Program-build differences within these tolerances are accepted.

### cclib-rich acceptance

The demonstration also requires:

- `success is True`
- `schema_name == "qcschema_output"`
- `schema_version == 1`
- `driver == "energy"`
- model method/basis `mp2`/`sto-3g`
- neutral singlet O/H/H molecule
- provenance creator `QChem`
- provenance routine identifying cclib's QCSchema writer
- `extras["cclib_harness"]["selector"] == "cclib-qchem"`
- representative flat extras including `atomcharges`, `atomcoords`, `atomnos`, `homos`, `moenergies`, `mosyms`, `mpenergies`, `scfenergies`, `scftargets`, and `scfvalues`
- Mulliken charges approximately `[-0.339215, 0.169607, 0.169607]` within `1e-5`
- seven MO energies matching the historical output within `1e-3` Hartree
- six SCF history rows

Exact provenance version, runtime metadata, stdout text, executable path, signed zeros, and JSON formatting are not compared.

## Documentation

Document:

- selectors `cclib-qchem` and `cclib-orca`
- continued availability of native `qchem`
- supported drivers, methods, basis, and molecule limitations
- Q-Chem flat keyword and ORCA structured keyword conventions
- optional cclib development dependency
- PATH-only executable discovery
- requirement to initialize proprietary program environments before Python starts
- schema conversion and cclib extras behavior
- local live-verification commands as examples, clearly separated from portable harness behavior

Add a changelog entry describing the new optional cclib-backed selectors.

## Planned file surface

- Create `qcengine/programs/cclib.py`
- Modify `qcengine/programs/base.py`
- Modify `qcengine/testing.py`
- Modify `pyproject.toml` for pytest markers
- Create `qcengine/programs/tests/test_cclib.py`
- Create `qchem_water_mp2.py`
- Update the relevant program/developer documentation and changelog

No cclib file is modified by this implementation.

## Acceptance criteria

The implementation is complete only when all of the following hold:

1. `cclib-qchem` and `cclib-orca` are registered without replacing `qchem`.
2. Both instances are the same `CCLibHarness` class configured by `program`.
3. Availability requires cclib plus the correctly identified and sufficiently new executable from `PATH`.
4. The unrelated `/usr/bin/orca` does not make `cclib-orca` available.
5. Generated inputs satisfy the approved driver/method/basis/resource contracts and reproduce the semantics of representative cclib examples.
6. Executed output is parsed exclusively through cclib and converted through `QCSchemaWriter`.
7. Results preserve flat cclib extras and truthful program/cclib provenance.
8. Q-Chem environment failures identify the actual missing or invalid variable without machine-specific instructions.
9. Unit and optional fixture-integration tests pass.
10. The live Q-Chem demonstration meets every structural, numerical, and cclib-rich criterion above.
11. A live ORCA water MP2 smoke calculation confirms the second configured instance uses the shared harness.
12. No cclib source change is required; any newly discovered cclib blocker is reported before implementation scope changes.
