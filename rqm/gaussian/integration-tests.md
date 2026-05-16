# Feature: Gaussian Harness --- Cross-Program Integration Tests <!-- rq-0f585dce -->

This document specifies the integration of the Gaussian harness into the existing
cross-program test suites in `qcengine/programs/tests/`. These tests validate that
Gaussian conforms to the same interfaces and conventions as other harnesses (CFOUR,
GAMESS, NWChem, Psi4).

Ghost atom tests (`test_ghost.py`) are out of scope for this document.

Related requirements files:
- `runner.md` --- GaussianHarness class and `_defaults`
- `input-builder.md` --- `.com` file construction
- `harvester.md` --- cclib-based output parsing

## Prerequisite: native_files Key Convention <!-- rq-bb7d680f -->

The existing `parse_output` method in `runner.py` constructs `native_files` with
program-specific keys (`"gaussian.com"`, `"gaussian.log"`). The cross-program
`test_protocol_native` test expects:

- An `"input"` key for the input file (matching CFOUR, GAMESS, NWChem, Psi4).
- Additional output files keyed by their natural names (e.g. `"gaussian.log"`).
- The `outfiles` dict from `execute()` already contains the `"input"` key (set in
  `compute()` at line 160 of `runner.py`).

**Required change**: `parse_output` must switch from a hardcoded `native_files` dict to
the `{k: v for k, v in outfiles.items() if v is not None}` pattern used by every other
harness. After popping `"stdout"` and `"stderr"`, the remaining keys (`"input"`,
`"gaussian.log"`) become the native files. The `"gaussian.com"` key is dropped in
favour of the standard `"input"` key.

## Canonical Methods Entry <!-- rq-e71357c5 -->

Gaussian must be added to the `_canonical_methods` list in `test_canonical_config.py`:

```python
pytest.param("gaussian", {"method": "hf", "basis": "STO-3G"}, {}, marks=using("gaussian")),
```

This entry is shared by `test_local_options_memory_gib`, `test_local_options_scratch`,
`test_local_options_ncores`, and `test_protocol_native`.

The `_get_molecule` helper returns `"hydrogen"` by default, which is suitable for
Gaussian HF/STO-3G.

## Test Integrations <!-- rq-8070c4fd -->

### `test_local_options_memory_gib` (`test_canonical_config.py`) <!-- rq-0edebcd3 -->

Validates that memory configuration flows through to the program.

The Gaussian harness has `managed_memory=True`, so a `stdout_ref` regex is required.
Gaussian writes memory to its `.com` file as `%Mem=<N>GB` and this text appears in
`ret.stdout` via the log or input echo. However, Gaussian's primary evidence of memory
is in the input file (Link 0 section), not in stdout. The log file echoes the route but
not the Link 0 section.

**Required `stdout_ref` entry**: A regex that matches the memory setting in the
Gaussian log output. Gaussian echoes a line like
`%NProc...` and memory information in the log. If the log does not reliably echo
`%Mem`, the harness may need to append the input file content to `stdout` (as CFOUR
does with a "freebie" pattern), or the `stdout_ref` should match a line that Gaussian
always prints showing the memory available. The exact pattern must be determined
empirically by inspecting a real Gaussian log.

**Required `memory_trickery` entry** (for the `"dsl"` parametrization): Gaussian does
not have a native keyword that sets memory in its route line (memory is a Link 0
command, not a route keyword), so the `memory_trickery` dict for Gaussian should be
empty `{}`, matching the Psi4 pattern:
```python
"gaussian": {},
```

### `test_local_options_scratch` (`test_canonical_config.py`) <!-- rq-815d5fee -->

Validates that scratch directory configuration flows through to the program.

The Gaussian harness has `scratch=True`.

**Required `scratch_sample` entry**: A glob pattern for a file left behind in the
scratch directory when `scratch_messy=True`. The Gaussian log file is the natural
choice:
```python
"gaussian": "*/gaussian.log",
```

**Required `stdout_ref` entry**: A regex that matches scratch-related output. Gaussian
does not print its scratch directory path in the log. Use a "freebie" pattern (a string
that always appears in Gaussian output), matching the CFOUR/GAMESS/NWChem convention:
```python
"gaussian": "Gaussian, Inc.",
```

### `test_local_options_ncores` (`test_canonical_config.py`) <!-- rq-3a18d6af -->

Validates that thread/core configuration flows through to the program.

The Gaussian harness has `thread_parallel=True`, so a `stdout_ref` regex is required.

**Required `stdout_ref` entry**: Gaussian writes `%NProcShared=<N>` to the Link 0
section of the `.com` file and the log file header echoes `Will use up to <N>
processors` or a similar line. The exact pattern must be determined empirically. A
candidate regex:
```python
"gaussian": rf"NProcShared={ncores}",
```
This requires the input echo or log to contain the NProcShared value. If it does not
appear in stdout, an alternative freebie pattern may be needed, with a comment
explaining why.

### `test_protocol_native` (`test_canonical_fields.py`) <!-- rq-c01c892e -->

Validates the `native_files` protocol (`"none"`, `"input"`, `"all"`).

**Depends on** the native_files key convention fix described above.

**Required `input_ref` entry**: A regex matching content in the Gaussian input file
(stored under the `"input"` key). The `.com` file always contains the route line:
```python
"gaussian": rf"#\s*HF/STO-3G",
```

**Required `all_ref` entry**: A tuple `(filename, regex)` for an ancillary output file
present when `native == "all"`. The Gaussian log is the natural choice:
```python
"gaussian": ("gaussian.log", rf"Normal termination of Gaussian"),
```

### `test_hf_alignment` (`test_alignment.py`) <!-- rq-81bc1d68 -->

Validates that output molecule geometry is correctly aligned to the input molecule
across energy, gradient, and hessian drivers, with randomised input frame scrambling.
Tests both RHF (closed-shell) and UHF (open-shell) references.

**Required `inp` parametrizations**:

RHF:
```python
pytest.param(
    {"call": "gaussian", "reference": "rhf", "fcae": "ae", "keywords": {}},
    id="hf  rhf ae: gaussian",
    marks=using("gaussian"),
),
```

UHF:
```python
pytest.param(
    {"call": "gaussian", "reference": "uhf", "fcae": "ae", "keywords": {}},
    id="hf  uhf ae: gaussian",
    marks=using("gaussian"),
),
```

No special keywords are needed: the Gaussian harness auto-prepends `"U"` to the method
for open-shell systems (multiplicity > 1), and basis set names are passed through
directly.

**Reference data**: The standard suite reference file
(`standard_suite_ref.py`) must include Gaussian HF energy and gradient/hessian
reference values for the test molecules (`hf`, `bh3p`, `h2o`, `nh2`) at the
`cc-pVDZ` and `aug-cc-pVDZ` basis sets. These values must be determined empirically
by running Gaussian.

### `test_sp_ccsd_t_rhf_full` (`test_standard_suite_ccsd(t).py`) <!-- rq-9aa3cac4 -->

Validates CCSD(T) single-point energy on water (RHF reference) against a known value.

**Required parametrization**:
```python
pytest.param("gaussian", "aug-cc-pVDZ", {}, marks=using("gaussian")),
```

**Reference value**: Gaussian will produce a CCSD(T)/aug-cc-pVDZ energy for water that
may differ slightly from the values produced by other programs (the existing reference
is -76.276030676767 Hartree). A Gaussian-specific reference value must be determined
empirically by running the calculation, then hard-coded into the test with the same
`atol = 1.0e-6` tolerance.

## Gherkin Scenarios <!-- rq-e7c13c1e -->

```gherkin
Feature: Gaussian cross-program integration tests

  Background:
    Given the Gaussian executable (g16 or g09) is available on PATH
    And cclib is importable
    And GaussianHarness is registered in _canonical_methods as ("gaussian", {"method": "hf", "basis": "STO-3G"}, {})

  # --- Canonical methods entry ---

  @rq-647b596c
  Scenario: Gaussian appears in _canonical_methods
    When _canonical_methods is inspected
    Then it contains a pytest.param for "gaussian" with model {"method": "hf", "basis": "STO-3G"}

  # --- native_files key convention ---

  @rq-f13e0d24
  Scenario: parse_output uses standard "input" key in native_files
    Given a successful HF/STO-3G energy calculation on hydrogen
    When GaussianHarness.compute() returns an AtomicResult
    Then native_files contains the key "input" (not "gaussian.com")
    And native_files["input"] contains the .com file content
    And native_files contains the key "gaussian.log"
    And native_files does not contain the keys "stdout" or "stderr"

  # --- Memory test ---

  @rq-65cb367e
  Scenario: Memory configuration is applied with 1.555 GiB total node memory
    Given a TaskConfig with memory = 1.555 GiB and ncores = 1
    When a HF/STO-3G energy calculation on hydrogen is run through qcng.compute
    Then the result is successful
    And the Gaussian input or output contains evidence of the memory setting

  @rq-2aa0e9e1
  Scenario: Memory configuration is not overridden by conflicting DSL keywords
    Given a TaskConfig with memory = 1.555 GiB
    And no conflicting memory keyword is provided (Gaussian memory is Link 0, not a route keyword)
    When a HF/STO-3G energy calculation on hydrogen is run
    Then the result is successful
    And the memory setting matches the TaskConfig value

  # --- Scratch test ---

  @rq-71824b99
  Scenario: Scratch directory is used and files are left behind when scratch_messy=True
    Given a TaskConfig with a custom scratch_directory and scratch_messy = True
    When a HF/STO-3G energy calculation on hydrogen is run through qcng.compute
    Then the result is successful
    And the file "*/gaussian.log" exists in the scratch directory

  # --- Ncores test ---

  @rq-519fbe49
  Scenario: Single-core execution succeeds
    Given a TaskConfig with ncores = 1
    When a HF/STO-3G energy calculation on hydrogen is run through qcng.compute
    Then the result is successful

  @rq-4252ac2d
  Scenario: Multi-core execution succeeds and ncores is reflected in output
    Given a TaskConfig with ncores = 3
    When a HF/STO-3G energy calculation on hydrogen is run through qcng.compute
    Then the result is successful
    And the Gaussian output contains evidence that 3 processors were requested

  # --- native_files protocol ---

  @rq-8df082e3
  Scenario: native_files protocol "none" returns empty native_files
    Given an AtomicInput with protocols.native_files = "none"
    When a HF/STO-3G gradient calculation on hydrogen is run through qcng.compute
    Then the result is successful
    And native_files is empty

  @rq-923b3212
  Scenario: native_files protocol "input" returns only the input file
    Given an AtomicInput with protocols.native_files = "input"
    When a HF/STO-3G gradient calculation on hydrogen is run through qcng.compute
    Then the result is successful
    And native_files contains exactly the key "input"
    And native_files["input"] matches the regex "#\s*HF/STO-3G"

  @rq-aa991a3a
  Scenario: native_files protocol "all" returns input and output files
    Given an AtomicInput with protocols.native_files = "all"
    When a HF/STO-3G gradient calculation on hydrogen is run through qcng.compute
    Then the result is successful
    And native_files contains the key "input"
    And native_files contains the key "gaussian.log"
    And native_files["gaussian.log"] matches the regex "Normal termination of Gaussian"

  @rq-d87ddeb8
  Scenario: native_files never contains "stdout" or "stderr"
    Given an AtomicInput with protocols.native_files = "all"
    When a HF/STO-3G gradient calculation on hydrogen is run through qcng.compute
    Then native_files does not contain "stdout"
    And native_files does not contain "stderr"
    And extras does not contain "outfiles"

  # --- HF Alignment (RHF) ---

  @rq-ba0f397f
  Scenario: RHF energy is invariant to input frame scrambling
    Given a closed-shell molecule (e.g. h2o) with randomly scrambled Cartesian frame
    And the input uses method "hf" with a cc-pVDZ or aug-cc-pVDZ basis
    And keywords = {}
    When qcng.compute is called with program "gaussian" and driver "energy"
    Then the result is successful
    And the return_result matches the reference energy within tolerance

  @rq-637d486d
  Scenario: RHF gradient is correctly aligned to the (scrambled) input frame
    Given a closed-shell molecule with randomly scrambled Cartesian frame and fix_com=True, fix_orientation=True
    And the input uses method "hf" with a cc-pVDZ basis
    When qcng.compute is called with program "gaussian" and driver "gradient"
    Then the result is successful
    And the return_result gradient matches the reference gradient (transformed to the scrambled frame) within tolerance

  @rq-87c8fb49
  Scenario: RHF hessian is correctly aligned to the (scrambled) input frame
    Given a closed-shell molecule with randomly scrambled Cartesian frame and fix_com=True, fix_orientation=True
    And the input uses method "hf" with a cc-pVDZ basis
    When qcng.compute is called with program "gaussian" and driver "hessian"
    Then the result is successful
    And the return_result hessian matches the reference hessian (transformed to the scrambled frame) within tolerance

  @rq-6a889dea
  Scenario: UHF energy is invariant to input frame scrambling
    Given an open-shell molecule (e.g. nh2, multiplicity=2) with randomly scrambled Cartesian frame
    And the input uses method "hf" with a cc-pVDZ or aug-cc-pVDZ basis
    And keywords = {}
    When qcng.compute is called with program "gaussian" and driver "energy"
    Then the result is successful
    And the return_result matches the reference UHF energy within tolerance

  @rq-349bf8c3
  Scenario: UHF gradient is correctly aligned to the (scrambled) input frame
    Given an open-shell molecule with randomly scrambled Cartesian frame and fix_com=True, fix_orientation=True
    And the input uses method "hf" with a cc-pVDZ basis
    When qcng.compute is called with program "gaussian" and driver "gradient"
    Then the result is successful
    And the return_result gradient matches the reference UHF gradient within tolerance

  @rq-5d9b0b40
  Scenario: UHF hessian is correctly aligned to the (scrambled) input frame
    Given an open-shell molecule with randomly scrambled Cartesian frame and fix_com=True, fix_orientation=True
    And the input uses method "hf" with a cc-pVDZ basis
    When qcng.compute is called with program "gaussian" and driver "hessian"
    Then the result is successful
    And the return_result hessian matches the reference UHF hessian within tolerance

  # --- CCSD(T) standard suite ---

  @rq-475c8796
  Scenario: CCSD(T)/aug-cc-pVDZ RHF energy on water matches Gaussian-specific reference
    Given the water molecule (H2O, charge=0, multiplicity=1) in Bohr coordinates
    And the model is {"method": "ccsd(t)", "basis": "aug-cc-pVDZ"}
    And keywords = {}
    When qcng.compute is called with program "gaussian" and driver "energy"
    Then the result is successful
    And provenance is present
    And return_result matches the Gaussian-specific CCSD(T) reference value within atol = 1.0e-6
```

## Implementation Notes <!-- rq-f28b3bed -->

### Empirical Reference Values <!-- rq-56e891e7 -->

The following reference values must be determined by running Gaussian before the
tests can be finalised:

1. **CCSD(T)/aug-cc-pVDZ water energy** --- for `test_sp_ccsd_t_rhf_full`.
2. **HF/{cc-pVDZ,aug-cc-pVDZ} energies/gradients/hessians** for `hf`, `bh3p`, `h2o`,
   `nh2` molecules --- for `test_hf_alignment` via `standard_suite_ref.py`.
3. **`stdout_ref` patterns** for memory and ncores tests --- determined by inspecting
   actual Gaussian log output for the relevant Link 0 echo lines.

### Files Modified <!-- rq-6423b98b -->

| File | Change |
|------|--------|
| `qcengine/programs/gaussian/runner.py` | Fix `parse_output` to use outfiles-passthrough pattern for `native_files` |
| `qcengine/programs/tests/test_canonical_config.py` | Add Gaussian to `_canonical_methods`; add `stdout_ref`, `scratch_sample`, `memory_trickery` entries |
| `qcengine/programs/tests/test_canonical_fields.py` | Add `input_ref` and `all_ref` entries for Gaussian |
| `qcengine/programs/tests/test_alignment.py` | Add RHF and UHF `inp` parametrizations for Gaussian |
| `qcengine/programs/tests/test_standard_suite_ccsd(t).py` | Add Gaussian parametrization with empirical reference value |
| `qcengine/programs/tests/standard_suite_ref.py` | Add Gaussian HF reference data for alignment test molecules |
