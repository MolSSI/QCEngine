# Feature: Gaussian Harness — End-to-End Coverage Tests <!-- rq-a9ab12b3 -->

This document specifies seven standalone end-to-end tests for the Gaussian harness,
all added to `qcengine/programs/tests/test_gaussian.py` and gated behind
`@using("gaussian")`. Each test exercises a code path that currently has no
real-Gaussian coverage (mocked unit tests aside) and that would silently regress
if cclib output formats, harness regex parsing, or method/driver routing changed.

These are not parametrisations of existing cross-program tests (those live in
`integration-tests.md`); they are additional standalone tests that fill gaps in
the per-harness coverage matrix.

Related files:
- `harvester.md` — harvester behaviour exercised by these tests
- `runner.md` — runner method/driver routing
- `dispersion.md` — dispersion-handling tests (the Hessian and gradient extensions of
  the existing `test_gaussian_b3lyp_d3_water` energy test)
- `integration-tests.md` — cross-program integration tests (out of scope here)

## Coverage Gaps Closed by This Feature <!-- rq-7472fc74 -->

The harness's end-to-end coverage matrix (method × driver) currently has the
following holes that this feature fills:

| Method | Energy | Gradient | Hessian | Properties |
|--------|--------|----------|---------|------------|
| HF     | ✓ | ✓ | **NEW (T1)** | **NEW (T2)** |
| Plain DFT (no dispersion) | **NEW (T3)** | — | — | — |
| DFT-D | ✓ | **NEW (T4)** | **NEW (T7)** | — |
| HF, single atom | **NEW (T6)** | — | — | — |
| UHF, open-shell radical | **NEW (T5)** | — | — | — |

(T1–T7 below.)

## Feature API <!-- rq-86480a7b -->

No new public functions or classes. Seven new pytest test functions in
`qcengine/programs/tests/test_gaussian.py`, each decorated with
`@using("gaussian")`, named:

- `test_gaussian_hessian_hf_water` (T1) <!-- rq-54afb7b8 -->
- `test_gaussian_properties_hf_water` (T2) <!-- rq-70b16c33 -->
- `test_gaussian_dft_b3lyp_water` (T3) <!-- rq-3727e53b -->
- `test_gaussian_b3lyp_d3_gradient_water` (T4) <!-- rq-66ca9a20 -->
- `test_gaussian_uhf_nh2` (T5) <!-- rq-d06d7266 -->
- `test_gaussian_single_atom_he` (T6) <!-- rq-cd65e68a -->
- `test_gaussian_b3lyp_d3_hessian_water` (T7) <!-- rq-b8065553 -->

## Reference Values <!-- rq-7c70da2a -->

All reference values below were captured against Gaussian 16 Revision C.02 with
the existing `water` fixture (charge 0, multiplicity 1; standard equilibrium
geometry — see `_make_water_mol` in the test file) unless otherwise noted.
Tolerances are stated per test.

### T1: HF Hessian on water <!-- rq-6385e33d -->

Method `HF`, basis `STO-3G`, driver `hessian`. Reference values (Hartree/Bohr²):

- `return_result.shape == (9, 9)`
- `numpy.trace(hess) ≈ 3.2315247000` (atol `1e-6`)
- `numpy.linalg.norm(hess) ≈ 2.2444460940` (atol `1e-6`)
- `hess[0, 0] ≈ -0.0560541000` (atol `1e-6`)
- `hess[3, 3] ≈ -0.0223351000` (atol `1e-6`)
- `hess[0, 1] ≈ 0.0` (atol `1e-10`)
- Matrix is symmetric within `1e-10`: `numpy.allclose(hess, hess.T, atol=1e-10)`

Hessian unit tests use a synthetic identity matrix; this test exercises the real
`_parse_force_constants` regex parse, the D-exponent → float conversion via
`_fortran_float`, the lower-triangle reconstruction, and the
`mill.align_hessian` frame transformation, none of which are otherwise exercised
end-to-end.

### T2: HF properties on water <!-- rq-e199db9a -->

Method `HF`, basis `STO-3G`, driver `properties`. Reference values:

- `result.return_result` is a float (the total energy)
- `qcvars["HF TOTAL ENERGY"] ≈ -74.962963201` (atol `1e-6`)
- `qcvars["DIPOLE MOMENT"] ≈ 1.7252` Debye (atol `1e-3`)
- `qcvars["MULLIKEN CHARGES"]` is a list of length 3
- `sum(qcvars["MULLIKEN CHARGES"]) ≈ 0.0` (atol `1e-4`)
- `qcvars["MULLIKEN CHARGES"][0] ≈ -0.366114` (oxygen, atol `1e-3`)

Exercises the `properties` driver branch (Pop=Full route emission) and the
harvester's `ccdata.moments[1]` and `ccdata.atomcharges["mulliken"]` extraction.

### T3: Plain B3LYP energy on water <!-- rq-1e76932b -->

Method `B3LYP`, basis `cc-pVDZ`, driver `energy`, no dispersion suffix. Reference values:

- `result.return_result ≈ -76.4203540604` (atol `1e-6`)
- `qcvars["HF TOTAL ENERGY"] ≈ -76.4203540604` (atol `1e-6` — the SCF energy IS the B3LYP energy when no dispersion)
- `qcvars["CURRENT ENERGY"] ≈ -76.4203540604` (atol `1e-6`)
- `"DISPERSION CORRECTION ENERGY"` is **not** present in qcvars
- `"B3LYP-D3 DISPERSION CORRECTION ENERGY"` is **not** present in qcvars

Confirms that the unknown-functional pass-through branch in `germinate.py`
emits a usable route line, that Gaussian recognises `B3LYP` as a functional,
and that the harvester does not emit spurious dispersion qcvars when no
dispersion suffix was supplied.

### T4: B3LYP-D3 gradient on water <!-- rq-80a8de50 -->

Method `b3lyp-d3`, basis `6-31g*`, driver `gradient`. Reference values
(Hartree/Bohr):

- `result.return_result.shape == (3, 3)` (natoms × 3)
- `result.return_result ≈`
  ```
  [[ 0.0,         0.0,        -0.01437422],
   [ 0.0,        -0.00800138,  0.00718711],
   [ 0.0,         0.00800138,  0.00718711]]
  ```
  with `atol = 1e-6`.
- `qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"] ≈ -0.0000078919` (atol `1e-6`)
- `qcvars["B3LYP FUNCTIONAL TOTAL ENERGY"] ≈ -76.4087263383` (atol `1e-6`)
- `qcvars["B3LYP-D3 TOTAL ENERGY"] ≈ -76.4087342302` (atol `1e-6`)

Confirms that `EmpiricalDispersion=GD3` combined with `Force=NoStep` (the
Gaussian gradient request emitted by `muster_modelchem`) produces a coherent
gradient and that cclib's `grads` attribute already contains the
dispersion-corrected total gradient.

### T5: NH₂· UHF on cc-pVDZ <!-- rq-c38f1b61 -->

Open-shell radical. Geometry (Å):
```
0 2
N  0.000000000  0.000000000  0.140030000
H  0.000000000  0.802300000 -0.490100000
H  0.000000000 -0.802300000 -0.490100000
```
(r(N–H) ≈ 1.020 Å, ∠HNH ≈ 103.7°.)

Method `HF`, basis `cc-pVDZ`, driver `energy`, charge 0, multiplicity 2.
Reference values:

- `result.return_result ≈ -55.5671259096` (atol `1e-6`)
- `qcvars["HF TOTAL ENERGY"] ≈ -55.5671259096` (atol `1e-6`)
- `native_files["input"]` (the .com file) contains the substring `"UHF/cc-pVDZ"`
  (the U-prefix was auto-prepended by `muster_modelchem` for multiplicity > 1)
- `native_files["input"]` contains the line `"0 2"` (charge/multiplicity)

Exercises the U-prefix logic on a real open-shell calculation (the existing
`test_hf_alignment` covers UHF but does not assert the .com-file U-prefix).

### T6: Single-atom He HF/cc-pVDZ <!-- rq-954dd6f3 -->

Method `HF`, basis `cc-pVDZ`, driver `energy`. Molecule: single He atom at
origin, charge 0, multiplicity 1. Reference values:

- `result.return_result ≈ -2.85516047724` (atol `1e-6`)
- `result.properties.nuclear_repulsion_energy == 0.0` (exact)
- `result.success is True`

Edge case: validates that a single-atom system (no NRE, no bonds, special
Gaussian handling for atoms) flows through the harness without geometry- or
NRE-related off-by-one errors.

### T7: B3LYP-D3 Hessian on water <!-- rq-dbf8a9a5 -->

Method `b3lyp-d3`, basis `6-31g*`, driver `hessian`. Reference values:

- `return_result.shape == (9, 9)`
- `numpy.trace(hess) ≈ 2.3626060600` (atol `1e-6`)
- `numpy.linalg.norm(hess) ≈ 1.5979941261` (atol `1e-6`)
- `qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"] ≈ -0.0000078919` (atol `1e-6`)
- Matrix is symmetric: `numpy.allclose(hess, hess.T, atol=1e-10)`

Cheap follow-on to T1 and T4: confirms the Hessian + dispersion coexistence
case (Gaussian's analytic Hessian with `EmpiricalDispersion=GD3` succeeds and
the harness handles the combined output).

## Gherkin Scenarios <!-- rq-051e6b78 -->

```gherkin
Feature: Gaussian harness end-to-end coverage tests

  Background:
    Given the Gaussian executable (g16 or g09) is available on PATH
    And cclib is importable
    And the standard `water` fixture is used unless otherwise noted

  # --- T1: HF Hessian ---

  @rq-89896430
  Scenario: HF/STO-3G Hessian on water has correct shape and reference values
    Given an AtomicInput for water with method "HF", basis "STO-3G", driver "hessian"
    When qcng.compute is called with program "gaussian"
    Then the result is successful
    And result.return_result is a 2D array of shape (9, 9)
    And numpy.trace(return_result) is approximately 3.2315247 (atol 1e-6)
    And numpy.linalg.norm(return_result) is approximately 2.244446094 (atol 1e-6)
    And return_result is symmetric to within 1e-10
    And return_result[0, 1] is approximately 0.0 (atol 1e-10)

  # --- T2: HF properties ---

  @rq-dc1bf912
  Scenario: HF/STO-3G properties on water populates dipole and Mulliken charges
    Given an AtomicInput for water with method "HF", basis "STO-3G", driver "properties"
    When qcng.compute is called with program "gaussian"
    Then the result is successful
    And qcvars["HF TOTAL ENERGY"] is approximately -74.962963201 (atol 1e-6)
    And qcvars["DIPOLE MOMENT"] is approximately 1.7252 Debye (atol 1e-3)
    And qcvars["MULLIKEN CHARGES"] is a list of length 3
    And sum(qcvars["MULLIKEN CHARGES"]) is approximately 0.0 (atol 1e-4)
    And qcvars["MULLIKEN CHARGES"][0] is approximately -0.366114 (atol 1e-3)

  # --- T3: Plain DFT (no dispersion) ---

  @rq-af1fd6a9
  Scenario: B3LYP/cc-pVDZ energy on water matches reference without emitting dispersion qcvars
    Given an AtomicInput for water with method "B3LYP", basis "cc-pVDZ", driver "energy"
    When qcng.compute is called with program "gaussian"
    Then the result is successful
    And result.return_result is approximately -76.4203540604 (atol 1e-6)
    And qcvars["HF TOTAL ENERGY"] equals result.return_result (within 1e-10)
    And qcvars["CURRENT ENERGY"] equals result.return_result (within 1e-10)
    And "DISPERSION CORRECTION ENERGY" is not present in qcvars
    And "B3LYP-D3 DISPERSION CORRECTION ENERGY" is not present in qcvars

  # --- T4: DFT-D gradient ---

  @rq-b10501a0
  Scenario: B3LYP-D3/6-31G* gradient on water matches reference and exposes dispersion qcvars
    Given an AtomicInput for water with method "b3lyp-d3", basis "6-31g*", driver "gradient"
    When qcng.compute is called with program "gaussian"
    Then the result is successful
    And result.return_result is a 2D array of shape (3, 3)
    And result.return_result matches the reference gradient (atol 1e-6)
    And qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"] is approximately -0.0000078919 (atol 1e-6)
    And qcvars["B3LYP FUNCTIONAL TOTAL ENERGY"] is approximately -76.4087263383 (atol 1e-6)
    And qcvars["B3LYP-D3 TOTAL ENERGY"] is approximately -76.4087342302 (atol 1e-6)

  # --- T5: UHF open-shell ---

  @rq-18a8cc19
  Scenario: UHF/cc-pVDZ energy on NH2 radical matches reference and emits UHF route line
    Given an AtomicInput for NH2 radical (charge=0, multiplicity=2) with method "HF", basis "cc-pVDZ", driver "energy"
    When qcng.compute is called with program "gaussian"
    Then the result is successful
    And result.return_result is approximately -55.5671259096 (atol 1e-6)
    And the input .com file (native_files["input"]) contains "UHF/cc-pVDZ"
    And the input .com file contains a charge/multiplicity line "0 2"

  # --- T6: Single-atom edge case ---

  @rq-5c04235d
  Scenario: HF/cc-pVDZ energy on a single He atom returns the expected value with NRE=0
    Given an AtomicInput for a single He atom at origin with method "HF", basis "cc-pVDZ", driver "energy"
    When qcng.compute is called with program "gaussian"
    Then the result is successful
    And result.return_result is approximately -2.85516047724 (atol 1e-6)
    And result.properties.nuclear_repulsion_energy equals 0.0 (exactly)

  # --- T7: DFT-D Hessian ---

  @rq-cedda71e
  Scenario: B3LYP-D3/6-31G* Hessian on water succeeds and emits the dispersion qcvar
    Given an AtomicInput for water with method "b3lyp-d3", basis "6-31g*", driver "hessian"
    When qcng.compute is called with program "gaussian"
    Then the result is successful
    And result.return_result is a 2D array of shape (9, 9)
    And numpy.trace(return_result) is approximately 2.3626060600 (atol 1e-6)
    And numpy.linalg.norm(return_result) is approximately 1.5979941261 (atol 1e-6)
    And return_result is symmetric to within 1e-10
    And qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"] is approximately -0.0000078919 (atol 1e-6)
```

## Implementation Notes <!-- rq-1ed899e3 -->

### Fixtures and helpers <!-- rq-c88cc8ec -->

- The existing `water` fixture (`_make_water_mol`) is reused for T1, T2, T3, T4, T7.
- A new `nh2_radical` fixture (or inline `Molecule.from_data(...)`) is added for T5.
- A new `helium_atom` fixture (or inline) is added for T6.
- All seven tests use the standard `@using("gaussian")` decorator and the existing
  `qcng.compute(..., raise_error=True)` pattern.

### Reference-value provenance <!-- rq-ce2ecf6a -->

All reference values were captured against **Gaussian 16 Revision C.02** during
planning. Minor patch-level Gaussian updates may shift values within the stated
tolerance (`atol = 1e-6` is generous enough to absorb sub-µHartree changes in
SCF, dispersion, and integration grids). If a future Gaussian release exceeds
this tolerance, the test should be updated rather than the tolerance loosened —
unexpected drift in reference values is itself worth investigating.

### Test runtime budget <!-- rq-d42168a0 -->

Estimated wall time for the full suite of 7 E2E tests on a single core,
Gaussian 16 C.02:

| Test | Approx. wall time |
|------|-------------------|
| T1 HF Hessian / STO-3G water | ~1 s |
| T2 HF properties / STO-3G water | <1 s |
| T3 B3LYP / cc-pVDZ water | ~2 s |
| T4 B3LYP-D3 gradient / 6-31G* water | ~3 s |
| T5 UHF / cc-pVDZ NH₂ | ~2 s |
| T6 HF / cc-pVDZ He | <1 s |
| T7 B3LYP-D3 Hessian / 6-31G* water | ~5 s |
| **Total** | **~15 s** |

This is comparable to the existing Gaussian integration test runtime
(test_simple_ghost + test_tricky_ghost + test_atom_labels currently take
~60 s) and fits well within CI budgets.

