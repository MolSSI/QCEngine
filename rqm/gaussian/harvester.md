# Feature: Gaussian Harness — Harvester <!-- rq-9d12489e -->

This document specifies the output-parsing layer for the Gaussian harness. All parsing
lives in `qcengine/programs/gaussian/harvester.py`. The harvester uses **cclib** as the
primary parsing backend, supplemented by regex extraction for quantities that cclib does
not expose (e.g. normal/abnormal termination, version identification, error messages).

## cclib Integration <!-- rq-bb0a5d96 -->

cclib is called once per log file via:

```python
import io
import cclib

ccdata = cclib.io.ccopen(io.StringIO(log_content)).parse()
```

The `ccdata` object exposes parsed attributes used as follows:

| cclib attribute | Units (cclib) | Conversion | QCEngine usage |
|-----------------|---------------|------------|----------------|
| `ccdata.scfenergies[-1]` | eV | `cclib.parser.utils.convertor("eV"→"hartree")` | SCF/HF/DFT total energy |
| `ccdata.mpenergies[-1][-1]` | eV | `cclib.parser.utils.convertor("eV"→"hartree")` | MP2 total energy |
| `ccdata.ccenergies[-1]` | eV | `cclib.parser.utils.convertor("eV"→"hartree")` | CCSD or CCSD(T) total energy |
| `ccdata.dispersionenergies[-1]` | eV | `cclib.parser.utils.convertor("eV"→"hartree")` | Empirical dispersion correction (see [[dispersion]]) |
| `ccdata.grads[-1]` | Hartree/Bohr | negate (forces → gradient) | Gradient array |
| `ccdata.atomcoords[-1]` | Ångström | × Bohr/Å | Geometry for output molecule (standard orientation) |
| `ccdata.atomnos` | atomic numbers | via qcel periodic table | Element symbols (no-ghost branch only) |
| `ccdata.natom` | integer | none | Atom count for the sanity check and downstream sizing (preferred over `len(atomnos)`) |
| `ccdata.nbasis` | integer | none | Number of basis functions |
| `ccdata.nmo` | integer | none | Number of molecular orbitals |
| `ccdata.homos` | integer array | +1 for count | N alpha/beta electrons |
| `ccdata.moments[1]` | Debye | none | Dipole moment vector |
| `ccdata.atomcharges["mulliken"]` | e | none | Mulliken partial charges |

### cclib Energy Precision <!-- rq-f6537c40 -->

cclib stores parsed energies in electron volts. For Gaussian, the parser
reads the Hartree value from a line such as `SCF Done: E(RHF) =
-76.0571491640 A.U.`, calls `utils.float()` on the token (just `float()` with
D-exponent handling), multiplies by cclib's Hartree→eV factor (`27.21138505`),
and stores the result on `ccdata.scfenergies` (and analogously on
`mpenergies`, `ccenergies`, `dispersionenergies`).

The harvester reverses this with `cclib.parser.utils.convertor(x, "eV",
"hartree")`, which uses the same `27.21138505` factor in both directions
inside cclib's `_convertor` dict. The round-trip precision is:

| Step | Loss |
|------|------|
| Gaussian text → Python `float` via `utils.float()` | none (the text fits within IEEE-754 double precision) |
| `× cclib's Hartree→eV factor` | none beyond 1 ULP of float multiplication |
| `× cclib's eV→Hartree factor` (via `cclib.parser.utils.convertor`) | none beyond 1 ULP |
| **End-to-end Ha → eV → Ha round-trip** | **≤ 1 ULP (~1×10⁻¹⁴ Ha for ~100-Ha energies); often exactly 0** |

This is 7+ orders of magnitude below any tolerance the test suite cares
about (e.g. `atol = 2e-7` for cross-program alignment).

The reverse step **must** use `cclib.parser.utils.convertor` rather than
`qcel.constants.conversion_factor("eV", "hartree")`: qcel and cclib carry
different CODATA Hartree↔eV factors (qcel: `27.21138602`; cclib:
`27.21138505`), and reversing cclib's eV-stored value with qcel's factor
introduces a `~3.6×10⁻⁶ Ha` bias for typical SCF energies that exceeds the
cross-program tolerance.

All SCF, MP*, CC, and dispersion energies are read directly from cclib's
parsed attributes (`ccdata.scfenergies`, `ccdata.mpenergies`,
`ccdata.ccenergies`, `ccdata.dispersionenergies`) via this converter. No
log-text regex is used for energies — see "Regex Use" below for the items
that warrant regex.

### cclib Limitation: Full Cartesian Hessian not extracted <!-- rq-a3c37e5d -->

cclib's Gaussian parser does populate `ccdata.vibfconsts` (per-normal-mode reduced
force constants, shape `(nmodes,)`) from the `"Frc consts --"` lines of a `Freq`
job. It does **not**, however, populate `ccdata.hessian` (the full `3N × 3N`
Cartesian Hessian matrix) for Gaussian logs — that attribute is only assigned by
parsers for programs whose output prints the full matrix in a format cclib
recognises. The "Force constants in Cartesian coordinates" lower-triangular block
that Gaussian prints is not parsed by cclib.

**Resolution:** The full Cartesian Hessian is parsed directly from the log via
`_parse_force_constants()` (see "Hessian Parsing" below). `ccdata.vibfconsts` is
not used by this harvester — the full matrix is required for QCSchema's
`return_result` and downstream consumers.

### Coordinate Frame <!-- rq-f4900ee1 -->

cclib's `atomcoords[-1]` returns geometry in Gaussian's **standard orientation**
(the internally reoriented frame). This affects frame alignment logic — see
"Output Molecule Construction" below.

**Important**: Gaussian reports atomic **forces** (= −dE/dR) in its output; cclib
preserves this sign convention in `ccdata.grads`. The gradient (= +dE/dR) is obtained by
**negating** the forces array:

```python
gradient = -ccdata.grads[-1]   # shape (natoms, 3)
```

## Energy Extraction <!-- rq-c20746b1 -->

Energies are read from cclib's parsed attributes and converted eV → Hartree
via `cclib.parser.utils.convertor`. The per-method attribute mapping is:

| Method | cclib source | qcvars set |
|--------|--------------|------------|
| `"hf"` / `"scf"` / DFT (default) | `ccdata.scfenergies[-1]` | `HF TOTAL ENERGY`, `SCF TOTAL ENERGY`, `CURRENT REFERENCE ENERGY`, `CURRENT ENERGY` |
| `"mp2"` | `ccdata.mpenergies[-1][-1]` | adds `MP2 TOTAL ENERGY`, `MP2 CORRELATION ENERGY`; `CURRENT ENERGY` = MP2 total |
| `"ccsd"` | `ccdata.ccenergies[-1]` | adds `CCSD TOTAL ENERGY`, `CCSD CORRELATION ENERGY`; `CURRENT ENERGY` = CCSD total |
| `"ccsd(t)"` | `ccdata.ccenergies[-1]` | adds `CCSD(T) TOTAL ENERGY`, `CCSD(T) CORRELATION ENERGY`; `CURRENT ENERGY` = CCSD(T) total |
| `"b3lyp-d3"` (etc., dispersion suffix) | `scfenergies[-1]` + `dispersionenergies[-1]` | per `dispersion.md` |

**`mpenergies[-1][-1]`** retrieves the highest-order Møller–Plesset energy
cclib stored for the last step (e.g. `[[MP2]]` for an MP2-only job).

**`ccenergies[-1]`** retrieves the highest-level coupled-cluster energy:
the CCSD total for a CCSD job, the CCSD(T) total for a CCSD(T) job. cclib
does not separately store the intermediate CCSD value inside a CCSD(T) job,
so no CCSD qcvar is surfaced inside a CCSD(T) run.

**Error handling:** if the cclib attribute the method requires is absent
(`hasattr` is false) or empty (zero-length array), the harvester raises
`ValueError(...)` with a message naming the missing attribute. The runner
wraps this in `UnknownError`.

**Method-string normalisation:** `"ccsd(t,full)" → "ccsd(t)"` and
`"mp2(full)" → "mp2"` are applied at the top of `harvest()` so the above
branching matches. Dispersion suffixes (e.g. `"b3lyp-d3"`) are stripped to
recover the functional name; the unknown-functional branch falls through to
the SCF/DFT case.

## Hessian Parsing <!-- rq-c5b27d26 -->

The Hessian is parsed from the "Force constants in Cartesian coordinates:" block in
the Gaussian Freq log. This block contains the lower triangle of the symmetric
force constant matrix in Fortran D-exponent format, printed in column groups of 5.

The parser (`_parse_force_constants`) reconstructs the full `(3N, 3N)` symmetric matrix.

**Coordinate frame for force constants:** Gaussian prints force constants in the
**input orientation** (the coordinates from the `.com` file), NOT in the standard
orientation. This is significant for frame alignment:

- When `fix_com=True` and `fix_orientation=True`: the return molecule is `in_mol`
  (input frame) and the Hessian is already in the input frame — no rotation needed.
- When `fix_com=False` or `fix_orientation=False`: the return molecule is in the
  standard/calc frame. The Hessian must be transformed from input→calc frame using
  `in_mol.align(calc_mol)`.

## Regex Use <!-- rq-0743e952 -->

Log-text regex is used only for the following non-energy items:

| Item | Pattern | Why regex (not cclib) |
|------|---------|-----------------------|
| Normal termination | `"Normal termination of Gaussian"` substring | cclib doesn't surface a success/failure flag |
| Abnormal termination | absence of the above line, or `"Error termination"` | same |
| Gaussian version (for provenance) | `r"Gaussian\s+(\d+),\s+Revision\s+([\w.]+),"` | extracted from the startup banner before any cclib parse |
| Hessian (Cartesian, full 3N×3N) | `"Force constants in Cartesian coordinates:"` block, parsed by `_parse_force_constants()` | cclib does not populate `ccdata.hessian` for Gaussian (see "cclib Limitation: Full Cartesian Hessian not extracted") |

The Fortran D-exponent helper `_fortran_float()` is retained for the Hessian
block.

## Feature API <!-- rq-1a01df8f -->

### `harvest(in_mol: Molecule, method: str, log_content: str) -> Tuple[PreservingDict, Optional[np.ndarray], Optional[np.ndarray], Molecule]` <!-- rq-d9434480 -->

**Purpose**: parses the Gaussian log file and returns all QC data needed to construct an
`AtomicResult`.

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `in_mol` | `qcelemental.models.Molecule` | The input molecule (used for alignment and NRE cross-check). |
| `method` | `str` | Lowercase QCSchema method string (e.g. `"hf"`, `"mp2"`, `"b3lyp"`, `"ccsd(t,full)"`). Gaussian-specific sub-options like `,full` are normalised internally (e.g. `"ccsd(t,full)"` → `"ccsd(t)"` for energy extraction). |
| `log_content` | `str` | Full text of `gaussian.log`. |

**Returns** a 4-tuple `(qcvars, grad, hess, out_mol)`:

| Element | Type | Description |
|---------|------|-------------|
| `qcvars` | `PreservingDict` | QC variables (see below). |
| `grad` | `Optional[np.ndarray]` | Shape `(3*natoms,)` gradient in Hartree/Bohr, or `None`. |
| `hess` | `Optional[np.ndarray]` | Shape `(3*natoms, 3*natoms)` Hessian in Hartree/Bohr², or `None`. |
| `out_mol` | `Molecule` | Output molecule aligned to `in_mol`. |

**Raises**:
- `ValueError` if no coordinate information can be extracted from the log. <!-- rq-4336c3cd -->
- `ValueError` if the nuclear repulsion energy of the parsed geometry deviates from <!-- rq-90604835 -->
  `in_mol.nuclear_repulsion_energy()` by more than `1.0e-3` Hartree.
- Any `cclib` exception is allowed to propagate (caller wraps in `UnknownError`).

### QC Variables Populated <!-- rq-929a262a -->

> **Priority note:** The `properties` field on `AtomicResult` (populated via
> `build_atomicproperties`) is the proper QCSchema interface for consumers.
> `extras["qcvars"]` is retained for diagnostics and compatibility with other
> harnesses, but `properties` mappings are the primary contract.  All energy,
> gradient, and hessian variables listed below have existing mappings in
> `qcvars_to_atomicproperties`.  `"MULLIKEN CHARGES"` and `"DIPOLE MOMENT"`
> do not have `properties` mappings (the dipole mapping expects a vector
> `"HF DIPOLE"`, not the scalar magnitude stored here).

The following variables are placed in `qcvars` based on the method:

**All calculations**:
- `"NUCLEAR REPULSION ENERGY"` — from parsed geometry. <!-- rq-0753cd03 -->
- `"HF TOTAL ENERGY"` — from `ccdata.scfenergies[-1]` converted via `cclib.parser.utils.convertor("eV", "hartree")`. <!-- rq-37f40e92 -->
- `"SCF TOTAL ENERGY"` — same value as `"HF TOTAL ENERGY"` (required by standard suite contracts). <!-- rq-f59a07b8 -->
- `"CURRENT REFERENCE ENERGY"` — same value as `"HF TOTAL ENERGY"` (required by standard suite contracts). <!-- rq-0d2711f1 -->
- `"CURRENT ENERGY"` — set to the highest-level energy available (see below). <!-- rq-c5165f1f -->
- `"N ATOMS"` — from `ccdata.natom`. <!-- rq-463ef385 -->
- `"N BASIS FUNCTIONS"` — from `ccdata.nbasis` (when available). <!-- rq-17e8c0f7 -->
- `"N MOLECULAR ORBITALS"` — from `ccdata.nmo` (when available). <!-- rq-7f2c22fd -->
- `"N ALPHA ELECTRONS"` — from `ccdata.homos[0] + 1`. <!-- rq-6453046a -->
- `"N BETA ELECTRONS"` — from `ccdata.homos[-1] + 1` (or same as alpha for closed-shell). <!-- rq-0cdd386d -->

**Method-specific**:

| Method | Variables set |
|--------|--------------|
| `"hf"` or `"scf"` or DFT | `"HF TOTAL ENERGY"`, `"CURRENT ENERGY"` |
| `"mp2"` | `"MP2 CORRELATION ENERGY"` (= MP2 total − HF total), `"MP2 TOTAL ENERGY"`, `"CURRENT ENERGY"` |
| `"ccsd"` | `"CCSD CORRELATION ENERGY"`, `"CCSD TOTAL ENERGY"`, `"CURRENT ENERGY"` |
| `"ccsd(t)"` | `"CCSD(T) CORRELATION ENERGY"`, `"CCSD(T) TOTAL ENERGY"`, `"CURRENT ENERGY"` (cclib does not separately store the CCSD intermediate inside a CCSD(T) job, so no `"CCSD TOTAL ENERGY"` is set) |

**Properties driver extras** (when cclib attributes are available):
- `"MULLIKEN CHARGES"` — from `ccdata.atomcharges["mulliken"]`. <!-- rq-c7fd08b5 -->
- `"DIPOLE MOMENT"` — magnitude of `ccdata.moments[1]` in Debye. <!-- rq-943856cc -->

**Gradient** (when present):
- The negated forces array is transformed via `mill.align_gradient()` and returned as
  `grad` with shape `(3*natoms,)`. The runner stores this in `qcvars` as
  `"CURRENT GRADIENT"` and `"<METHOD> TOTAL GRADIENT"` with shape `(natoms, 3)`.

**Hessian** (when present):
- Parsed from the log's "Force constants in Cartesian coordinates" block (see
  "Hessian Parsing" above), transformed to the output molecule frame as needed,
  and returned as `hess` with shape `(3*natoms, 3*natoms)`. The runner sets
  `"CURRENT HESSIAN"` and `"<METHOD> TOTAL HESSIAN"` in `qcvars`.

### Output Molecule Construction and Frame Alignment <!-- rq-8d6fc024 -->

The output molecule is constructed from `ccdata.atomnos` and `ccdata.atomcoords[-1]`
(converted from Ångström to Bohr) as `calc_mol`. Note that cclib returns coordinates in
Gaussian's **standard orientation** (internally reoriented), not the input orientation.

Frame handling follows the same pattern as the GAMESS harvester:

- **`fix_com=True` and `fix_orientation=True`**: the return molecule is `in_mol`
  (preserving the input frame). An alignment mill is computed via
  `calc_mol.align(in_mol)` and applied to the gradient to transform it from standard
  orientation to the input frame. The Hessian (from direct log parsing) is already in
  the input frame and requires no rotation.

- **Otherwise** (free frame): the return molecule is produced by
  `in_mol.align(calc_mol)` (adopting the standard orientation frame). The gradient
  mill is identity (no transformation needed since cclib's gradient is already in
  standard orientation). The Hessian must be transformed from input→standard frame
  via `in_mol.align(calc_mol)`.

The `mill.align_gradient()` and `mill.align_hessian()` methods handle atom reordering
and rotation as needed.

## Gherkin Scenarios <!-- rq-058819da -->

```gherkin
Feature: Gaussian harvester

  Background:
    Given an AtomicInput for water (3 atoms: O, H, H) with charge=0, multiplicity=1

  # --- Energy parsing ---

  @rq-71b27e68
  Scenario: Parse HF total energy from a closed-shell log
    Given a Gaussian log for HF/STO-3G single-point on water
    And ccdata.scfenergies[-1] is -2041.3 eV
    When harvest(in_mol, "hf", log_content) is called
    Then qcvars["HF TOTAL ENERGY"] equals -2041.3 converted to Hartree
    And qcvars["CURRENT ENERGY"] equals qcvars["HF TOTAL ENERGY"]

  @rq-6bff3228
  Scenario: Parse DFT total energy (treated as HF SCF energy by cclib)
    Given a Gaussian log for B3LYP/6-31G* single-point on water
    And ccdata.scfenergies[-1] is -2082.1 eV
    When harvest(in_mol, "b3lyp", log_content) is called
    Then qcvars["HF TOTAL ENERGY"] equals -2082.1 converted to Hartree
    And qcvars["CURRENT ENERGY"] equals qcvars["HF TOTAL ENERGY"]

  @rq-578bb785
  Scenario: Parse MP2 total and correlation energies
    Given a Gaussian log for MP2/cc-pVDZ on water
    And ccdata.scfenergies[-1] is -2041.3 eV (HF energy)
    And ccdata.mpenergies[-1][-1] is -2044.8 eV (MP2 total energy)
    When harvest(in_mol, "mp2", log_content) is called
    Then qcvars["HF TOTAL ENERGY"] equals -2041.3 converted to Hartree
    And qcvars["MP2 TOTAL ENERGY"] equals -2044.8 converted to Hartree
    And qcvars["MP2 CORRELATION ENERGY"] equals the difference (in Hartree)
    And qcvars["CURRENT ENERGY"] equals qcvars["MP2 TOTAL ENERGY"]

  @rq-1b7ae62e
  Scenario: Parse CCSD total energy
    Given a Gaussian log for CCSD/cc-pVDZ on water
    And ccdata.scfenergies[-1] is -2041.3 eV
    And ccdata.ccenergies[-1] is -2045.2 eV (CCSD total energy)
    When harvest(in_mol, "ccsd", log_content) is called
    Then qcvars["CCSD TOTAL ENERGY"] equals -2045.2 converted to Hartree
    And qcvars["CCSD CORRELATION ENERGY"] equals the difference from HF (in Hartree)
    And qcvars["CURRENT ENERGY"] equals qcvars["CCSD TOTAL ENERGY"]

  @rq-eddef743
  Scenario: Parse CCSD(T) total energy
    Given a Gaussian log for CCSD(T)/cc-pVDZ on water
    And ccdata.ccenergies[-1] is the CCSD(T) total energy
    When harvest(in_mol, "ccsd(t)", log_content) is called
    Then qcvars["CCSD(T) TOTAL ENERGY"] is set
    And qcvars["CURRENT ENERGY"] equals qcvars["CCSD(T) TOTAL ENERGY"]

  # --- Gradient parsing ---

  @rq-36736115
  Scenario: Parse gradient and negate forces
    Given a Gaussian log for HF/STO-3G Force on water
    And ccdata.grads[-1] has shape (3, 3) with values representing forces (Hartree/Bohr)
    When harvest(in_mol, "hf", log_content) is called
    Then grad equals -ccdata.grads[-1] flattened to shape (9,)
    And each element of grad equals the negation of the corresponding force element

  @rq-3b0adf96
  Scenario: No gradient attribute in ccdata when running energy-only calculation
    Given a Gaussian log for HF/STO-3G single-point (no Force keyword)
    And ccdata has no grads attribute
    When harvest(in_mol, "hf", log_content) is called
    Then grad is None

  # --- Hessian parsing ---

  @rq-f2e32162
  Scenario: Parse Hessian from frequency calculation
    Given a Gaussian log for HF/STO-3G Freq on water
    And ccdata.hessian has shape (9, 9)
    When harvest(in_mol, "hf", log_content) is called
    Then hess has shape (9, 9)
    And hess equals ccdata.hessian

  @rq-344f3432
  Scenario: No Hessian when running energy-only or gradient calculation
    Given a Gaussian log for HF/STO-3G without Freq
    And ccdata has no hessian attribute
    When harvest(in_mol, "hf", log_content) is called
    Then hess is None

  # --- Properties parsing ---

  @rq-4c4a1542
  Scenario: Parse dipole moment for properties driver
    Given a Gaussian log with Pop=Full for HF/STO-3G on water
    And ccdata.moments[1] is a vector with magnitude 1.85 Debye
    When harvest(in_mol, "hf", log_content) is called
    Then qcvars["DIPOLE MOMENT"] is approximately 1.85

  @rq-ffae3fa6
  Scenario: Parse Mulliken charges for properties driver
    Given a Gaussian log with Pop=Full for HF/STO-3G on water
    And ccdata.atomcharges["mulliken"] is a list of 3 floats
    When harvest(in_mol, "hf", log_content) is called
    Then qcvars["MULLIKEN CHARGES"] is a list of length 3

  @rq-53f8ed12
  Scenario: Missing Mulliken charges do not raise an error
    Given a Gaussian log without Pop=Full
    And ccdata.atomcharges is empty or absent
    When harvest(in_mol, "hf", log_content) is called
    Then no error is raised
    And "MULLIKEN CHARGES" is not present in qcvars

  # --- Geometry and molecule output ---

  @rq-a76c3d7b
  Scenario: Output molecule geometry is constructed from ccdata
    Given ccdata.atomcoords[-1] contains the final Cartesian coordinates (Angstrom)
    And ccdata.atomnos contains [8, 1, 1] for water
    When harvest(in_mol, "hf", log_content) is called
    Then out_mol is a valid Molecule with symbols ["O", "H", "H"]

  @rq-a9194afb
  Scenario: Nuclear repulsion energy cross-check passes for consistent geometry
    Given the parsed geometry and in_mol have matching nuclear repulsion energies
    When harvest(in_mol, "hf", log_content) is called
    Then no ValueError is raised

  @rq-0b4cf843
  Scenario: Nuclear repulsion energy mismatch raises ValueError
    Given the parsed geometry has a nuclear repulsion energy differing from in_mol by > 1e-3 Hartree
    When harvest(in_mol, "hf", log_content) is called
    Then ValueError is raised

  # --- Error conditions ---

  @rq-b7c7b8cb
  Scenario: Log with no coordinate data raises ValueError
    Given a Gaussian log file from which cclib cannot extract atomic coordinates
    When harvest(in_mol, "hf", log_content) is called
    Then ValueError is raised

  @rq-bb775969
  Scenario: cclib parse failure propagates as-is
    Given a log_content string that cclib.io.ccopen cannot parse
    When harvest(in_mol, "hf", log_content) is called
    Then the cclib exception propagates to the caller

  # --- Missing cclib energy attributes ---

  @rq-28b5b1b3
  Scenario: Missing scfenergies for an HF/DFT job raises ValueError
    Given a cclib parse where ccdata has no "scfenergies" attribute (or it is empty)
    When harvest(in_mol, "hf", log_content) is called
    Then ValueError is raised
    And the message names "scfenergies"

  @rq-69740b1b
  Scenario: Missing mpenergies for an MP2 job raises ValueError
    Given a cclib parse where scfenergies is populated but mpenergies is absent
    When harvest(in_mol, "mp2", log_content) is called
    Then ValueError is raised
    And the message names "mpenergies"

  @rq-2abaea4c
  Scenario: Missing ccenergies for a CCSD job raises ValueError
    Given a cclib parse where scfenergies is populated but ccenergies is absent
    When harvest(in_mol, "ccsd", log_content) is called
    Then ValueError is raised
    And the message names "ccenergies"

  @rq-e05f5cfe
  Scenario: Missing ccenergies for a CCSD(T) job raises ValueError
    Given a cclib parse where scfenergies is populated but ccenergies is absent
    When harvest(in_mol, "ccsd(t)", log_content) is called
    Then ValueError is raised

  # --- Numerical consistency ---

  @rq-2816c5f6
  Scenario: SCF energy round-trips through cclib's convertor exactly
    Given a Gaussian log with "SCF Done: E(RHF) = -76.0571491640 A.U."
    And cclib stores ccdata.scfenergies[-1] = convertor(-76.0571491640, "hartree", "eV")
    When harvest(in_mol, "hf", log_content) is called
    Then qcvars["HF TOTAL ENERGY"] equals -76.0571491640 within 1 ULP (~1e-14 Hartree)

  @rq-dd7024df
  Scenario: MP2 total energy is read from mpenergies[-1][-1]
    Given a Gaussian MP2 log where ccdata.mpenergies[-1] has one entry MP2
    When harvest(in_mol, "mp2", log_content) is called
    Then qcvars["MP2 TOTAL ENERGY"] equals convertor(mpenergies[-1][-1], "eV", "hartree")

  @rq-4d653775
  Scenario: CCSD(T) total energy is read from ccenergies[-1]
    Given a Gaussian CCSD(T) log where ccdata.ccenergies[-1] is the CCSD(T) total
    When harvest(in_mol, "ccsd(t)", log_content) is called
    Then qcvars["CCSD(T) TOTAL ENERGY"] equals convertor(ccenergies[-1], "eV", "hartree")
    And qcvars does NOT contain "CCSD TOTAL ENERGY"

  @rq-44289c40
  Scenario: Dispersion energy is read from ccdata.dispersionenergies
    Given a B3LYP-D3 log where ccdata.dispersionenergies[-1] is the dispersion contribution (eV)
    When harvest(in_mol, "b3lyp-d3", log_content) is called
    Then qcvars["DISPERSION CORRECTION ENERGY"] equals convertor(dispersionenergies[-1], "eV", "hartree")

  # --- Source-level structure ---

  @rq-7edf875d
  Scenario: No log-text regex is used for SCF, MP*, CC, or dispersion energies
    Given an inspection of the harvester source
    Then the module contains no `_SCF_DONE_RE`, `_MP2_ENERGY_RE`, `_CCSD_ENERGY_RE`,
        `_CCSD_T_ENERGY_RE`, or `_DISPERSION_ENERGY_RE` constants
    And the Fortran D-exponent helper `_fortran_float()` is retained for Hessian parsing

  # --- Normal termination check (regex) ---

  @rq-90141f40
  Scenario: Normal termination line is present in a successful log
    Given a log file containing "Normal termination of Gaussian"
    When the log is inspected for normal termination
    Then the termination is classified as normal

  @rq-12c41d5e
  Scenario: Absent normal termination line is classified as abnormal
    Given a log file that does not contain "Normal termination of Gaussian"
    When the log is inspected for normal termination
    Then the termination is classified as abnormal
```
