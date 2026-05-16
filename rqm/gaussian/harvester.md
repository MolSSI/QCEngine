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
| `ccdata.scfenergies[-1]` | eV | `cclib.parser.utils.convertor("eV"→"hartree")` (fallback only) | SCF/HF/DFT total energy |
| `ccdata.mpenergies[-1][-1]` | eV | `cclib...convertor` (fallback only) | MP2 total energy |
| `ccdata.ccenergies[-1]` | eV | `cclib...convertor` (fallback only) | CCSD or CCSD(T) total energy |
| `ccdata.grads[-1]` | Hartree/Bohr | negate (forces → gradient) | Gradient array |
| `ccdata.atomcoords[-1]` | Ångström | × Bohr/Å | Geometry for output molecule (standard orientation) |
| `ccdata.atomnos` | atomic numbers | via qcel periodic table | Element symbols (no-ghost branch only) |
| `ccdata.natom` | integer | none | Atom count for the sanity check and downstream sizing (preferred over `len(atomnos)`) |
| `ccdata.nbasis` | integer | none | Number of basis functions |
| `ccdata.nmo` | integer | none | Number of molecular orbitals |
| `ccdata.homos` | integer array | +1 for count | N alpha/beta electrons |
| `ccdata.moments[1]` | Debye | none | Dipole moment vector |
| `ccdata.atomcharges["mulliken"]` | e | none | Mulliken partial charges |

### cclib Limitations <!-- rq-f6537c40 -->

Two significant cclib limitations require direct log parsing as the primary strategy
for energies and the Hessian:

**1. Energy precision loss when reversing cclib's eV storage with a foreign factor:**

cclib internally stores energies in electron volts. It parses the Hartree value
from the Gaussian log, multiplies by its internal Hartree→eV factor
(`27.21138505`), and stores the result on `ccdata.scfenergies`, `ccdata.mpenergies`,
`ccdata.ccenergies`, and `ccdata.dispersionenergies`. If the harvester reverses
this with qcelemental's factor (`qcel.constants.conversion_factor("eV", "hartree")`
= `0.036749322481536`, corresponding to `27.21138602` Hartree→eV — the values
disagree at the 8th decimal place because they were derived from different
CODATA cycles), the round-trip introduces a relative error of approximately
`3.6e-9`. For typical SCF energies (~100 Hartree) that's an absolute error of
~3.6e-6 Hartree, exceeding the cross-program tolerances (`atol = 2e-7`) and
the empirical-dispersion test tolerances.

**Resolution (two layers):**

1. *Primary path:* energies are parsed directly from the Gaussian log via regex
   (see "Direct Energy Parsing" below). This preserves Gaussian's full printed
   precision (up to 12 significant digits) and avoids cclib's eV detour entirely.
2. *Fallback path:* when the regex misses, the harvester reads the eV-stored
   value from `ccdata` and converts back via `cclib.parser.utils.convertor(x,
   "eV", "hartree")` — cclib's own converter, which uses the same `27.21138505`
   factor in both directions. The forward/reverse pair is structurally exact
   (single factor, both directions in the same `_convertor` dict), so the
   round-trip introduces zero round-trip error and is forward-compatible with
   any future CODATA refresh inside cclib.

**2. Full Cartesian Hessian not extracted by cclib:**

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

## Direct Energy Parsing <!-- rq-c20746b1 -->

Energies are parsed directly from the Gaussian log in Hartree to avoid the cclib eV
round-trip. The following regex patterns are used:

| Energy | Regex | Format |
|--------|-------|--------|
| SCF/HF/DFT | `r"SCF Done:\s+E\(\S+\)\s*=\s*([-\d.]+)"` | Plain float |
| MP2 | `r"EUMP2\s*=\s*([-\d.DE+]+)"` | Fortran D-exponent |
| CCSD(T) | `r"CCSD\(T\)=\s*([-\d.DE+]+)"` | Fortran D-exponent |
| CCSD | `r"DE\(Corr\)=\s*[-\d.]+\s+E\(CORR\)=\s*([-\d.]+)"` (last match) | Plain float |

Fortran D-exponent values (e.g. `-0.76276030623D+02`) are converted by replacing
`D` with `E` before calling `float()`.

If a regex does not match, the harvester falls back to cclib's eV-based energy with
the qcelemental conversion factor.

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

## Other Regex Items <!-- rq-0743e952 -->

| Item | Pattern |
|------|---------|
| Normal termination | `"Normal termination of Gaussian"` in log |
| Abnormal termination | absence of the above line, or `"Error termination"` |
| Gaussian version (for provenance) | `r"Gaussian\s+(\d+),\s+Revision\s+([\w.]+),"` |

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
- `"HF TOTAL ENERGY"` — from direct log parsing (SCF Done line), falling back to `ccdata.scfenergies[-1]` converted from eV. <!-- rq-37f40e92 -->
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
| `"ccsd(t)"` | `"CCSD CORRELATION ENERGY"`, `"CCSD TOTAL ENERGY"`, `"CCSD(T) CORRELATION ENERGY"`, `"CCSD(T) TOTAL ENERGY"`, `"CURRENT ENERGY"` |

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
