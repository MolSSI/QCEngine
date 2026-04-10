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
| `ccdata.scfenergies[-1]` | eV | ÷ eV/Hartree | SCF/HF/DFT total energy |
| `ccdata.mpenergies[-1][-1]` | eV | ÷ eV/Hartree | MP2 total energy |
| `ccdata.ccenergies[-1]` | eV | ÷ eV/Hartree | CCSD or CCSD(T) total energy |
| `ccdata.grads[-1]` | Hartree/Bohr | negate (forces → gradient) | Gradient array |
| `ccdata.hessian` | Hartree/Bohr² | none | Hessian matrix |
| `ccdata.atomcoords[-1]` | Ångström | × Bohr/Å | Geometry for output molecule |
| `ccdata.atomnos` | atomic numbers | via qcel periodic table | Element symbols |
| `ccdata.atommasses` | amu | none | Atomic masses |
| `ccdata.moments[1]` | Debye | none | Dipole moment vector |
| `ccdata.atomcharges["mulliken"]` | e | none | Mulliken partial charges |

Energy conversion factor: `qcel.constants.conversion_factor("eV", "hartree")`.

**Important**: Gaussian reports atomic **forces** (= −dE/dR) in its output; cclib
preserves this sign convention in `ccdata.grads`. The gradient (= +dE/dR) is obtained by
**negating** the forces array:

```python
gradient = -ccdata.grads[-1]   # shape (natoms, 3) → flatten to (3*natoms,)
```

## Regex Fallback <!-- rq-c20746b1 -->

Regex is used for items cclib does not expose:

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
| `method` | `str` | Lowercase QCSchema method string (e.g. `"hf"`, `"mp2"`, `"b3lyp"`). |
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

The following variables are placed in `qcvars` based on the method:

**All calculations**:
- `"NUCLEAR REPULSION ENERGY"` — from parsed geometry. <!-- rq-0753cd03 -->
- `"HF TOTAL ENERGY"` — always set from `ccdata.scfenergies[-1]` (in Hartree). <!-- rq-37f40e92 -->
- `"CURRENT ENERGY"` — set to the highest-level energy available (see below). <!-- rq-c5165f1f -->

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
- The negated forces array (shape `(3*natoms,)`) is returned as `grad`; the runner sets
  `"CURRENT GRADIENT"` and `"<METHOD> TOTAL GRADIENT"` in `qcvars`.

**Hessian** (when present):
- The Hessian matrix (shape `(3*natoms, 3*natoms)`) is returned as `hess`; the runner
  sets `"CURRENT HESSIAN"` and `"<METHOD> TOTAL HESSIAN"` in `qcvars`.

### Output Molecule Construction <!-- rq-8d6fc024 -->

The output molecule is constructed from `ccdata.atomnos` and `ccdata.atomcoords[-1]`
(converted from Ångström to Bohr), then aligned to `in_mol` using
`in_mol.align(calc_mol, ...)` — the same pattern used by the GAMESS harvester.

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
