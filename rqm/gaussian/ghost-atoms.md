# Feature: Gaussian Harness — Ghost Atom Support <!-- rq-73c1f919 -->

This document specifies ghost-atom support for the Gaussian harness. A "ghost atom" is
an atom whose basis functions are placed at a position in space without a nucleus or
electrons — used primarily for counterpoise (BSSE) corrections and basis-set studies.
In a QCSchema `Molecule`, ghost atoms are marked by `real=False` at the matching index
in the `real` array; their entry in `symbols` gives the element whose basis is to be
placed at that position.

This feature replaces the prior policy in [[input-builder]] that raised `InputError`
on any molecule containing a ghost atom. The corresponding scenario in `input-builder.md`
(formerly `@rq-5053204b`) is removed when this feature is implemented; the runner
precheck at `qcengine/programs/gaussian/runner.py:183-184` is replaced by the build
logic specified below.

Related files:
- `input-builder.md` — `.com` file construction (modified by this feature).
- `harvester.md` — cclib-based output parsing (extended by this feature).
- `integration-tests.md` — integration test catalogue (extended by this feature).

## Background: Gaussian Ghost-Atom Syntax <!-- rq-dceb51f5 -->

Gaussian (16 and 09) supports ghost atoms in the molecule specification by appending
the suffix `-Bq` to the element symbol. The resulting atom contributes its basis
functions but has zero nuclear charge and contributes no electrons. The route line
and `model.basis` (e.g. `HF/6-31G*`) are unchanged; the element-based basis assignment
applies to the ghost atom exactly as if it were real.

Example atom block for an He···Ne system with a ghost Ne:

```
He    0.0000000000    0.0000000000    0.0000000000
Ne-Bq 2.5000000000    0.0000000000    0.0000000000
```

In Gaussian's standard/input-orientation table, ghost atoms appear with their
underlying element's atomic number in the "Atomic Number" column and a marker value
of `1000` in the "Atomic Type" column (real atoms have type `0`). cclib parses the
"Atomic Number" column but does **not** expose the "Atomic Type" column; consequently
`ccdata.atomnos` contains the underlying-element atomic numbers for both real and
ghost atoms (e.g. `[2, 10]` for He···Ne(ghost), not `[2, 0]`). Ghost atoms are
otherwise present in `ccdata.atomcoords` and `ccdata.grads` (with zero forces) in
input order.

Because cclib does not differentiate ghosts from real atoms, the harvester treats
`in_mol.real` as the sole source of truth for ghost designation.

## Feature API <!-- rq-a78c68a8 -->

This feature does not add new public functions or classes. It modifies the behavior of
existing entry points:

### Modified: `GaussianHarness.build_input(input_model, config, template=None)` <!-- rq-bd6cd9a4 -->

Located at `qcengine/programs/gaussian/runner.py`.

**Removed behavior**: the unconditional `InputError` for any molecule with
`real=False` (rq-5053204b) is removed.

**New behavior**:

- The atom block (see [[input-builder]] §"Molecule Block Construction") is constructed
  by iterating jointly over `mol.symbols`, `mol.real`, and the per-atom Ångström
  coordinates.
- For each atom `i`:
  - If `mol.real[i]` is `True`, the line is `"<symbol> <x> <y> <z>"` (unchanged).
  - If `mol.real[i]` is `False`, the line is `"<symbol>-Bq <x> <y> <z>"`.
- Coordinate formatting (≥10 decimal places) is unchanged.
- All other route-line and Link-0 behavior is unchanged.

**Validation**: no additional pre-validation is performed. Molecules that consist
entirely of ghost atoms, or that combine ghosts with method/driver combinations
unsupported by Gaussian, are passed through to Gaussian; failures are surfaced via the
existing log-pattern error dispatch (see [[runner]] §"Error dispatch").

### Modified: `harvest(in_mol, method, log_content)` <!-- rq-f1d712bb -->

Located at `qcengine/programs/gaussian/harvester.py`.

The harvester is extended so that the parsed output molecule preserves the ghost
designation present in `in_mol`:

- **Output molecule construction**: `out_mol.real` is taken from `in_mol.real` rather
  than defaulted to all-real. When ghosts are present, `out_mol.symbols` is taken
  from `in_mol.symbols` rather than from cclib's `atomnos` (since `in_mol` is
  authoritative for both ghost designation and element identity).

- **Frame alignment** (cf. [[harvester]] §"Output Molecule Construction and Frame
  Alignment"): the `align()` calls between `in_mol` and the cclib-parsed geometry use
  `generic_ghosts=True`, matching the pattern already used by the NWChem harvester
  (`qcengine/programs/nwchem/harvester.py:1260`).

- **Nuclear-repulsion-energy cross-check** (cf. [[harvester]] §"Feature API → Raises",
  rq-90604835): when ghost atoms are present, the cross-check compares NREs computed
  over real atoms only. Concretely, construct real-only molecules from `in_mol` and
  from the parsed geometry, compute the NRE of each via
  `Molecule.nuclear_repulsion_energy()`, and apply the existing `1.0e-3` Hartree
  tolerance. The behavior for molecules with no ghosts is unchanged.

- **Atom-count sanity check**: after parsing, when ghost atoms are present in
  `in_mol`, the harvester verifies that `ccdata.natom == len(in_mol.symbols)`
  (using cclib's `natom` integer attribute rather than `len(ccdata.atomnos)`).
  A mismatch indicates that Gaussian dropped or duplicated atoms during input parsing,
  or that cclib's parse drifted; the harvester raises `ValueError` with a descriptive
  message. (This check is performed only when ghosts are present, since it is the
  cheapest cross-check available given that cclib cannot itself flag ghost atoms.)

- **Gradient and Hessian transformations** are unchanged in form but operate over the
  full atom array (real + ghost). Gaussian prints zero forces for ghost atoms; cclib
  preserves this. Existing `mill.align_gradient()` / `mill.align_hessian()` calls
  handle the all-atom arrays correctly when alignment is computed with
  `generic_ghosts=True`.

## Method and Driver Compatibility <!-- rq-4091a959 -->

Gaussian itself imposes restrictions on which method+driver combinations may be used
with ghost atoms (e.g. some post-HF analytic frequency methods do not support ghost
atoms). This harness does **not** pre-validate such combinations. When Gaussian
rejects a ghost-bearing input at runtime, the existing error dispatch in
`GaussianHarness.compute()` will surface the failure as `InputError` or `UnknownError`
based on the log pattern, per [[runner]].

## Gherkin Scenarios <!-- rq-84b093c1 -->

```gherkin
Feature: Gaussian harness ghost-atom support

  Background:
    Given a TaskConfig with ncores=4, memory=8.0, scratch_directory=None, scratch_messy=False

  # --- Input builder: atom-block formatting ---

  @rq-2229f087
  Scenario: Real atom emits a plain element symbol in the atom block
    Given an AtomicInput whose molecule has one atom with symbol "H" and real=True
    When build_input() is called
    Then the gaussian.com atom block contains a line beginning with "H " (not "H-Bq ")

  @rq-916e4919
  Scenario: Ghost atom emits a "-Bq" suffix in the atom block
    Given an AtomicInput whose molecule has one atom with symbol "H" and real=False
    When build_input() is called
    Then the gaussian.com atom block contains a line beginning with "H-Bq "

  @rq-2a015edf
  Scenario: Mixed real and ghost atoms in one molecule
    Given an AtomicInput for He at origin (real=True) and Ne at (2.5, 0, 0) (real=False)
    When build_input() is called
    Then the atom block contains a line beginning with "He "
    And the atom block contains a line beginning with "Ne-Bq "
    And the atoms appear in input order

  @rq-a100502c
  Scenario: Route line is unchanged when ghost atoms are present
    Given an AtomicInput for HF/6-31G* with two atoms, one of which has real=False
    When build_input() is called
    Then the route line is exactly "#P HF/6-31G*"

  @rq-ce76f384
  Scenario: All-ghost molecule produces an atom block with every line "-Bq"-suffixed
    Given an AtomicInput whose molecule has two atoms, both with real=False
    When build_input() is called
    Then every atom-block coordinate line contains "-Bq"
    And build_input() does not raise

  @rq-ceef73bd
  Scenario: Coordinate precision is preserved for ghost atoms
    Given an AtomicInput whose molecule has a ghost atom at (1.23456789012, 0, 0) Å
    When build_input() is called
    Then the ghost atom's coordinate line contains at least 10 decimal digits

  # --- Harvester: output molecule construction ---

  @rq-b27d4aff
  Scenario: out_mol.real is copied from in_mol when ghosts are present
    Given an AtomicInput for He (real=True) and Ne (real=False)
    And a Gaussian log produced by the corresponding HF/STO-3G calculation
    When harvest(in_mol, "hf", log_content) is called
    Then out_mol.real equals [True, False]

  @rq-731f01ba
  Scenario: out_mol.symbols are copied from in_mol for ghost positions
    Given an AtomicInput for He (real=True) and Ne (real=False)
    And a Gaussian log in which cclib reports atomnos == [2, 10]
    When harvest(in_mol, "hf", log_content) is called
    Then out_mol.symbols equals ["He", "Ne"]

  @rq-7411efa0
  Scenario: Atom count sanity check passes when cclib atom count matches in_mol
    Given an in_mol with exactly one ghost atom (length 2)
    And a Gaussian log where ccdata.natom == 2
    When harvest(in_mol, "hf", log_content) is called
    Then no ValueError is raised

  @rq-66067019
  Scenario: Atom count mismatch raises ValueError when ghosts are present
    Given an in_mol with exactly one ghost atom (length 2)
    And a Gaussian log where ccdata.natom != 2
    When harvest(in_mol, "hf", log_content) is called
    Then ValueError is raised

  # --- Harvester: nuclear repulsion energy cross-check ---

  @rq-a69d0898
  Scenario: NRE cross-check excludes ghost atoms when present
    Given an in_mol containing real atom He and ghost atom Ne separated by 2.5 Å
    And a Gaussian log whose standard-orientation geometry matches in_mol
    When harvest(in_mol, "hf", log_content) is called
    Then the NRE cross-check uses only the real atoms of both molecules
    And no ValueError is raised

  @rq-b57c882b
  Scenario: NRE cross-check still catches geometry mismatches with ghosts present
    Given an in_mol containing one real atom and one ghost atom
    And a Gaussian log whose parsed real-atom geometry differs from in_mol
    And the real-atom-only NRE of the parsed geometry differs by > 1e-3 Hartree from in_mol's
    When harvest(in_mol, "hf", log_content) is called
    Then ValueError is raised

  @rq-e5a98119
  Scenario: NRE cross-check is unchanged for all-real molecules
    Given an in_mol with all real atoms and a matching parsed geometry
    When harvest(in_mol, "hf", log_content) is called
    Then the NRE cross-check is performed on the full molecules (no real-atom filtering needed)
    And no ValueError is raised

  # --- Harvester: gradient with ghost atoms ---

  @rq-d97cc860
  Scenario: Gradient array spans real and ghost atoms in input order
    Given an AtomicInput for He···Ne (He real, Ne ghost) with driver=gradient
    And a Gaussian Force log in which forces are printed for all 2 atoms
    When harvest(in_mol, "hf", log_content) is called
    Then grad has shape (6,)
    And the gradient on the ghost atom is approximately zero

  # --- Harvester: alignment with generic_ghosts ---

  @rq-f8ce079b
  Scenario: align() calls use generic_ghosts=True when ghosts are present
    Given an in_mol with at least one ghost atom
    When harvest(in_mol, "hf", log_content) is called
    Then every Molecule.align() invocation in the harvester passes generic_ghosts=True

  # --- Round-trip via compute() ---

  @rq-04028610
  Scenario: AtomicResult.molecule preserves ghost designation
    Given an AtomicInput for He (real) and Ne (ghost) with HF/aug-cc-pVDZ energy
    When GaussianHarness.compute() runs and terminates normally
    Then the returned AtomicResult.molecule.real equals [True, False]
    And the returned AtomicResult.molecule.symbols equals ["He", "Ne"]

  @rq-24d73182
  Scenario: Energy of a ghost-bearing system matches a reference value
    Given the He···Ne system from test_simple_ghost (He real at origin, Ne ghost at (2.5, 0, 0) Å)
    And HF/aug-cc-pVDZ with nocom/noreorient enforced
    When GaussianHarness.compute(driver="energy") is called
    Then return_result is approximately -2.8557143339397539 Hartree (atol 1e-6)
    And properties.nuclear_repulsion_energy is approximately 0.0 (atol 1e-6)

  @rq-04ed684d
  Scenario: Gradient of a ghost-bearing system matches a reference vector
    Given the same He···Ne system
    When GaussianHarness.compute(driver="gradient") is called
    Then return_result has shape (2, 3)
    And the gradient on the real He atom matches the reference vector (atol 1e-6)
    And the gradient on the ghost Ne atom is approximately zero (atol 1e-6)

  # --- test_tricky_ghost: eneyne + Ne, MP2 ---

  @rq-c40101e6
  Scenario: Method "mp2(full)" maps to "MP2(Full)" in the Gaussian route line
    Given AtomicInput.specification.model.method is "mp2(full)"
    When muster_modelchem("mp2(full)", "energy", 1) is called
    Then the returned method_string is "MP2(Full)"

  @rq-a9026c9f
  Scenario: Method "mp2(full)" is normalised back to "mp2" by the harvester
    Given a Gaussian log produced by a method "mp2(full)" calculation
    When harvest(in_mol, "mp2(full)", log_content) is called
    Then the MP2-specific qcvars (MP2 TOTAL ENERGY, MP2 CORRELATION ENERGY) are populated
    And qcvars["CURRENT ENERGY"] equals MP2 TOTAL ENERGY

  @rq-8fa7c4c3
  Scenario: Gaussian "NoSymm" keyword forces C1 symmetry in the route line
    Given user_keywords = {"NoSymm": ""}
    When build_route_line is invoked
    Then the route line contains the bare token "NoSymm"

  @rq-3257a6f5
  Scenario: MP2/6-31G* energy on eneyne dimer matches consensus reference
    Given the dimer subject from eneyne_ne_qcschemamols() (10 real atoms, 0 ghosts)
    And keywords = {"NoSymm": ""} and method = "mp2(full)"
    When qcng.compute is called with program "gaussian", method "mp2", basis "6-31g*", driver "energy"
    Then the result is successful
    And properties.nuclear_repulsion_energy matches bimol_ref["eneyne"]["nre"]["dimer"] (atol 1e-4)
    And properties.calcinfo_nbasis equals 72
    And properties.calcinfo_nmo equals 72
    And properties.return_energy matches bimol_ref["eneyne"]["mp2"]["dimer"] (atol 3e-6)

  @rq-609268f8
  Scenario: MP2/6-31G* gradient on eneyne dimer matches consensus reference
    Given the dimer subject from eneyne_ne_qcschemamols()
    And keywords = {"NoSymm": ""} and method = "mp2(full)"
    When qcng.compute is called with driver "gradient"
    Then return_result equals bimol_ref["eneyne"]["mp2_gradient"]["dimer"] (atol 3e-6)

  @rq-d4e60dd2
  Scenario: MP2/6-31G* energy on mAgB (real ethylene + ghost ethyne) matches reference
    Given the mAgB subject (6 real + 4 ghost atoms)
    And keywords = {"NoSymm": ""} and method = "mp2(full)"
    When qcng.compute is called with driver "energy"
    Then properties.return_energy matches bimol_ref["eneyne"]["mp2"]["mAgB"] (atol 3e-6)

  @rq-31334859
  Scenario: MP2/6-31G* energy on gAmB (ghost ethylene + real ethyne) matches reference
    Given the gAmB subject (4 real + 6 ghost atoms)
    And keywords = {"NoSymm": ""} and method = "mp2(full)"
    When qcng.compute is called with driver "energy"
    Then properties.return_energy matches bimol_ref["eneyne"]["mp2"]["gAmB"] (atol 3e-6)

  @rq-a5593d5c
  Scenario: MP2/6-31G* on mB (linear ethyne alone) returns the correct energy
    Given the mB subject (4 real atoms, linear)
    And keywords = {"NoSymm": ""} and method = "mp2(full)"
    When qcng.compute is called with driver "energy"
    Then the result is successful
    And properties.return_energy matches bimol_ref["eneyne"]["mp2"]["mB"] (atol 3e-6)

  @rq-7e1b0bb4
  Scenario: test_tricky_ghost accepts Gaussian's linear PG for the mB subject
    Given qcprog is "gaussian" and subject is "mB"
    And the test has parsed pg "C*v" from atres.stdout
    When the test reaches the point-group assertion
    Then the Gaussian/mB branch is taken (assert pg == "C*v"), not the default ref branch

  @rq-6121b3e1
  Scenario: test_tricky_ghost uses the default ref branch for non-mB Gaussian subjects
    Given qcprog is "gaussian" and subject is "dimer"
    And the test has parsed pg "C2v" from atres.stdout
    When the test reaches the point-group assertion
    Then the default branch is taken (assert pg in ref["pg"][subject])
```

## Integration Tests <!-- rq-c4439b82 -->

### `test_simple_ghost` (He···Ne) <!-- rq-95dab09b -->

`qcengine/programs/tests/test_ghost.py::test_simple_ghost` already parametrises a He
real / Ne ghost system across cfour, gamess, nwchem, and psi4. As part of this feature,
Gaussian is added to the parametrisation:

```python
pytest.param("gaussian", "aug-cc-pvdz", {}, marks=using("gaussian")),
```

The integration test asserts:
- `nuclear_repulsion_energy` ≈ 0
- `nbasis == 32` and `nmo == 32`
- `return_energy` ≈ `-2.8557143339397539` (atol `1e-6`)
- `return_result` for the gradient driver matches the reference gradient (atol `1e-6`)

No new integration-test infrastructure is required; the existing parametrisation is
extended in place.

### `test_tricky_ghost` (eneyne+Ne, MP2) <!-- rq-2cb2c5c0 -->

`qcengine/programs/tests/test_ghost.py::test_tricky_ghost` exercises MP2/6-31G* on the
ethylene–ethyne–Ne system across five "subjects" (`dimer`, `mA`, `mB`, `mAgB`, `gAmB`)
that vary which fragments are real vs. ghost, and across `energy` and `gradient`
drivers. As part of this feature, Gaussian is added to the `qcprog` parametrisation:

```python
pytest.param(
    "gaussian",
    "6-31g*",
    {"NoSymm": ""},
    id="gaussian",
    marks=using("gaussian"),
),
```

**Per-program method override**: the test overrides the QCSchema method for
Gaussian only (other programs continue to use `"mp2"`):

```python
method = "mp2(full)" if qcprog == "gaussian" else "mp2"
```

**Keyword rationale**:

- `"NoSymm": ""` — Gaussian uses symmetry to standardise the molecular orientation
  by default; `NoSymm` (emitted as a bare token by `build_route_line`) tells
  Gaussian not to use symmetry in the calculation. With `fix_com=True` and
  `fix_orientation=True` already imposed on the input molecule, `NoSymm` keeps
  Gaussian's internal frame aligned with the input frame for clean integration
  accuracy. Note that `NoSymm` does **not** suppress the geometric point-group
  printout — Gaussian still prints e.g. `Full point group  C*V` for linear
  systems regardless.

**Method-spec rationale (`mp2(full)`)**:

Gaussian's MP2 freezes core orbitals by default (`NFC=4` for the
ethylene+ethyne system, ~9.5 mHa offset from the all-electron reference). Two
syntaxes exist for forcing all-electron MP2:

1. `MP2=Full` as a separate route keyword. Works for energy, but Gaussian
   refuses to compute analytic gradients with this form:
   `"Post-SCF densities or gradients only with Real MP2, MP3, MP4SDQ, CI, CCD, and QCI."`
2. `MP2(Full)/basis` with the option attached to the method spec inline. Works
   for both energy and gradient.

This feature adds a new `_METHOD_MAP` entry to `germinate.py` that maps the
lower-case method string `"mp2(full)"` to the Gaussian token `"MP2(Full)"`,
mirroring the existing `"ccsd(t,full)" → "CCSD(T,Full)"` convention. The
harvester's `_METHOD_NORMALISE` dict in `harvester.py` is also extended so
`"mp2(full)"` is normalised to `"mp2"` for energy-extraction purposes (the
EUMP2 regex matches both forms).

**Reference values** are reused from the existing `bimol_ref["eneyne"]` dict with
the existing tolerances (`atol = 1.0e-4` for NRE, `atol = 3.0e-6` for energy and
return-result comparison). If empirical verification during implementation reveals
a systematic Gaussian-specific offset > 3.0e-6, the test should be amended to use a
Gaussian-specific reference (with a comment explaining why); otherwise the existing
references stand.

**Point-group assertion**: the `pgline` dict at the bottom of `test_tricky_ghost`
gains a Gaussian regex:

```python
"gaussian": r"Full point group\s+(?P<pg>\S+)",
```

The `\S+` form (rather than `\w+`) is required because Gaussian writes
non-word characters (e.g. `*`) in linear point-group labels such as `C*V`.

Gaussian's natural PG labels match the existing `bimol_ref["eneyne"]["pg"]`
list for four of the five subjects (after the `pg.capitalize()` normalisation
that already runs in the test):

| subject | Gaussian | ref |
|---------|----------|------|
| dimer | C2V → C2v | "C2v" ✓ |
| mA | C2V → C2v | ["D2h", "C2v"] ✓ |
| mB | C*V → C*v | ["D4h", "C4v", "C2v"] ✗ |
| mAgB | C2V → C2v | "C2v" ✓ |
| gAmB | C2V → C2v | "C2v" ✓ |

Only the `mB` subject (ethyne alone, linear) is problematic: Gaussian reports
the natural linear point group `C*v`, whereas the existing consensus list
records the abelian subgroups other programs reduce it to. The test grows a
single Gaussian-specific branch that accepts `C*v` for `mB`:

```python
if qcprog == "gamess":
    assert pg == "C1"
elif qcprog == "gaussian" and subject == "mB":
    assert pg == "C*v"
else:
    assert pg in ref["pg"][subject]
```

This preserves the original PG check for all other Gaussian subjects while
avoiding a Gaussian-specific edit to the cross-program `bimol_ref["eneyne"]["pg"]`
dict.

**Out-of-scope**: no Gaussian-specific reference-data file is created. No harness
changes are required — the existing ghost-atom support (this document) plus the
existing user-keyword pass-through (`input-builder.md`) suffice.

## Migration Notes <!-- rq-df7eb0c3 -->

- The runner check at `qcengine/programs/gaussian/runner.py:183-184` (the
  `InputError("The Gaussian harness does not support ghost atoms (real=False).")`
  block) is removed in its entirety. The corresponding rq-5053204b scenario in
  `input-builder.md` is removed at the same time.
- `input-builder.md` §"Molecule Block Construction" is updated to specify the
  `<symbol>-Bq` formatting for ghost atoms.
- `harvester.md` is updated so that the output-molecule construction, NRE cross-check,
  and alignment sections accommodate ghost atoms as specified above.
- No public API surface changes; consumers who previously avoided ghost atoms will see
  no behavioral change.
