# Feature: Gaussian Harness — Atom-Label Tolerance <!-- rq-0dab9c40 -->

This document specifies the Gaussian harness's behavior when the input molecule carries
non-empty `atom_labels` — QCSchema metadata that distinguishes same-element atoms
(e.g. `"H"` vs `"H_other"` vs `"H_4sq"`). It also adds Gaussian to the
`test_atom_labels` integration test in `test_ghost.py`.

In QCSchema, `Molecule.atom_labels` is a per-atom list of strings parallel to
`Molecule.symbols`. Labels are derived from inputs like `H5` → element `"H"`, label
`"5"`; `H_other` → element `"H"`, label `"_other"`. The labels carry no chemistry on
their own — atoms with the same element get the same basis from the route-line
`<method>/<basis>` specification, regardless of label — but downstream consumers and
some QC programs use them to break symmetry-driven solution degeneracies or to apply
per-atom keyword overrides.

The chosen scope for this feature: **the Gaussian harness treats `atom_labels` as
informational metadata that does not affect input generation**. The `.com` atom block
continues to emit only the bare element symbol (with the `-Bq` suffix for ghosts).
Atom labels are not propagated into the `.com` file. The harness must, however,
*accept* molecules that carry non-empty labels without error, and the resulting
calculation must produce results comparable to what other harnesses produce.

Related files:
- `input-builder.md` — `.com` file construction (no changes; this feature confirms
  the existing behavior is correct for labeled inputs).
- `ghost-atoms.md` — ghost atom support (separate feature; this one is independent
  but uses the same atom-block formatting layer).
- `integration-tests.md` — cross-program integration test catalogue (declared out
  of scope for `test_ghost.py`; this feature extends `test_ghost.py::test_atom_labels`
  specifically).

## Design Rationale <!-- rq-5d5df0b4 -->

Gaussian itself accepts arbitrary suffixes on atom-spec lines (e.g. `H_other` is a
valid atom-spec, parsed as element `H` with type identifier `_other`), but our test
target — `test_atom_labels` — does not exercise that distinction. All 4 hydrogen
atoms get the same `aug-cc-pVDZ` basis regardless of label. Emitting bare `H` lines
is therefore sufficient and avoids two classes of risk:

1. **Label-syntax edge cases**: QCElemental allows labels containing characters
   (`_`, digits) that may or may not parse cleanly across Gaussian versions and
   release patch levels.
2. **Symmetry-equivalence surprises**: with bare element symbols, Gaussian sees four
   chemically-identical hydrogens and may apply its own symmetry detection. With
   distinct labels (e.g. `H1`, `H2`), Gaussian may instead apply per-type treatment
   that subtly changes the SCF guess. Keeping labels out of the `.com` keeps
   results comparable to the existing cross-program references.

This matches the de facto convention in the rest of QCEngine: neither the CFOUR,
NWChem, nor Psi4 harness propagates `atom_labels` into its respective input format
(grep for `atom_label` in those harnesses returns nothing).

## Feature API <!-- rq-824f4a96 -->

No new public functions or classes are introduced. The behavior of existing entry
points is documented (not changed):

### `GaussianHarness.build_input` — confirmed behavior for labeled molecules <!-- rq-03c1f71c -->

- The atom block iterates `mol.symbols` and `mol.real`; `mol.atom_labels` is
  **not** read.
- For a molecule with `symbols = ["H", "H", "H", "H"]`,
  `atom_labels = ["", "5", "_other", "_4sq"]`, and all `real=True`, the atom
  block contains four lines each beginning with `"H "` and no occurrences of the
  label strings.
- No `InputError` is raised on non-empty `atom_labels`.

### `harvest` — confirmed behavior for labeled molecules <!-- rq-a203a314 -->

- The output molecule (`out_mol`) returned by `harvest()` carries whatever
  `atom_labels` the input molecule had when `fix_com=True` and
  `fix_orientation=True` (because `out_mol` is then `in_mol`).
- When `fix_com=False` or `fix_orientation=False`, `out_mol` is reconstructed from
  cclib output and carries default (empty) `atom_labels`. This matches the existing
  behavior for non-labeled molecules and does not change with this feature.

## Integration Test <!-- rq-4f50add1 -->

### `test_atom_labels` (`test_ghost.py`) <!-- rq-4061edf7 -->

`qcengine/programs/tests/test_ghost.py::test_atom_labels` runs MP2/aug-cc-pVDZ on
a 4-hydrogen system (`H`, `H5`, `H_other`, `H_4sq`) at Bohr coordinates
`(0,0,0)`, `(5,0,0)`, `(0,5,0)`, `(5,5,0)`. Gaussian is added to the `qcprog`
parametrisation:

```python
pytest.param("gaussian", "aug-cc-pvdz", {}, marks=using("gaussian")),
```

**Keyword rationale**: empty `keywords = {}` should suffice. No frozen-core
concern (H has no core electrons). Default Gaussian SCF/MP2 convergence is at the
1e-8/1e-6 level, well within the test's `atol = 3.0e-6` tolerance.

**Reference values**: the test already encodes two acceptable solutions:
- "Primary" SCF solution: `scf = -1.656138508`, `mp2 = -1.7926264513`
- "Alternative" solution (currently labeled `nwchem`): `scf_alt = -1.705577613`,
  `mp2_alt = -1.870251459939`

Gaussian must converge to one of these two solutions. The existing `try`/`except`
block at the bottom of `test_atom_labels` already accepts both, but its `else`
branch (the alternative-solution path) gates on `qcprog == "nwchem"`. As part of
this feature, that gate is relaxed to also accept `qcprog == "gaussian"`:

```python
except AssertionError as exc:
    if qcprog in ("nwchem", "gaussian"):
        # accept alternative solution
        ...
    else:
        raise AssertionError from exc
```

If empirical Gaussian results land cleanly on the primary solution, this
relaxation may turn out to be unused; it is added defensively because 4-H
rectangle systems are notoriously sensitive to the SCF guess (the symmetric and
broken-symmetry solutions are nearly degenerate).

**`nre` and `nmo` references**: the existing `nre = 1.0828427` and `nmo = 36`
values are program-independent (NRE is a property of the geometry, `nmo = 4 atoms
× 9 orbitals/atom = 36` is a property of the basis) and apply unchanged to
Gaussian.

**No harness changes**: this test passes by virtue of the parametrisation
addition alone.

## Gherkin Scenarios <!-- rq-4bd911ec -->

```gherkin
Feature: Gaussian harness atom-label tolerance

  Background:
    Given a TaskConfig with ncores=1, memory=1.0, scratch_directory=None

  # --- Input builder: atom-block formatting with labels ---

  @rq-20063525
  Scenario: Labeled hydrogen atom emits a bare "H " line (no label suffix)
    Given an AtomicInput whose molecule has one atom with symbol "H" and atom_label "_other"
    When build_input() is called
    Then the gaussian.com atom block contains a line beginning with "H "
    And no atom-block line contains the substring "_other"

  @rq-da7df6e9
  Scenario: Multiple same-element atoms with distinct labels each emit bare element lines
    Given an AtomicInput whose molecule has 4 H atoms with labels ["", "5", "_other", "_4sq"]
    When build_input() is called
    Then the gaussian.com atom block contains exactly 4 lines beginning with "H "
    And no atom-block line contains any of the substrings "H5", "H_other", "H_4sq"

  @rq-7a592160
  Scenario: Empty atom_labels (all-default) preserves existing behavior
    Given an AtomicInput whose molecule has atom_labels = ["", "", "", ""]
    When build_input() is called
    Then the gaussian.com atom block content is identical to what it would be if atom_labels were omitted entirely

  @rq-090d144e
  Scenario: Labels do not appear in the route line
    Given an AtomicInput with non-empty atom_labels
    When build_input() is called
    Then the route line of the resulting .com file does not contain any atom-label substring

  @rq-a1bc15af
  Scenario: build_input() accepts a labeled molecule without raising
    Given an AtomicInput with non-empty atom_labels
    When build_input() is called
    Then no exception is raised

  # --- Atom-label and ghost-atom coexistence ---

  @rq-b604feb4
  Scenario: A ghost atom with a non-empty label still emits "<symbol>-Bq" (label dropped)
    Given an AtomicInput whose molecule has one atom with symbol "H", atom_label "ghost_a", real=False
    When build_input() is called
    Then the atom block line for that atom begins with "H-Bq "
    And the atom block line does not contain "ghost_a"

  # --- Harvester: round-trip with labeled molecule ---

  @rq-5c529f02
  Scenario: out_mol preserves atom_labels when fix_com and fix_orientation are True
    Given an in_mol with atom_labels ["", "5", "_other", "_4sq"]
    And in_mol has fix_com=True and fix_orientation=True
    When harvest(in_mol, "hf", log_content) is called
    Then list(out_mol.atom_labels) equals ["", "5", "_other", "_4sq"]

  # --- Integration test: test_atom_labels ---

  @rq-3dd2d934
  Scenario: test_atom_labels runs successfully with Gaussian
    Given the labeled-hydrogen molecule from test_atom_labels (4 H atoms, charge=0, multiplicity=1)
    And keywords = {} and method "mp2" with basis "aug-cc-pvdz"
    When qcng.compute is called with program "gaussian" and driver "energy"
    Then the result is successful
    And atres.properties.nuclear_repulsion_energy matches 1.0828427 (atol 1e-4)
    And atres.properties.calcinfo_nmo equals 36

  @rq-351864ac
  Scenario: Gaussian SCF and MP2 on the labeled-H system match one of the two accepted solutions
    Given the labeled-hydrogen molecule and method "mp2" with basis "aug-cc-pvdz"
    When qcng.compute is called with program "gaussian" and driver "energy"
    Then atres.properties.scf_total_energy matches either -1.656138508 or -1.705577613 (atol 3e-6)
    And atres.return_result matches either -1.7926264513 or -1.870251459939 (atol 3e-6)
    And the chosen pair is consistent (primary scf with primary mp2, OR alternative scf with alternative mp2)

  @rq-b074f046
  Scenario: The "alternative SCF solution" branch in test_atom_labels accepts Gaussian
    Given qcprog is "gaussian"
    And Gaussian converged to the alternative SCF solution
    When the test reaches the try/except AssertionError branch
    Then the alternative-solution branch is executed (no AssertionError propagates)
```

## Method and Driver Compatibility <!-- rq-c3ebee03 -->

This feature places no method/driver restrictions. The Gaussian harness already
supports MP2 (cf. `input-builder.md` §`muster_modelchem`); the labels are
metadata that bypass the input builder entirely.

## Out-of-Scope <!-- rq-c2a6aabb -->

The following are explicitly **not** part of this feature and would require a
separate plan:

- **Per-atom basis assignment via labels.** Gaussian supports per-type basis
  overrides via `Gen` basis specifications. This would couple the harness to the
  `model.basis` type (string vs BasisSet) in a much deeper way than the current
  string-only contract.
- **Per-atom isotope or spin specifications.** Gaussian's `<element>(Iso=N)` and
  `<element>(Spin=N)` modifiers could be plumbed from `mol.atom_labels` or from
  separate QCSchema fields, but doing so requires a dedicated label-format
  contract that no QCEngine harness currently provides.
- **Echo-back of atom labels in the parsed output.** cclib does not parse atom
  labels from Gaussian output; reconstructing them in `out_mol` when
  `fix_com=False`/`fix_orientation=False` would require additional regex parsing
  of Gaussian's "Symbolic Z-matrix" section. Out of scope.

## Migration Notes <!-- rq-1601c435 -->

- No harness source files are modified by this feature.
- No requirements files other than this one are modified by this feature.
- The `test_atom_labels` test's `except AssertionError` branch (currently gated
  on `qcprog == "nwchem"`) is widened to include `"gaussian"` as part of the
  implementation. This is a one-line edit; see the Gherkin scenario tagged
  "alternative SCF solution branch."
