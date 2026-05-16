# Feature: Gaussian Harness — Empirical Dispersion Support <!-- rq-17959824 -->

This document specifies empirical-dispersion support for the Gaussian harness, covering
both the input side (translating method strings like `b3lyp-d3` into Gaussian route
keywords) and the output side (harvesting the `Dispersion energy=` line and emitting
the standard QCEngine dispersion qcvars).

Gaussian's empirical dispersion is invoked by the route-section keyword
`EmpiricalDispersion=<G..>`, where the value selects the dispersion type:

| QCEngine alias | Gaussian keyword value | Formal label | dashcoeff key |
|----------------|------------------------|--------------|---------------|
| `d2`, `d` | `GD2` | `D2` | `d2` |
| `d3`, `d3zero`, `d32b`, `d3zero2b` | `GD3` | `D3` | `d3zero2b` |
| `d3bj`, `d3bj2b` | `GD3BJ` | `D3(BJ)` | `d3bj2b` |

Gaussian does not natively support `D3M`, `D3M(BJ)`, `D3OP`, `D4`, `-NL`, `CHG`,
`DAS2009`, or `DAS2010`. Requests for these levels are rejected by the harness
with `InputError`.

Non-local-correlation DFT functionals (e.g. `wB97M-V`, `B97M-V`, `wB97X-V`) are
**not** empirical-dispersion methods — their NLC contribution is part of the SCF
(stored in cclib's `scfenergies` already) — and are handled by the existing
pass-through DFT branch unchanged.

Related files:
- `input-builder.md` — `germinate.py` (`muster_modelchem` extended by this feature) and
  `keywords.py` (`build_route_line`, no functional change but documented behaviour).
- `harvester.md` — `harvester.py` (energy extraction extended by this feature).
- `runner.md` — `runner.py` (conflict detection between method-suffix and
  user-keyword `EmpiricalDispersion`).
- `ghost-atoms.md`, `atom-labels.md` — unrelated.

## Background: Gaussian Dispersion Output <!-- rq-29751017 -->

When Gaussian runs a calculation with `EmpiricalDispersion=GD3` (or a GD2/GD3BJ
analogue), the log contains a single line of the form:

```
 Dispersion energy=        -0.0046617305 Hartrees.
```

The numeric value is in Hartrees. cclib parses this line into
`ccdata.dispersionenergies` (as a list, in eV), separately from `ccdata.scfenergies`
(which remains the pure SCF without the dispersion contribution).

The `SCF Done:` line that the existing harvester regex (`_SCF_DONE_RE`) matches is
also the **pure SCF** value, *not* the dispersion-corrected total.

For analytic gradients, Gaussian writes a single forces table that already includes
the dispersion contribution; no special handling is required on the gradient side.

## Feature API <!-- rq-7d1fb16e -->

### Extended: `muster_modelchem(method, driver, multiplicity)` <!-- rq-0c1b5066 -->

Located in `qcengine/programs/gaussian/germinate.py`.

**New behaviour**:

- If `method` (lower-cased) ends in `"-<alias>"` for some `alias` in
  `qcengine.programs.empirical_dispersion_resources.get_dispersion_aliases()`:
  - The functional name is the prefix before `"-<alias>"`.
  - The canonical dashlevel key is `get_dispersion_aliases()[alias]`.
  - The canonical key is looked up in the harness's `_GAUSSIAN_DISPERSION_MAP`
    (defined in `germinate.py`):
    ```python
    _GAUSSIAN_DISPERSION_MAP = {
        "d2":        ("GD2",   "D2"),
        "d3zero2b":  ("GD3",   "D3"),
        "d3bj2b":    ("GD3BJ", "D3(BJ)"),
    }
    ```
    A canonical key not in this map raises `InputError` with a message naming the
    requested level and the levels Gaussian supports.
  - The returned dict gains two keys:
    - `"dispersion_level"`: the formal label (e.g. `"D3"`, `"D3(BJ)"`, `"D2"`).
    - `"functional_name"`: the upper-cased functional prefix (e.g. `"B3LYP"`).
  - The Gaussian `EmpiricalDispersion=<G..>` token is added to the returned
    `extra_keywords` dict.
  - The returned `method_string` is the functional alone (e.g. `"B3LYP"`), **not**
    including the dispersion suffix, because Gaussian's `<method>/<basis>` token
    must name only the functional.

- If `method` has no recognised dispersion suffix, all existing behaviour is
  preserved (no `EmpiricalDispersion` added; `dispersion_level` and
  `functional_name` keys are absent from the returned dict).

**Open-shell handling**: the U prefix (for multiplicity > 1) is applied to the
**functional** before the dispersion suffix is stripped/added back. E.g.
`method="b3lyp-d3"` with multiplicity 2 produces `method_string="UB3LYP"` and
`extra_keywords={"EmpiricalDispersion": "GD3"}`.

### Extended: `build_input(input_model, config, template=None)` <!-- rq-59b816b7 -->

Located in `qcengine/programs/gaussian/runner.py`.

**Conflict check**: if `muster_modelchem` returned a non-None `dispersion_level`
AND `input_model.specification.keywords` contains an `EmpiricalDispersion` key
(case-insensitive), raise `InputError` with a message indicating that the
dispersion type is overspecified.

**Pass-through**: if `muster_modelchem` did **not** detect a dispersion suffix but
the user explicitly passed `keywords["EmpiricalDispersion"] = "GD3"` (or similar),
the harness emits the route line as before with the keyword pass-through; the
harvester (see below) detects dispersion from the log itself, so qcvars are still
emitted correctly.

### Extended: `harvest(in_mol, method, log_content)` <!-- rq-1d5da737 -->

Located in `qcengine/programs/gaussian/harvester.py`.

**Method normalisation**: before energy extraction, strip any recognised dispersion
suffix from the (lowercased) method. The resulting bare method drives the existing
energy-extraction branch. The original method (with suffix) is retained so qcvars
can be tagged with the correct formal label.

**Dispersion-energy parsing**: after the existing SCF and post-HF energy
extraction:

- A new regex `_DISPERSION_ENERGY_RE = re.compile(r"Dispersion energy=\s*([-\d.]+)\s*Hartrees?")`
  is matched against the log. If a match exists, `disp_hartree = float(match.group(1))`.
- If the regex does not match, the cclib fallback is
  `disp_hartree = cclib.parser.utils.convertor(float(ccdata.dispersionenergies[-1]), "eV", "hartree")`
  when the attribute exists and is non-empty. The cclib converter is used
  (rather than `* qcel.constants.conversion_factor("eV", "hartree")`) so the
  eV → Hartree round-trip is exact relative to how cclib stored the value
  — see [[harvester]] §"cclib Limitations" for the bias-vs-exactness analysis.
- If neither source yields a value, `disp_hartree` is `None` and dispersion qcvars
  are not set.

**qcvars emitted when `disp_hartree is not None`** (Psi4/QCEngine convention; see
`qcvar_identities_resources.py`):

| Variable | Value |
|----------|-------|
| `"DISPERSION CORRECTION ENERGY"` | `disp_hartree` |
| `"<FCTL>-<FORMAL> DISPERSION CORRECTION ENERGY"` | `disp_hartree` |
| `"<FCTL>-<FORMAL> TOTAL ENERGY"` | `scf_hartree + disp_hartree` |
| `"<FCTL> FUNCTIONAL TOTAL ENERGY"` | `scf_hartree` (alias for clarity) |
| `"CURRENT ENERGY"` (DFT branch only) | `scf_hartree + disp_hartree` |

Where `<FCTL>` is the upper-cased functional name and `<FORMAL>` is the formal
dash label.

**HF TOTAL ENERGY / SCF TOTAL ENERGY semantics** are unchanged: both refer to the
pure SCF without dispersion. This is consistent with the QCEngine/Psi4 convention
in which `"SCF TOTAL ENERGY"` excludes the empirical dispersion correction.

**Gradient**: unchanged. Gaussian's analytic gradient already includes dispersion
contributions; cclib's `grads` and the existing parsing path produce the correct
total gradient.

**Hessian**: same — cclib's `hessian` is the total. (Note: Gaussian supports
analytic Hessians with D2/D3 but with restrictions for some method/basis combos;
this is enforced at runtime by Gaussian and surfaced via the existing log-pattern
dispatch.)

## Method-Name Recognition Details <!-- rq-f357f8db -->

The recognition step in `muster_modelchem` is:

```python
from qcengine.programs.empirical_dispersion_resources import get_dispersion_aliases, dashcoeff

method_lower = method.lower()
dispersion_level = None
functional_name = None
emp_disp_value = None

aliases = get_dispersion_aliases()
# Sort by descending length to prefer e.g. "d3bj" over "d3" + "bj"
for alias in sorted(aliases, key=len, reverse=True):
    suffix = f"-{alias}"
    if method_lower.endswith(suffix):
        canonical_key = aliases[alias]
        if canonical_key not in _GAUSSIAN_DISPERSION_MAP:
            raise InputError(
                f"Dispersion level '{alias}' (canonical: {canonical_key}) is not "
                f"natively supported by Gaussian. Supported levels: "
                f"{sorted({fmt for _, fmt in _GAUSSIAN_DISPERSION_MAP.values()})}."
            )
        emp_disp_value, dispersion_level = _GAUSSIAN_DISPERSION_MAP[canonical_key]
        functional_name = method[: -len(suffix)]
        break
```

The "sort by descending length" rule ensures `b3lyp-d3bj` matches the `d3bj` alias
(not `d3` followed by a stray `bj` mismatch). For two aliases of the same length
that both match a tail, ties are broken by the first one encountered in the
sorted iteration (consistent across runs because the alias dict is stable).

## Open Cases <!-- rq-f237caf7 -->

- Pure HF with dispersion (e.g. `hf-d3`): supported. The functional name is `HF`
  and the route line is `#P HF/<basis> EmpiricalDispersion=GD3`.
- Methods named with embedded hyphens that aren't dispersion suffixes (e.g.
  hypothetical `cam-b3lyp`): the recognition step only strips a *trailing* hyphen
  segment, and `b3lyp` is not in the alias map, so `cam-b3lyp` is left intact. No
  false stripping.
- Repeat tokens (e.g. `b3lyp-d3-d3`): the trailing `-d3` is stripped, leaving
  `b3lyp-d3` as the functional. Gaussian will see `#P B3LYP-D3/...` which it does
  not recognise, and will fail; the existing log-pattern dispatch surfaces this
  as an `InputError`. We do not pre-validate.

## Gherkin Scenarios <!-- rq-dea9f7c9 -->

```gherkin
Feature: Gaussian harness empirical-dispersion support

  Background:
    Given a TaskConfig with ncores=1, memory=1.0, scratch_directory=None

  # --- germinate.muster_modelchem: method recognition ---

  @rq-1bdc1372
  Scenario: muster_modelchem strips "-d3" and emits EmpiricalDispersion=GD3
    When muster_modelchem("b3lyp-d3", "energy", 1) is called
    Then method_string equals "B3LYP"
    And extra_keywords contains the entry "EmpiricalDispersion" -> "GD3"
    And the returned dispersion_level equals "D3"
    And the returned functional_name equals "B3LYP"

  @rq-6a93ae95
  Scenario: muster_modelchem strips "-d3bj" and emits EmpiricalDispersion=GD3BJ
    When muster_modelchem("b3lyp-d3bj", "energy", 1) is called
    Then method_string equals "B3LYP"
    And extra_keywords contains the entry "EmpiricalDispersion" -> "GD3BJ"
    And dispersion_level equals "D3(BJ)"

  @rq-af4f9bf0
  Scenario: muster_modelchem strips "-d2" and emits EmpiricalDispersion=GD2
    When muster_modelchem("b3lyp-d2", "energy", 1) is called
    Then method_string equals "B3LYP"
    And extra_keywords["EmpiricalDispersion"] equals "GD2"
    And dispersion_level equals "D2"

  @rq-7c3a72b9
  Scenario: muster_modelchem recognises canonical key alias "d3zero2b"
    When muster_modelchem("b3lyp-d3zero2b", "energy", 1) is called
    Then dispersion_level equals "D3"
    And functional_name equals "B3LYP"

  @rq-95aa6edf
  Scenario: muster_modelchem recognises the bare alias "d" as D2
    When muster_modelchem("b3lyp-d", "energy", 1) is called
    Then dispersion_level equals "D2"

  @rq-f4cbb515
  Scenario: muster_modelchem leaves a non-dispersion method unchanged
    When muster_modelchem("b3lyp", "energy", 1) is called
    Then method_string equals "B3LYP"
    And extra_keywords does not contain the key "EmpiricalDispersion"
    And the returned dispersion_level is absent (or None)

  @rq-4262d5f2
  Scenario: muster_modelchem leaves a hyphenated non-dispersion method unchanged
    When muster_modelchem("cam-b3lyp", "energy", 1) is called
    Then method_string equals "CAM-B3LYP"
    And no "EmpiricalDispersion" entry is added
    And dispersion_level is absent

  @rq-11602838
  Scenario: HF with dispersion produces an HF functional
    When muster_modelchem("hf-d3", "energy", 1) is called
    Then method_string equals "HF"
    And extra_keywords["EmpiricalDispersion"] equals "GD3"

  # --- germinate.muster_modelchem: open-shell + dispersion ---

  @rq-879110d8
  Scenario: Doublet adds U prefix to functional before stripping dispersion
    When muster_modelchem("b3lyp-d3", "energy", 2) is called
    Then method_string equals "UB3LYP"
    And extra_keywords["EmpiricalDispersion"] equals "GD3"
    And dispersion_level equals "D3"
    And functional_name equals "B3LYP"

  @rq-e1150329
  Scenario: ROHF + dispersion does not double-prefix U
    When muster_modelchem("rohf-d3", "energy", 2) is called
    Then method_string equals "ROHF"
    And extra_keywords["EmpiricalDispersion"] equals "GD3"

  # --- germinate.muster_modelchem: unsupported levels ---

  @rq-d96c128f
  Scenario: muster_modelchem raises InputError for D3M
    When muster_modelchem("b3lyp-d3m", "energy", 1) is called
    Then InputError is raised
    And the error message mentions Gaussian's supported dispersion levels

  @rq-95535d13
  Scenario: muster_modelchem raises InputError for D3M(BJ)
    When muster_modelchem("b3lyp-d3mbj", "energy", 1) is called
    Then InputError is raised

  @rq-09dd261d
  Scenario: muster_modelchem raises InputError for D4
    When muster_modelchem("b3lyp-d4", "energy", 1) is called
    Then InputError is raised

  @rq-e9a472db
  Scenario: muster_modelchem raises InputError for OP damping
    When muster_modelchem("b3lyp-d3op", "energy", 1) is called
    Then InputError is raised

  @rq-2de0fabe
  Scenario: muster_modelchem raises InputError for -NL non-local dispersion alias
    When muster_modelchem("b3lyp-nl", "energy", 1) is called
    Then InputError is raised

  # --- germinate.muster_modelchem: precedence rules ---

  @rq-ef982ff3
  Scenario: muster_modelchem prefers the longer suffix (d3bj over d3)
    When muster_modelchem("b3lyp-d3bj", "energy", 1) is called
    Then dispersion_level equals "D3(BJ)"
    And functional_name equals "B3LYP"

  # --- runner.build_input: conflict detection ---

  @rq-540caf62
  Scenario: Conflict between method suffix and user keyword raises InputError
    Given an AtomicInput with method "b3lyp-d3" and keywords {"EmpiricalDispersion": "GD3BJ"}
    When build_input() is called
    Then InputError is raised
    And the message indicates that EmpiricalDispersion is overspecified

  @rq-44319e8c
  Scenario: Conflict detection is case-insensitive on the keyword name
    Given an AtomicInput with method "b3lyp-d3" and keywords {"empiricaldispersion": "GD3"}
    When build_input() is called
    Then InputError is raised

  @rq-c68dde2e
  Scenario: User keyword alone (no method suffix) passes through unchanged
    Given an AtomicInput with method "b3lyp" and keywords {"EmpiricalDispersion": "GD3"}
    When build_input() is called
    Then the route line contains "EmpiricalDispersion=GD3"
    And no InputError is raised

  @rq-96507cf6
  Scenario: Method suffix alone (no keyword) produces a single EmpiricalDispersion route token
    Given an AtomicInput with method "b3lyp-d3" and keywords {}
    When build_input() is called
    Then the route line contains exactly one "EmpiricalDispersion=GD3" occurrence
    And the method/basis token is "B3LYP/<basis>", not "B3LYP-D3/<basis>"

  # --- harvester: dispersion energy extraction ---

  @rq-a83b7d5f
  Scenario: Dispersion energy is parsed via regex when the log contains the line
    Given a Gaussian log containing " Dispersion energy=       -0.0046617305 Hartrees."
    When harvest(in_mol, "b3lyp-d3", log_content) is called
    Then qcvars["DISPERSION CORRECTION ENERGY"] equals -0.0046617305 (atol 1e-12)

  @rq-3136523b
  Scenario: Dispersion energy falls back to ccdata.dispersionenergies when regex misses
    Given a Gaussian log without a "Dispersion energy=" line
    And ccdata.dispersionenergies equals [-0.1269192] (in eV)
    When harvest(in_mol, "b3lyp-d3", log_content) is called
    Then qcvars["DISPERSION CORRECTION ENERGY"] equals cclib.convertor(-0.1269192, "eV", "hartree") (atol 1e-8)

  @rq-6f52a6d9
  Scenario: Functional-specific dispersion qcvar uses the formal label
    Given a log with "Dispersion energy=  -0.005 Hartrees." for a B3LYP-D3 calculation
    When harvest(in_mol, "b3lyp-d3", log_content) is called
    Then qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"] equals -0.005

  @rq-f5d05419
  Scenario: Functional total energy excludes dispersion
    Given an SCF energy of -76.0 Ha and a dispersion energy of -0.005 Ha
    When harvest(in_mol, "b3lyp-d3", log_content) is called
    Then qcvars["B3LYP FUNCTIONAL TOTAL ENERGY"] equals -76.0
    And qcvars["B3LYP-D3 TOTAL ENERGY"] equals -76.005
    And qcvars["CURRENT ENERGY"] equals -76.005

  @rq-e89c4144
  Scenario: BJ-damped variant uses parenthesised formal label
    Given a B3LYP-D3(BJ) calculation log with "Dispersion energy=  -0.007 Hartrees."
    When harvest(in_mol, "b3lyp-d3bj", log_content) is called
    Then qcvars["B3LYP-D3(BJ) DISPERSION CORRECTION ENERGY"] equals -0.007
    And qcvars["B3LYP-D3(BJ) TOTAL ENERGY"] equals scf + disp

  @rq-7211fb49
  Scenario: SCF TOTAL ENERGY and HF TOTAL ENERGY exclude dispersion
    Given a B3LYP-D3 calculation
    When harvest is called
    Then qcvars["HF TOTAL ENERGY"] equals the SCF value alone (no dispersion)
    And qcvars["SCF TOTAL ENERGY"] equals qcvars["HF TOTAL ENERGY"]

  @rq-f49de098
  Scenario: Methods without dispersion suffix do not get dispersion qcvars even if log has the line
    Given a Gaussian log containing a "Dispersion energy=" line (e.g. user passed EmpiricalDispersion keyword)
    And the method passed to harvest is "b3lyp" (no dispersion suffix)
    When harvest(in_mol, "b3lyp", log_content) is called
    Then qcvars["DISPERSION CORRECTION ENERGY"] is still populated
    And qcvars["CURRENT ENERGY"] equals scf + disp
    But no functional-specific "<FCTL>-<DASH>" qcvars are set (we cannot infer the dash level)

  @rq-df49b3b8
  Scenario: No dispersion in log produces no dispersion qcvars
    Given a Gaussian log without a "Dispersion energy=" line
    And ccdata has no dispersionenergies attribute (or it's empty)
    When harvest(in_mol, "b3lyp", log_content) is called
    Then no DISPERSION* qcvars are set
    And qcvars["CURRENT ENERGY"] equals qcvars["SCF TOTAL ENERGY"]

  @rq-93777ba8
  Scenario: NLC method (wB97M-V) does not emit dispersion qcvars
    Given a wB97M-V calculation log (NLC self-consistent, no "Dispersion energy=" line)
    When harvest(in_mol, "wb97m-v", log_content) is called
    Then no DISPERSION qcvars are set
    And the existing pass-through DFT behaviour is preserved

  @rq-68840f83
  Scenario: Gradient is unchanged when dispersion is present
    Given a B3LYP-D3 gradient calculation
    When harvest is called with driver "gradient"
    Then the returned gradient is the full gradient (functional + dispersion)
    And no separate "dispersion gradient" qcvar is exposed

  # --- Integration test ---

  @rq-0fe273a8
  Scenario: End-to-end B3LYP-D3/6-31G* energy on water with real Gaussian
    Given an AtomicInput for water (H2O) with method "b3lyp-d3", basis "6-31g*", driver "energy"
    When qcng.compute is called with program "gaussian"
    Then the result is successful
    And qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"] is approximately equal to the empirical reference value (atol 1e-6)
    And qcvars["B3LYP-D3 TOTAL ENERGY"] is approximately equal to the empirical reference value (atol 1e-6)
    And qcvars["B3LYP FUNCTIONAL TOTAL ENERGY"] is approximately equal to the empirical reference value (atol 1e-6)
```

## Integration Test <!-- rq-ce8f4843 -->

A real-Gaussian integration test (`@using("gaussian")`) is added to
`qcengine/programs/tests/test_gaussian.py`:

```python
@using("gaussian")
def test_gaussian_b3lyp_d3_water(water):
    """B3LYP-D3/6-31G* energy on water; verify dispersion qcvars are populated."""
    inp = AtomicInput(
        molecule=water,
        specification=AtomicSpecification(
            program="gaussian", driver="energy",
            model={"method": "b3lyp-d3", "basis": "6-31g*"},
        ),
    )
    res = qcng.compute(inp, "gaussian", raise_error=True)
    qcvars = res.extras["qcvars"]
    # Verify the standard dispersion qcvars are present and consistent.
    assert "B3LYP-D3 DISPERSION CORRECTION ENERGY" in qcvars
    assert "B3LYP-D3 TOTAL ENERGY" in qcvars
    assert "B3LYP FUNCTIONAL TOTAL ENERGY" in qcvars
    # SCF + dispersion == total
    scf = float(qcvars["B3LYP FUNCTIONAL TOTAL ENERGY"])
    disp = float(qcvars["B3LYP-D3 DISPERSION CORRECTION ENERGY"])
    tot = float(qcvars["B3LYP-D3 TOTAL ENERGY"])
    assert abs((scf + disp) - tot) < 1.0e-10
```

**Reference values** for the magnitude check on the dispersion contribution must
be determined empirically by running Gaussian once during implementation; they
are then hard-coded into the test with `atol = 1.0e-6`. The reference values
depend on the Gaussian patch level (D3 parameter set may vary), but are stable
within a major Gaussian revision.

## Migration Notes <!-- rq-ea37df39 -->

- No existing tests are altered; all current Gaussian tests use methods without
  empirical dispersion and remain green.
- The `_METHOD_MAP` in `germinate.py` is unchanged. Dispersion-bearing method
  strings flow through the new branch *before* the existing `_METHOD_MAP` lookup,
  so the existing HF/MP2/CCSD entries continue to apply to the stripped
  functional name. (E.g. `mp2-d3` is recognised as MP2 + D3.)
- The `qcvar_identities_resources.py` solver already declares the relationship
  `<FCTL>-<DASH> TOTAL ENERGY = <FCTL> FUNCTIONAL TOTAL ENERGY + <FCTL>-<DASH> DISPERSION CORRECTION ENERGY`
  for `B3LYP`, `B3LYP5`, `PBE`, `B97`, `BLYP`, `BP86`, `PBE0`, `WPBE` × `-D2`,
  `-D3`, `-D3(BJ)`, `-D3M`, `-D3M(BJ)`. Functionals outside this list still get
  the qcvars emitted but the solver won't auto-derive missing siblings.
- Public API surface is unchanged. No new functions or classes are introduced;
  the changes are confined to extending the existing entry points.

## Out-of-Scope <!-- rq-06ce6129 -->

- **D4 support** — requires a separate executable (`dftd4`) and the existing
  `dftd_ng.py` harness handles it. Gaussian itself does not invoke D4.
- **D3M / D3OP** — Gaussian does not implement these dispersion-damping
  variants. Users can run the SCF in Gaussian and apply the dispersion
  correction post-hoc via the dedicated `dftd3`/`mp2d` harness.
- **NLC-functional plumbing** — wB97M-V et al. are pass-through DFT methods
  that already work with the existing harness; no dispersion qcvars are
  emitted because there is no empirical dispersion contribution.
- **Counterpoise + dispersion** — Gaussian supports `Counterpoise=N` together
  with `EmpiricalDispersion`, but the dispersion contribution behaves differently
  for ghost atoms (depending on the dispersion model). Out of scope for this
  feature; the existing ghost-atom support remains valid for non-dispersion
  calculations.
- **Hessian with dispersion** — Gaussian's analytic Hessian supports some
  dispersion levels but not all. We do not pre-validate; runtime failures
  surface via the existing log-pattern dispatch.

## References <!-- rq-29bb67d8 -->

- Grimme et al., J. Chem. Phys. 132, 154104 (2010) — D3 paper.
- Grimme, Ehrlich, Goerigk, J. Comput. Chem. 32, 1456 (2011) — D3(BJ) paper.
- QCEngine `empirical_dispersion_resources.py` — the canonical alias table and
  parameter database used here.
- Psi4's psivars convention for dispersion — the qcvar naming scheme followed
  by this feature.
