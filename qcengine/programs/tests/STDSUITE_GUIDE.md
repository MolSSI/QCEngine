# Standard Suite ("stdsuite") Guide

This guide is for agents and humans extending QCEngine's standard-suite tests: new methods, new routes
through a program's modules, new programs, or new harness capabilities. Read it before editing the files
below. It records conventions that aren't obvious from the code.

psi4 keeps its own stdsuite (psi4 `tests/pytests/test_standard_suite.py`) that shares this directory's
reference data and contracts; see psi4's guide for that copy. The QCEngine copy was greatly extended by
porting qcdb's stdsuite (qcdb is archived), so its row conventions follow qcdb's.

## 1. What the stdsuite tests

For every method x driver (energy/gradient/Hessian) x reference (RHF/UHF/ROHF) x frozen core (AE/FC) x
program x module route, a test checks that the program, run through QCEngine,

* runs, or raises the *expected* error for combinations the program can't do;
* routes to the expected module (`AtomicResult.provenance.module`);
* produces the right numbers: return energy/gradient/Hessian, and the full set of QCVariables an ideal
  implementation of the method should produce (the "contract"), each in `AtomicResult.extras["qcvars"]`
  and, where `AtomicProperties` has a field for it, in `AtomicResult.properties`;
* agrees with conventional-integral reference values (to a looser tolerance) when it runs DF/CD;
* returns sensible generics (nbasis, nmo, nalpha, nbeta) and the `CURRENT *` aliases.

The stdsuite runs QCSchema v2 only (`schema_versions_v2` fixture). Schema v1 and v1<->v2 conversion are
covered by other tests that use the `schema_versions` fixture.

## 2. File map

| file | role |
|---|---|
| `test_standard_suite.py` | ~60 `test_<method>_<driver>_module` functions (plus `test_ccsd_energy_default`); all parametrization, error/wrong tuples, convergence dicts |
| `standard_suite_runner.py` | `runner_asserter()`: derives scf/corl types, picks reference blocks and tolerances, runs QCEngine, asserts |
| `standard_suite_contracts.py` | `contractual_<method>()` generators: which QCVariables are expected, with per-program/module excusals |
| `standard_suite_ref.py` | `std_molecules`, `std_refs`, reference blocks `std_suite` keyed by `answer_hash` |
| `test_alignment.py` | HF tests with scrambled (shifted/rotated/atom-reordered) input and fixed/free frames; reuses `_processor` |
| `../qcvar_identities_resources.py` | identities that fill in derived QCVariables (totals, doubles, MP2.5, corrections, ...) in every harness's `build_out` |
| `../<program>/germinate.py` | `muster_modelchem()`: method (+ module) -> program input. New methods for a program start here |
| `../<program>/harvester.py` | output parsing -> QCVariables, and the `module` that ran |
| `pyproject.toml` | marker registration; `qcengine/testing.py` `_programs` for `using("<program>")` |

## 3. Reference data model (`standard_suite_ref.py`)

* `std_molecules`: three-molecule sets, one per reference. Each `basis, subjects` param row lists
  `[rhf_mol, uhf_mol, rohf_mol]`: `["hf", "bh3p", "bh3p"]` for cc-pvdz, `["h2o", "nh2", "nh2"]` for
  aug-cc-pvdz and cfour-qz2p.
* Blocks are keyed by `answer_hash(system, basis, scf_type, reference, fcae, corl_type, sdsc)`. The runner
  derives `corl_type` from the method and `*_type` keywords (psi4 defaults are hard-coded there) and
  `scf_type` from keywords or the natural pairing (conv <-> pk).
* `sdsc`: standard (`sd`) vs semicanonical (`sc`) orbitals; matters only for ROHF. ROHF defaults to `sc`
  unless the row sets `"sdsc": "sd"`. Programs differ (e.g., CFOUR VCC ROHF-CCSD is `sc`, CFOUR with
  `orbitals: 0` is `sd`, psi4 ccenergy ROHF-CCSD is `sd`).
* Only primary quantities are stored; derived ones come from `compute_derived_qcvars` (refs) and
  `qcvar_identities` (harness output). If a contract expects a QCVariable that a program doesn't print but
  that follows from ones it does, add the identity to `qcvar_identities_resources.py`, not the harvester.

## 4. How one test executes

```
test_<mtd>_<drv>_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2)
  -> _processor(): picks the molecule by reference; resolves per-dertype or per-basis "error"/"wrong";
                   per-dertype "xptd" values and "marks"; sets model.basis (program-specific spelling via
                   _trans_key), driver, scf_type="pk", corl_type="conv", and qc_module_in
                   ("<prog>-<keywords.qc_module|cc_program|mp2__code>")
  -> runner_asserter():
       * skips (e.g., cfour-qz2p basis outside CFOUR/psi4)
       * corl_type/scf_type -> chash (exact block) and chash_conv (pk/conv block)
       * tolerances: atol_e/g/h 2e-7/5e-7/1e-5; DFT gradients 6e-6; "xptd": {"fd": True} rows (3-point
         internal findif) 1e-4 (dz) or 2e-5 gradients; atol_conv 1e-4 vs CONV refs
       * qcng.compute(...)
       * "error" rows: assert exception type + regex, then return (test passes)
       * routing: provenance module vs qc_module requested (keywords) and expected ("xptd": {"qc_module": ...})
       * qcvar_assertions(): for each contract QCVariable, compare to exact block and to CONV block, in
         extras.qcvars and (if mapped) properties; excused ones must be *absent*
       * "wrong" rows: expect an AssertionError naming the QCVariable, then xfail
       * CURRENT aliases, return_result/properties, generics
```

## 5. Param rows

```python
pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"fd": True}, "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"}, "error": {...}}, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
```

* `call`: QCEngine program name. `reference`, `fcae`: select the molecule and reference block; they don't set
  keywords. Keywords must set the reference and frozen core *natively* (`reference: uhf` for CFOUR/psi4,
  `scf__uhf: True` for NWChem; `dropmo: [1]` for CFOUR, `freeze_core: True` for psi4, `<module>__freeze: 1`
  for NWChem).
* `keywords`: **native, unprefixed** program keywords, exactly as a user would pass them to QCEngine.
* Convergence: rows spread a shared dict so references reproduce: `_c4_tight` (CFOUR), `_p4_conv` (psi4),
  `_nwc_conv` (NWChem). `_p4_fd` tightens psi4's own findif (5-point, no projection) for rows where psi4 does
  findif internally. DFT rows also spread grid dicts (`_psi_grid`, `_nwc_grid`).
* Module routing:
  * psi4: `qc_module` keyword (`occ`, `fnocc`, `ccenergy`, `detci`, `mrcc`, ...).
  * CFOUR: `cc_program` keyword (`vcc`, `ecc`, `ncc`). CFOUR's own default is VCC.
  * NWChem: `qc_module` is a module *switch*: `True`/`"tce"` -> TCE; `"directmp2"` -> semidirect MP2;
    absent -> the default module (`mp2grad` for MP2, `cc` for CCSD, ...). Document the default module with
    `"xptd": {"qc_module": "mp2grad"}` rather than a keyword.
* `"xptd"` (expectations, not inputs): `qc_module` (routing check), `fd` (`False`: internal findif known to
  be accurate; `True`: internal 3-point findif, loosen tolerances). Values may be per-dertype dicts,
  e.g., `"xptd": {"fd": {1: True}}`.
* `"sdsc"`: `"sd"`/`"sc"` to select the ROHF reference block.
* `"error": {dertype_or_basis: _qNN}`: `_qNN = (ExceptionType, regex, human_reason)`. The program (or its
  QCEngine harness) can't do this combination; the test passes if that exception is raised.
* `"wrong": {dertype_or_basis: _wNN}`: `_wNN = (QCVariable, reason)`. The program runs but gives a
  nonstandard answer; the test xfails at that QCVariable.
* `"marks": {dertype: pytest.mark.xfail(...)}`: runtime marks. Selection marks (`long`, program marks) must
  go on the `pytest.param` instead, since they're needed at collection.
* `marks=using("<program>")` on every row; rows needing two codes use `marks=[*using("psi4"), *using("mrcc")]`.
  Rows needing a minimum program version add a versioned entry, e.g., `*using("psi4_112")` (see `_programs` in
  `qcengine/testing.py`, where a not-yet-released version may be short-circuited to `True`).
* `@pytest.mark.nonroutine` (on a test function, row, or param): skipped unless `--runnonroutine` is given
  (option defined in `qcengine/conftest.py`). All method tests except HF, MP2 energy & gradient, CCSD energy, and B3LYP carry it on
  the function, so a default run stays near its pre-qcdb-port length. Mark new method tests nonroutine unless
  they're core.
* Commented `# DEBUG` rows are alternate routes kept for investigation, not regular coverage.

Error handling categories (also in the file header): E1 can't run -> `"error"`; E2 runs but some QCVariables
not harvestable -> excuse in the contract; E3 parameter not available (e.g., basis) -> skip in runner;
E4 runs but wrong -> `"wrong"`.

## 6. Running

```bash
pytest qcengine/programs/tests/test_standard_suite.py -n 30                                # routine: HF, MP2 ene/grd, CCSD ene, B3LYP
pytest qcengine/programs/tests/test_standard_suite.py -m "cfour or nwchem" --runnonroutine -n 30  # all methods; ~25 min on 30 cores
pytest qcengine/programs/tests/test_standard_suite.py -m psi4 -k "dz and not adz" --runnonroutine -n 30  # psi4 subset; full psi4 is hours
pytest qcengine/programs/tests/test_standard_suite.py -k "test_ccsd_energy and nwchem" -m "not long"
pytest qcengine/programs/tests/test_standard_suite.py --runnonroutine -k test_fci_          # a nonroutine method
```

* Test ids contain spaces and parentheses. To rerun a list of failures, write node ids to a file, one per
  line, and run `pytest @failed.txt` (pytest >= 8.2).
* With `-n`, set `COLUMNS=400` and use `-rf --tb=no` to get the full assertion message per failure; the
  messages name the QCVariable, the store (`extras.qcvars` vs `properties.<field>`), and whether the miss is
  vs exact ref or vs CONV.

## 7. Workflows

### 7a. Adding a row for an existing method and program
Copy the nearest row; set reference/frozen core natively; spread the program's convergence dict; run all
three bases. If the program can't do it, find the actual exception text and add an `_qNN`.

### 7b. Adding a method to a program (harness work)
1. `germinate.py`: map the method to the program's input (CALC_LEVEL, TCE flag, task line). Raise
   `InputError` for unsupported method/module combinations rather than letting the harness fail later.
2. `harvester.py`: parse the method's energies and set `module`. Check the routing assertion passes.
3. `qcvar_identities_resources.py`: add identities for derived QCVariables.
4. `standard_suite_contracts.py`: `contractual_<method>`; excusals keyed precisely on `qc_module`/
   reference/method, never in the runner.
5. Runner: `corl_natural_values` and a branch in `qcvar_assertions` if the method is new to the stdsuite.

### 7c. Adding a program
Add `using("<prog>")` support (`qcengine/testing.py` `_programs` and a marker in `pyproject.toml`), the basis
spelling in `_trans_key`, and rows for the six rhf/uhf/rohf x ae/fc combinations of each method.

## 8. Triage: failure message -> likely meaning

| symptom | first suspects |
|---|---|
| `X expected but missing from extras.qcvars` | harvester doesn't parse it, or a `qcvar_identities` rule is missing, or the contract shouldn't expect it for this module |
| `X expected but missing from properties.<field>` | the harness's `build_atomicproperties` mapping, or the program's own QCSchema output omits it |
| `<var> wrongly present` | contract excuses it but the program now provides it: update the contract |
| `QC_MODULE used (prog) != requested (prog-mod)` | harvester didn't detect the module for this output |
| `QC_MODULE used (prog-a) != expected (prog-b)` | default routing differs from the row's `"xptd": {"qc_module": ...}` |
| `Not found: <Exc> '<regex>'` | the program's error message changed; update the `_qNN` regex |
| `DID NOT RAISE` | the program gained the capability (or now falls back, e.g., psi4 findif): convert the row to a passing one |
| gradient off ~1e-5 for internal findif | 3-point findif: mark `"xptd": {"fd": True}`, or for psi4 spread `_p4_fd` |
| `KeyError` on `std_suite[...]` | no reference block for the derived hash (check scf/corl type and sdsc) |
| NWChem TCE `non-Abelian symmetry not permitted` on cc-pvdz | linear HF molecule; rows carry `"error": {"cc-pvdz": _q45}` |
