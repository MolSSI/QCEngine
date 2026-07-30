# CCLib Harness Review Refactor Design

**Date:** 2026-07-30  
**Branch:** `cclib_harness`  
**Status:** Approved for implementation planning

## Purpose

Address the review of the initial cclib-backed Q-Chem and ORCA harnesses by replacing the monolithic implementation with an extensible package, removing runtime synthetic compatibility checks, passing native methods through to the selected program, consolidating tests, and aligning documented driver support with verified behavior.

## Package architecture

The implementation will use program-specific harness subclasses over a shared base:

```text
qcengine/programs/cclib_programs/
├── __init__.py
├── base.py
├── cclib_qchem.py
├── cclib_orca.py
└── tests/
    └── test_cclib.py
```

The existing `qcengine/programs/cclib.py` and `qcengine/programs/tests/test_cclib.py` will be removed. `qcengine/programs/base.py` will register `QChemCCLibHarness()` and `ORCACCLibHarness()` directly.

### Shared base

`cclib_programs/base.py` will contain:

- the lazy cclib API loader;
- shared job and execution-result data structures;
- a documented program-definition interface;
- managed execution and bounded failure reporting;
- cclib parser invocation and parser-identity checking;
- QCSchema v1 validation and v2 conversion;
- native-file protocol handling;
- the shared `CCLibHarness` implementation.

The shared base will not contain native Q-Chem or ORCA input syntax, executable probes, environment preflight logic, or program-name branches. A program-specific harness subclass will provide one immutable program definition consumed by the inherited harness methods.

The lazy cclib API will expose only the interfaces required for real parsed results, including cclib's version, auto-detecting opener, and `QCSchemaWriter`. It will not expose `ccData`. Each program module will lazily provide its expected parser class so that importing QCEngine does not import cclib.

### Q-Chem module

`cclib_programs/cclib_qchem.py` will own:

- Q-Chem keyword validation and deterministic `$rem` rendering;
- Q-Chem input generation;
- Q-Chem environment preflight;
- executable identity/version probing;
- output selection and termination identity;
- `QChemCCLibHarness`.

Q-Chem will support energy, gradient, and Hessian drivers. Native method strings will be emitted without a QCEngine method allowlist.

### ORCA module

`cclib_programs/cclib_orca.py` will own:

- ORCA simple-keyword and block validation;
- deterministic ORCA input generation;
- executable identity/version probing;
- stdout selection and termination identity;
- `ORCACCLibHarness`.

ORCA will support energy only. Live verification with the available ORCA 6.1.1 and development cclib checkout showed that gradient and Hessian calculations execute but fail during cclib parsing, so those drivers will be rejected before execution until the parser path supports them.

## Compatibility and availability

The runtime synthetic `ccData` object, writer compatibility probe, compatibility cache, and associated global state will be removed. `QCSchemaWriter` is responsible for accepting real parser output; QCEngine will validate its real output during conversion rather than fabricate parser data during availability checks.

Availability will:

1. lazily import the cclib opener and writer;
2. resolve the selected executable through `PATH`;
3. run the program-specific preflight;
4. probe the real executable for identity and version;
5. enforce the selected program's minimum version.

Executable versions remain cached by resolved executable path, matching established QCEngine harness practice. Unit tests will not fabricate executable banners to test version parsing. Live, availability-gated tests will call `get_version()` against installed software.

## Input and result behavior

The harness will retain validation for constraints it owns:

- the requested driver must be supported by that cclib program module;
- basis must be a non-empty native string;
- all atoms must be real;
- native keyword structures must be safe and valid for generated input;
- harness-owned resource and geometry fields cannot be overridden.

The harness will not restrict methods to a QCEngine-maintained list and will not maintain method aliases. It will pass the requested method spelling to the downstream program, which decides whether it is supported.

Parser identity, parsed driver, and parsed basis will continue to be checked to prevent accepting output from the wrong calculation. Parsed method metadata will be preserved as emitted by cclib rather than rejected through a QCEngine alias table. Existing writer-owned extras collision protection, portable temporary-file handling, stage-aware errors, and protocol behavior remain.

## Documentation

`docs/source/programs_cclib.rst` will be updated to:

- document Q-Chem energy/gradient/Hessian support;
- document ORCA energy-only support and remove unsupported gradient/Hessian claims;
- state that methods pass through to the downstream program;
- describe the `cclib_programs` extension layout and responsibilities required to add another cclib-backed program;
- retain keyword, resource, installation, result, and native-file contracts that remain accurate.

Every shared and program-specific helper will receive a concise docstring explaining its responsibility. The module and class docstrings will describe how future cclib-backed program modules integrate with the shared base.

## Test design

Tests will move to:

```text
qcengine/programs/cclib_programs/tests/test_cclib.py
```

The suite will be reduced rather than merely relocated. It will remove:

- synthetic cclib compatibility and cache tests;
- mocked executable-banner/version tests;
- redundant permutations that assert the same implementation detail;
- low-value demonstration serialization internals already covered by live paths.

Focused package tests will retain coverage for:

- registration, distinct harness instances, and lazy optional cclib imports;
- representative Q-Chem and ORCA input rendering;
- program-specific driver acceptance and rejection;
- malformed keyword and reserved-field rejection;
- Q-Chem environment preflight;
- managed execution, primary output selection, scratch behavior, and representative failure classification;
- parser lifecycle and parser-identity rejection;
- real writer conversion, required-field validation, extras collision handling, identity checks, native-file protocols, and public schema conversion;
- cclib fixture parsing where `CCLIB_SOURCE_ROOT` is available;
- the two live demonstrations.

`qcengine/tests/test_harness_canonical.py` will include:

- `cclib-qchem` energy;
- `cclib-qchem` gradient;
- `cclib-orca` energy.

A focused availability-gated package test will exercise Q-Chem Hessian end to end. Availability-gated tests for both selectors will call `get_version()` and enforce their real minimum versions. No machine-local executable paths or setup scripts will appear in reusable source or unit tests.

## Verification

Implementation completion requires:

1. focused package tests;
2. canonical harness tests for both selectors;
3. the non-addon QCEngine suite;
4. cclib fixture conversion tests with `CCLIB_SOURCE_ROOT`;
5. live Q-Chem energy, gradient, and Hessian execution;
6. live ORCA energy execution;
7. both demonstration scripts;
8. Sphinx with warnings as errors;
9. compile checks, `git diff --check`, and an audit for machine-local paths;
10. independent review of the complete branch diff.

## Review-comment disposition

1. Driver claims become end-to-end tested and ORCA's unsupported claims are removed.
2. Production code no longer explicitly depends on `ccData`.
3. Synthetic parser data is removed from production.
4. The synthetic compatibility global cache is removed; only standard executable-path version caching remains.
5. The runtime synthetic compatibility probe is removed.
6. The monolith is replaced by the documented base/Q-Chem/ORCA package split with helper docstrings.
7. Method allowlists and aliases are removed; native programs decide method support.
8. Tests are consolidated, reduced, relocated, and both selectors join canonical harness tests.
9. The synthetic missing-cclib compatibility-probe test is removed with that probe.
10. Mocked version transcripts are replaced by availability-gated checks of installed software.
