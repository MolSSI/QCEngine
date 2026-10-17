# QCEngine QCSchema examples

These directories contain JSON representations of QCSchema instances created
by the QCEngine test suite. The model family is encoded in the path:

```
v1/<Model>/qcengine-<test>.json
v2/<Model>/qcengine-<test>.json
```

The initial corpus comes from the Psi4 cases in the HF and CCSD(T) standard
suite tests and contains both the engine-facing `AtomicInput` and returned
`AtomicResult` for each case.

Generate and then validate the corpus with two separate commands:

```
pytest --qcschema-examples qcengine/programs/tests/test_standard_suite_hf.py 'qcengine/programs/tests/test_standard_suite_ccsd(t).py'
pytest --validate-qcschema-examples
```

Only the direct `as_v1` and `as_v2` schema pathways are written. Generation
removes stale JSON first, and ordinary test runs do not write example files.
Stored results omit program stdout and runner-specific host, user, resource,
and timing metadata. Result floats are stored to 12 significant digits to
remove insignificant last-bit variation from threaded engines, with absolute
values below 1e-12 normalized to zero.
`manifest.json` records source revisions, package versions, and counts for the
most recently generated corpus.
