Changelog
=========

.. vX.Y.0 / 2026-MM-DD (Unreleased)
.. --------------------
..
.. Breaking Changes
.. ++++++++++++++++
..
.. New Features
.. ++++++++++++
..
.. Enhancements
.. ++++++++++++
..
.. Bug Fixes
.. +++++++++
..
.. Misc.
.. +++++
..
.. MUST (Unmerged)
.. +++++++++++++++
..
.. WIP (Unmerged)
.. ++++++++++++++
.. - UNMERGED (:pr:`421`) GAMESS - error handling and memory @taylor-a-barnes
.. - UNSOLVED (:issue:`397`) extras failed


v0.50.0rc1 / 2026-01-31 (aka "next" aka "QCSchema v2 available") (Prerelease)
-----------------------------------------------------------------------------

Breaking Changes
++++++++++++++++
- (:pr:`486`) Deps - Require Python >=3.9. Practically for conda, >=3.10.
- (:pr:`453`) Deps - Require Pydantic v2 dependency (don't worry, this isn't
  changing QCEngine's role as QCSchema I/O runner. Also require pydantic-settings
  for CLI. @loriab
- (:pr:`455`) API - As promised 2 years ago for >=v0.30, `local_options` has
  been removed in favor of `task_config` in `compute` and `compute_procedure`.
  Note that Psi4 v1.6 will need an older qcel or a sed to work (see GHA). The
  `qcengine.MDIEngine` is on notice (probably not user-facing. @loriab
- (:pr:`455`) API - `qcengine.compute` and `qcengine.compute_procedure` have been
  merged in favor of the former. Also, treat the second argument (e.g., "nwchem"
  or "geometric") as a positional argument, rather than keyword argument with key
  "program" or "procedure". @loriab
- (:pr:`455`) API - `compute` learned an optional argument `return_version` to
  specify the schema_version of the returned model or dictionary. By default it'll
  return the input schema_version. If not determinable, it will return v1. @loriab

New Features
++++++++++++
- QCSchema v2 implemented and used internally for harnesses.
- Both QCSchema v1 and v2 may be used for input and requested as output.
- Note that QCSchema v2 field layout isn't finalized until v0.60, so there may
  be lock-step advancements of the QCArchive stack.
- Python 3.14 is now useable. Use of QCSchema v1 with 3.14 is limited (by Pydantic).
  Generally, input (as dict or json) will work. v1 output is technically forbidden
  but see envvar below for returning as dict. See table:compute_result_schver in
  compute.py for details. Note that when CMS community codes aren't 3.14-ready,
  their presence may interfere with normal QCEngine functionality (like `qcengine info`).
- See QCElemental for a map of QCSchema v1 and v2, a schema changelog, and a migration guide.

Enhancements (General)
++++++++++++++++++++++
- (:pr:`458`) ``qcengine run`` learned new argument ``--return-version`` analogous
  to ``qcengine.compute(..., return_version=1|2)`` so CLI matches API capabilities.
  Note *not* ported to phasing-out ``qcengine run-procedure``.
- API — Standardize failed job handling toward models for QCSchema v2.
  - (:pr:`458`) When qcengine.compute() fails and forms a fop = FailedOperation
    (raise_error=T), with v2, `fop.input_data` will be an <>Input model (when
    possible; if the error was in forming the model, it'll still be a dict), not
    always a dict like v1.
  - (:pr:`458`) When <executor>.compute() fails and collects the input for
    processing, with v2 it now uses the <>Input model passed to the executor,
    not the model-or-dict passed into compute().
  - (:pr:`458`) The net result of the two above is that whereas fop.input_data
    in v1 was reliably a dict and its contents would reflect whether a model or
    dict was passed to qcengine.compute(), now in v2, fop.input_data is a model
    whenever possible (to mirror <>Result.input_data) regardless of model or dict
    passed to qcengine.compute(); the only case where it's a dict is if the error
    was in forming the model.
  - (:pr:`486`) When the requested type of return can't be formed for Python 3.14,
    it will form the nearest approximation in QCSchema v2.FailedOperation. See
    table:compute_result_schver in compute.py for details.
- (:pr:`454`) Testing - Tests check QCSchema v1 and v2. @loriab
- (:pr:`453`) Maint - Convert internal (non-QCSchema) pydantic classes to
  pydantic v2 API, namely `NodeDescriptor`, `TaskConfig`, `ProgramHarness`,
  `ProcedureHarness`. @loriab
- (:pr:`460`) Testing — Integrate ``AtomicInput.specification`` into harnesses
  (as part of migrating Program harnesses to use QCSchema v2 internally) and
  show what new inputs look like in tests.
- (:pr:`460`) Schema — If you're missing something from AtomicResult.extras (or
  any <>Result.extras), check AtomicResult.input_data.extras in case it was passed
  in on input. Result extras are no longer merged with Input extras.
- (:pr:`468`) Schema — lightly adapt harnesses for Molecule v3 (remember, Molecule
  schema_version is one greater than general QCSchema schema_version. This doesn't
  change the layout, just the version field. Layout not finalized until v0.60.
- (:pr:`486`) Schema — Adapt more harnesses to use ``model_dump`` instead of ``dict``
  when serializing Pydantic models to clear away warnings. Pydantic v1 models use
  the latter and v2 the former but QCSchema aliases each.
- (:pr:`486`) Deps — Adapt harnesses to allow Py 3.14 through small syntax changes.
- (:pr:`486`) API — With :envvar:`QCNG_USE_V1V2_SHIM=1` in Python 3.14, QCSchema
  v1 inputs and outputs (atomic flavor) can be used as dictionaries (never models).
  See the table:compute_result_schver in compute.py for details.

Enhancements (Harnesses)
++++++++++++++++++++++++
- (:pr:`458`) DFTD3 & DFTD4 (new intf) - intercept ``v1.AtomicResult`` with
  ``success=False`` and ``error`` fields set from QCSchema interfaces and return
  ``FailedOperation``s. Someday when upstream switches to v2, request packages
  return FaileOp directly and use ``input_error`` rather than ``input error``.
- (:pr:`459`) OpenMM — Gained AtomicResult.properties.return_gradient
- (:pr:`459`) GCP, MP2D — several got properties.return_energy, retunr_gradient
- (:pr:`461`) Optking — now fills in ``v2.OptimizationResult.stdout``. Through
  QCSchema v2, one can alter gradient protocols in an optimization.
- (:pr:`461`) nwchemdriver — allow nwchemdriver w/o driver=energy. provenance now
  "nwchemdriver" not "nwchemrelax".
- (:pr:`461`) Pyberny — berny harness rewritten in v2. optking and geometric
  natively speak v1, so adapted as well as can be.
- (:pr:`461`) TorsionDrive — torsiondrive rewritten in QCSchema v2.
- (:pr:`468`) Torsiondrive — now accepts protocols. use ``protocols={"scan_results": "all"}``
  if going to be converted to v1.
- (:pr:`468`) RDKIT, MRChem — RDKit, store ``AtomicResult.properties.return_gradient`` and
  ``calcinfo_natom``. MRChem, store ``AtomicResult.properties.return_gradient``.
- (:pr:`468`) TorsionDrive — Adapt harnesses for TD.initial_molecules -> TD.initial_molecule
  and TD.optimization_history -> TD.scan_results
- (:pr:`490`) TorchANI - Updated species handling for newer Pythons, ANI syntax. @loriab
- (:pr:`486`) QCManyBody — Adapt harness to work with QCSchema v1 and v2 and to
  work with QCManyBody v0.50 (in preparation) as well as current v0.5.1.
- (:pr:`486`) Psi4 — Add another generation to the Psi4 harness (to become active
  ~v1.11) to pass/receive QCSchema v2 to/from Psi4.

Bug Fixes
+++++++++
- (:pr:`453`) Maint - Fix a warning thrown by `execute` about unclosed files. @loriab

Misc.
+++++
- (:pr:`452`) Maint - Set up pre-commit and run over repository. @loriab
- (:pr:`453`) CI - Dropped Entos/QCore and Psi4 v1.5 as too hard to solve with
  pydantic v2 and modern python versions. @loriab
- (:pr:`468`) Deps — use packaging instead of setuptools to provide version parsing.
- (:pr:`469`) Maint — setup.py replaced by pyproject.toml with setuptools backend
  and replaced versioneer with setuptools-scm. With removal of versioneer, remove
  __git_revision__ and get_information("git_revision").


v0.34.0 / 2026-01-16
--------------------

Misc.
+++++
- (:pr:`483`) Maint - update CI lanes for smirnoff-Frosst and pkg_resources import changes and Windows challenges. @loriab
- (:pr:`487`) Maint - update CI lanes for spacing and old Psi4 needing defaults channel. @loriab
- (:pr:`489`) Maint - encode in setup.py that QCEngine cannot work with Python 3.14. @loriab


v0.33.0 / 2025-07-31
--------------------

New Features
++++++++++++
- (:pr:`477`) DFTD4 - added a new dispersion type d4bjeeqtwo so parameters can be defined
  for fittings without 3-body dispersion, as for SAPT0-D4 in Psi4. @Awallace3

Misc.
+++++
- (:pr:`476`, :pr:`478`) Maint - repair some CI lanes that needed extra specs to operate correctly. @loriab


v0.32.0 / 2025-05-03
--------------------

Breaking Changes
++++++++++++++++
- (:pr:`472`) Maint - dftd3/gcp/entos/qcore harnesses are all deprecated. See `qcengine info` for details. @loriab

New Features
++++++++++++

Enhancements
++++++++++++
- (:pr:`472`) CLI - Add a `qcengine info bulletin` subcommand (also runs with `qcengine info`)
  for notices about additions and deprecations. @loriab

Bug Fixes
+++++++++
- (:pr:`472`) DFTD4 - New Windows builds of the latest version have a bug with the QCSchema interface
  (https://github.com/dftd4/dftd4/pull/292). QCEngine now compensates so v3.7.0 can be used as-is. @loriab
- (:pr:`472`) geomeTRIC - New v1.1 has been failing one of the mocked "retries" tests here. QCEngine
  now excuses it until upstream fix available (https://github.com/leeping/geomeTRIC/pull/222). @loriab
- (:pr:`472`) Config - Fixed some tests of Config.scratch_directory so they can run cleanly
  on Windows when HOME envvar isn't set. No change to Config itself. @loriab

Misc.
+++++
- (:pr:`472`) Maint - Use new ``packaging`` dependency instead of deprecated ``pkg_resources``. @loriab
- (:pr:`472`) Testing - Drop an unsolvable py38 environment and start testing py313. @loriab


v0.31.0 / 2025-01-17
--------------------

Breaking Changes
++++++++++++++++

New Features
++++++++++++
- (:pr:`466`) QCManyBody - new procedure for computing interaction energies or truncation or full
  many-body expansions with no, counterpoise, or Valiron-Mayer function counterpoise correction
  for basis set superposition error. @loriab

Enhancements
++++++++++++
- (:pr:`463`) Maint - Pin to QCElemental <0.70 since we now know QCSchema v2 release schedule.
- (:pr:`463`) MACE - New v0.3.9 release yields a pytorch error, so recommend pymace=0.3.6 .
- (:pr:`464`, :issue:`447`) CFOUR - Allow CC-PVDZ alias basis specification. Also fix the PwCVXZ
  basis keyword. @philipmnel
- (:pr:`440`) TorsionDrive - Support other geometric-style constraints by not overwriting those
  already present. @jthorton

Bug Fixes
+++++++++
- (:pr:`451`, :issue:`450`) Psi4 - Fixes bug in Psi4 detection when command `python` not available.
  @susilehtola, @topazus
- (:pr:`466`) CFOUR - fix error collecting molecule when it's a single atom with two-letter symbol. @loriab

Misc.
+++++
- (:pr:`457`) Docs - Fix auto-numbering on a documentation page. @mikemhenry


v0.30.0 / 2024-06-25
--------------------

New Features
++++++++++++
- (:pr:`441`) MACE - Added harness for MACE-OFF23 and local MACE models. @jthorton
- (:pr:`443`) AIMNET2 - Added harness for AIMNET2 NN ML models. @jthorton

Misc.
+++++
- (:pr:`445`) CI - fix some test regex issues.
- (:pr:`449`) Maint - bump the QCElemental compatibility range.


v0.29.0 / 2023-10-31
--------------------

Bug Fixes
+++++++++
- (:pr:`427`) Config - Once again, expand environment variables (now more flexibly) and newly expand ``~``
  passed into TaskConfig. Particularly relevant for scratch setting. @yueyericardo
- (:pr:`428`) MDI - Ensure that molecule orientation remains fixed for MDI. @taylor-a-barnes
- (:pr:`405`, :issue:`415`, :pr:`417`) Config - change default ``jobs_per_node`` from 2 to more expected 1
  so a single job fills the node. Alter CPU count formula to return physical cores on Hyperthreading
  machines, affecting default ``ncores``. The net effect (both changes) for default cores running on
  Hyperthreading machines is unchanged. Nevertheless, fixes some Windows problems. @cvsik, @loriab
- (:pr:`433`) Turbomole, Q-Chem - Use raw strings when needed to avoid py312 warnings. @loriab
- (:pr:`435`) GAMESS - Collect the correct MP2 module in parsing for newer versions (>2021,<=2023). @loriab

Misc.
+++++
- (:pr:`433`) CI - Check py312 and some Windows lanes. @loriab


v0.28.1 / 2023-08-18
--------------------

Bug Fixes
+++++++++
- (:pr:`426`) Psi4 - fix ``get_version`` on Windows where whole path and command were getting passed to version parser. @loriab


v0.28.0 / 2023-08-15
--------------------

Breaking Changes
++++++++++++++++

New Features
++++++++++++
- (:pr:`400`) Config - task configuration can now be set via CLI (`qcengine run -h` for details) or
  by environment variables beginning with `QCENGINE_`. @bennybp
- (:pr:`393`, :issue:`392`) MCTC-GCP - Adds b973c and r2scan3c methods to the gcp (mctc only, not classic) harness. @hokru
- (:pr:`393`) DFTD4 - Allows ga, gc, wf parameters to be tweaked (needed for r2scan-3c). This feature requires dftd4 3.5.0. @hokru

Enhancements
++++++++++++
- (:pr:`410`, :issue:`408`) TorsionDrive - silence warnings by using the ``task_config`` argument internally. @jthorton
- (:pr:`409`) Psi4 - improve no-valid-error message so classifies as a RandomError and is eligible for
  restart. @jthorton
- (:pr:`405`) Turbomole - correctly enable OpenMP and environment passing. Pass SCF convergence and
  maximum iterations to define. @cvsik
- (:pr:`403`, :issue:`402`) PyBerny - fix optimizer to respect the task_config options. @q-posev
- (:pr:`386`) CI - turn on formerly LGTM now GitHub CodeQL analysis. @lgtm-migrator
- (:pr:`388`) MRChem - more detailed info about the parallel setup saved to output provenance. @robertodr
- (:pr:`424`) testing - update SVWN Hessian reference values from Psi4. @loriab
- (:pr:`423`, :issue:`377`) NWChem - allow two answers for test ``test_atom_labels[nwchem]`` to accommodate SCF
  solutions in different versions. @loriab

Bug Fixes
+++++++++
- (:pr:`401`) MDI - fix bug in the shape of the MDI forces structure. @taylor-a-barnes
- (:issue:`399`, :pr:`401`) MPI - remove MPI setup for MDI. This eliminates a bug where interfering
  MPI environment variables were getting set upon ``import qcengine`` when pymdi and mpi4py packages
  were present. @awvwgk, @taylor-a-barnes
- (:pr:`418`, :pr:`389`, :issue:`292`) Psi4 - make Psi4 exe/module detection and version parsing more robust. @Flamefire, @coltonbh, @loriab

Misc.
+++++
- (:pr:`419`) CI - remove disabled LGTM and update badges. @loriab
- (:pr:`422`) CI - turn on crontab CI running to better notice external trouble. @loriab


v0.27.0 / 2023-08-02
--------------------

Bug Fixes
+++++++++
- (:pr:`414`) Import `pydantic.v1` from pydantic v2 so that QCEngine can work with any >=1.8.2 pydantic
  until QCEngine is updated for v2. If using v2, use QCElemental >=v0.26.0 that has a similar change.
  QCEngineRecords received similar treatment. @Lnaden, @loriab
- (:pr:`414`) Versioneer - update so works with Python 3.12.
- (:pr:`414`) Maintenance
   - Sphinx - fix build errors.
   - Lint - pin black to 2022 format.
   - GHA - switch to mamba solver. @loriab
- (:pr:`394`) Entos/Qcore - updated model environments. @loriab


v0.26.0 / 2022-11-30
--------------------

Breaking Changes
++++++++++++++++

- (:pr:`385`) Dispersion - the dispersion parameters resources file has been altered so that for D3 variants there's a
  2b set (e.g., d3bj2b) that is pure 2-body and doesn't accept s9 (effectively fixed at 0.0) and a atm set (e.g.,
  d3zeroatm) that does accept s9 (by default 1.0 but user-variable). Previous D3 levels are aliased to 2b. Only
  downstreams that call the dispersion resources directly should be affected, and retrofits are in place for the known
  victim/instigator (Psi4). @loriab

New Features
++++++++++++

Enhancements
++++++++++++
- (:pr:`380`) MRChem - added gradient and thus geometry optimizations support. @robertodr
- (:pr:`385`) dftd3 - the classic interface now accepts e.g., ``d3mbj2b`` as a level hint. @loriab
- (:pr:`385`) s-dftd3 - added keyword ``apply_qcengine_aliases`` that when True and ``level_hint`` present allows the
  levels and aliases in the dispersion resources (e.g., ``d3``, ``d3atm``, ``d32b``) to be given as ``level_hint``. The
  resource parameters are passed to s-dftd3 as param_tweaks. @loriab

Bug Fixes
+++++++++
- (:pr:`383`) yaml - uses safe loading. @mbanck, @loriab
- (:pr:`385`) dftd3 - the pairwise analysis requested through ``AtomicInput.keywords["pair_resolved"] = True`` and
  returned in ``AtomicResult.extras["qcvars"]["2-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS"]`` was elementwise too
  large by a factor of 2. It now matches the ``s-dftd3`` harness and fulfills that the sum of the array equals the
  2-body dispersion energy. @loriab


v0.25.0 / 2022-11-11
--------------------

Breaking Changes
++++++++++++++++
- (:pr:`376`) GAMESS - slight breaking changes of (1) ROHF MP2 ZAPT quantities now stored in "ZAPT" variables, not "MP2"
  variables; and (2) "HF TOTAL ENERGY" no longer stores DFT energy in DFT computation. @loriab
- (:pr:`376`) testing - reference quantities now indexed by "standard" or "semicanonical" orbitals since program defaults
  differ (mostly in CCSD ROHF FC). Downstream projects using the stdsuite interface will need to add an extra argument to query
  reference data. @loriab

New Features
++++++++++++

Enhancements
++++++++++++
- (:pr:`376`) Cfour - added parsing for BCCD and BCCD(T) methods. @loriab
- (:pr:`376`) NWChem - B2PLYP double-hybrid can now be run and parsed. Added CC2 parsing. @loriab
- (:pr:`376`) testing - added parsing contracts for ZAPT2, CEPA(1), CEPA(3), ACPF, AQCC, BCCD, BCCD(T), CC2, CC3, and DH-DFT. Added conventional references for most. @loriab
- (:pr:`378`) OpenFF - Support OpenFF Toolkit v0.11+. @Yoshanuikabundi

Bug Fixes
+++++++++


v0.24.1 / 2022-08-16
--------------------

Enhancements
++++++++++++
- (:pr:`375`) testing - in standard suite, add reference values for occ, a-ccsd(t), olccd grad, remp2, omp2, omp2.5, omp3, oremp2, density fitted ccsd, ccsd(t), a-ccsd(t). @loriab


v0.24.0 / 2022-07-08
--------------------

Upcoming Breaking Changes
+++++++++++++++++++++++++
- (:pr:`372`) QCSchema - changes are planned to schema layout and QCEngine API that will be outlined in an issue. These are not expected to involve detailed changes to the harnesses, and update helper functions will be provided. In preparation, QCEngine is pinned to a maximum compatible QCElemental v0.25.0 (current release). Projects using QCSchema through QCElemental are advised to pin to maximum v0.25.0 qcel and v0.24.0 to avert trouble, since this is our first experience with schema increments. @loriab

New Features
++++++++++++
- (:pr:`343`) DFT-D3 - added the ``SDFTD3Harness`` to handle DFT-D3 via a Python API. This has native QCSchema support and programmatic access to the parameter database. @awvwgk
- (:pr:`353`) TeraChem - added the ``TeraChemFrontEndHarness`` to handle file I/O in extension to the protocol buffer ``TeraChemPBSHarness`` interface. @coltonbh

Enhancements
++++++++++++
- (:pr:`350`) Rename the ``compute(..., local_options)`` argument to ``compute(..., task_config)``. Former still works and will for a while. @coltonbh
- (:pr:`361`) testing - in standard suite, add references for Hartree--Fock density-fitten Hessians. @loriab
- (:pr:`362`) docs - update setup with theme and fuller information on Pydantic models. @loriab
- (:pr:`363`) CFOUR - learned not to set ``DERIV_LEVEL`` when ``atomicinput.driver=properties``. Helps properties like DBOC. @loriab
- (:pr:`363`) Allow directory structure in ``execute(..., infiles)`` argument, not just flat-level files. @loriab
- (:pr:`364`) CFOUR - learned to harvest gradients when ghost atoms involved. Any CFOUR job with ghost atoms involves a hack that may go amiss when Xenon atoms in target molecule. @loriab
- (:pr:`364`) NWChem - learned to handle keyword ``geometry__autosym`` to tighten or loosen automatic symmetry detection. @loriab
- (:pr:`372`) testing - 2022 OpenMopac now actively tested in GHA. Note fields and output slightly different since 2019 harness. @awvwgk, @loriab

Bug Fixes
+++++++++
- (:pr:`301`, :pr:`367`) PyBerny - learned how to fail informatively when something goes wrong instead of assuming all is well and failing misleadingly while processing success. @coltonbh
- (:pr:`333`) NWChem - learned to skip writing the original ``atomicinput.molecule`` geometry to the input file only when both (1) the job is known to be part of a restart and (2) the job originates from the NWChem "driver" (that is, the optimizer). Previously, the geometry writing was skipped under (1) circumstances, so single-point e/g/h didn't have a geometry to work from. @WardLT
- (:pr:`349`) Turbomole - learned to correctly parse Hessian files when molecule contains more than 33 atoms. @eljost

Misc.
+++++
- (:pr:`354`, :pr:`356`, :pr:`361`, :pr:`366`, :pr:`368`) CI updates and fixes and changelog. @coltonbh, @loriab


v0.23.0 / 2022-03-10
--------------------

Enhancements
++++++++++++
- (:pr:`351`) Torsiondrive procedure refactored to make it easier for users to implement a parallel version via subclassing and overwriting the `_spawn_optimizations` method. @jthorton


v0.22.0 / 2022-01-25
--------------------

Bug Fixes
+++++++++
- (:pr:`338`) Correctly export version to tarballs created by git-archive. @mbanck, @loriab
- (:pr:`339`) QCEngine now tolerant of `cpuinfo` failure to populate `brand_raw`, `brand`. @dotsdl, @loriab, @WardLT


v0.21.0 / 2021-11-22
--------------------

Enhancements
++++++++++++
- (:pr:`321`) CFOUR, GAMESS, NWChem, Psi4, DFTD3, MP2D, gCP - learned to return certain native text
  files under control of the ``native_files`` protocol. GAMESS users are strongly advised to at
  least set ``protocols.native_files = "input"`` so that the job is reproducible. @loriab
- (:pr:`325`) Torsiondrive - learned to use multiple molecules as input to torsiondrives. @jthorton
- (:pr:`327`) TorchANI - learned to use GPUs if available. @kexul
- (:pr:`330`, :pr:`332`) NWChem - learned to restart from existing scratch if QCEngine is killed. @WardLT


v0.20.1 / 2021-10-08
--------------------

Bug Fixes
+++++++++

- (:pr:`322`) Psi4 - allowed more test cases with gradients and Hessians after a compatibility PR started
  saving them. @loriab
- (:pr:`323`) Turbomole - learned to store calcinfo_natom so that gradients and Hessians can be computed
  after QCElemental started using that quantity for shape checking in
  [v0.22.0](https://github.com/MolSSI/QCElemental/blob/master/docs/source/changelog.rst#0220--2021-08-26)
  @eljost


v0.20.0 / 2021-10-01
--------------------

New Features
++++++++++++
- (:pr:`305`) TorsionDrive - new procedure to automate constrained optimizations along a geometry
  grid. Akin to the longstanding QCFractal TorsionDrive service. @SimonBoothroyd

Enhancements
++++++++++++
- (:pr:`307`) NWChem - learns to automatically increase the number of iterations when SCF, CC, etc.
  fails to converge. @WardLT
- (:pr:`309`) ``qcengine info`` learned to print the location of found CMS programs, and geometric,
  OpenMM, and RDKit learned to return their versions. @loriab
- (:pr:`311`) CFOUR, GAMESS, NWChem harnesses learned to notice which internal module performs a calc
  (e.g., tce/cc for NWChem) and to store it in ``AtomicResult.provenance.module``. Psi4 already does
  this. @loriab
- (:pr:`312`) CFOUR, GAMESS, NWChem harnesses learned to run and harvest several new methods in the
  MP, CC, CI, DFT families. @loriab
- (:pr:`316`) Config - ``TaskConfig`` learned a new field ``scratch_messy`` to instruct a
  ``qcng.compute()`` run to not clean up the scratch directory at the end. @loriab
- (:pr:`316`) GAMESS - harness learned to obey ncores and scratch_messy local_config options. When
  ``ncores > 1``, the memory option is partitioned into replicated and non after exetyp=check trials. @loriab
- (:pr:`316`) Psi4 - harness learned to obey scratch_messy and memory local_config options. Memory
  was previously off by a little (GB vs GiB). @loriab
- (:pr:`316`) CFOUR - harness learned to obey scratch_messy and memory local_config options. Memory
  was previously off by a little. @loriab
- (:pr:`316`) NWChem - harness learned to obey scratch_messy and memory local_config options. Memory
  was previously very off for v7. @loriab
- (:pr:`315`) CFOUR, GAMESS, NWChem -- learned to return in AtomicInput or program native orientation
  depending on fix_com & fix_orientation= T or F. Psi4 already did this. Previously these three
  always returned AtomicInput orientation. Note that when returning program native orientation, the
  molecule is overwritten, so AtomicResult is not a superset of AtomicInput. @loriab
- (:pr:`315`) CFOUR, GAMESS, NWChem -- learned to harvest gradients and Hessians. @loriab
- (:pr:`317`) Docs - start "new harness" docs, so contributors have a coarse roadmap. @loriab
- (:pr:`318`) Docs - documentation is now served from https://molssi.github.io/QCEngine/ and built
  by https://github.com/MolSSI/QCEngine/blob/master/.github/workflows/CI.yml .
  https://qcengine.readthedocs.io/en/latest/ will soon be retired. @loriab
- (:pr:`320`) CFOUR, NWChem -- learned to run with ghost atoms, tentatively. @loriab

Bug Fixes
+++++++++
- (:pr:`313`, :pr:`319`) OpenMM - accommocate both old and new simtk/openmm import patterns. @dotsdl


v0.19.0 / 2021-05-16
--------------------

New Features
++++++++++++
- (:pr:`290`) MCTC-GCP - harness for new implementation of gCP, `mctc-gcp`, whose cmdline interface is drop-in replacement. @loriab
- (:pr:`291`) DFTD4 - new harness for standalone DFT-D4 executable. @awvwgk
- (:pr:`289`) TeraChem - new harness for TeraChem Protocol Buffer Server mode. @coltonbh

Enhancements
++++++++++++
- (:pr:`288`) GAMESS, Cfour, NWChem - add calcinfo harvesting, HF and MP2 gradient harvesting. @loriab

Bug Fixes
+++++++++
- (:pr:`288`) Avert running model.basis = BasisSet schema even though they validate. @loriab
- (:pr:`294`) NWChem - fixed bug where was retrieving only the first step in a geometry relaxation with line-search off. @WardLT
- (:pr:`297`) MDI - Update interface for v1.2. @loriab


v0.18.0 / 2021-02-16
--------------------

New Features
++++++++++++
- (:pr:`206`) OptKing - new procedure harness for OptKing optimizer. @AlexHeide
- (:pr:`269`) MRChem - new multiresolution chemistry program harness. @robertodr
- (:pr:`277`) ADCC - new program harness for ADC-connect. (Requires Psi4 for SCF.) @maxscheurer
- (:pr:`278`) gCP - new program harness for geometric counterpoise. @hokru
- (:pr:`280`) Add framework to register identifying known outfile errors, modify input schema, and rerun. @WardLT
- (:pr:`281`) NWChem - new procedure harness to use NWChem's DRIVER geometry optimizer with NWChem's program harness gradients. @WardLT
- (:pr:`282`) DFTD3 - added D3m and D3m(bj) parameters for SAPT0/HF. Allow pairwise analysis to be returned. @jeffschriber

Enhancements
++++++++++++
- (:pr:`274`) Entos/Qcore - renamed harness and updated to new Python bindings. @dgasmith
- (:pr:`283`) OpenMM - transition harness from `openforcefield` packages on omnia channel to `openff.toolkit` packages on conda-forge channel. @SimonBoothroyd
- (:pr:`286`, :pr:`287`) CI - moves from Travis-CI to GHA for open-source testing. @loriab

Bug Fixes
+++++++++
- (:pr:`273`) TeraChem - fixed bug of missing method field. @stvogt


v0.17.0 / 2020-10-02
--------------------

New Features
++++++++++++
- (:pr:`262`) Add project authors information. @loriab

Enhancements
++++++++++++
- (:pr:`264`) Turbomole - add analytic and finite difference Hessians. @eljost
- (:pr:`266`) Psi4- error messages from Psi4Harness no longer swallowed by `KeyError`. @dotsdl

Bug Fixes
+++++++++
- (:pr:`264`) Turbomole - fix output properties handling. @eljost
- (:pr:`265`) xtb - ensure extra tags are preserved in XTB harness. @WardLT
- (:pr:`270`) TorchANI - now lazily loads models as requested for compute. @dotsdl


v0.16.0 / 2020-08-19
--------------------

New Features
++++++++++++

Enhancements
++++++++++++
- (:pr:`241`) NWChem - improved performance by turning on ``atoms_map=True``, which does seem to be true. @WardLT
- (:pr:`257`) TorchANI - learned the ANI2x model and to work with v2. @farhadrgh
- (:pr:`259`) Added MP2.5 & MP3 energies and HF, MP2.5, MP3, LCCD gradients reference data to stdsuite. @loriab
- (:pr:`261`) Q-Chem - learned to return more informative Provenance, learned to work with v5.1. @loriab
- (:pr:`263`) NWChem - learned how to turn off automatic Z-Matrix coordinates with ``geometry__noautoz = True``. @WardLT

Bug Fixes
+++++++++
- (:pr:`261`) Molpro - learned to error cleanly if version too old for XML parsing. @loriab
- (:pr:`261`) Q-Chem - learned to extract version from output file instead of ``qchem -h`` since command isn't available
  from a source install. @loriab


v0.15.0 / 2020-06-26
--------------------

New Features
++++++++++++
- (:pr:`232`) PyBerny - new geometry optimizer procedure harness. @jhrmnn
- (:pr:`238`) Set up testing infrastructure, "stdsuite", where method reference values and expected results names (e.g.,
  total energy and correlation energy from MP2) are stored here in QCEngine but may be used from anywhere (presently,
  Psi4). Earlier MP2 and CCSD tests here converted to new scheme, removing ``test_standard_suite_mp2.py`` and ``ccsd``.
- (:pr:`249`, :pr:`254`) XTB - new harness for xtb-python that natively speaks QCSchema. @awvwgk

Enhancements
++++++++++++
- (:pr:`230`) NWChem - improved dipole, HOMO, LUMO harvesting.
- (:pr:`233`) ``qcng.util.execute`` learned argument ``exit_code`` above which to fail, rather than just ``!= 0``.
- (:pr:`234`) MDI - harness updated to support release verion v1.0.0 .
- (:pr:`238`) Cfour, GAMESS, NWChem -- harnesses updated to collect available spin components for MP2 and CCSD.
  Also updated to set appropriate ``qcel.models.AtomicProperties`` from collected QCVariables.
- (:pr:`239`) OpenMM - OpenMM harness now looks for cmiles information in the
  molecule extras field when typing. Also we allow for the use of gaff
  forcefields. @jthorton
- (:pr:`243`) NWChem - more useful stdout error return.
- (:pr:`244`) Added CCSD(T), LCCD, and LCCSD reference data to stdsuite. @loriab
- (:pr:`246`) TorchANI - harness does not support v2 releases.
- (:pr:`251`) DFTD3 - added D3(0) and D3(BJ) parameters for PBE0-DH functional.

Bug Fixes
+++++++++
- (:pr:`244`) Psi4 - fixed bug in ``extras["psiapi"] == True`` mode where if calc failed, error not handled by QCEngine. @loriab
- (:pr:`245`) Added missing import to sys for ``test_standard_suite.py``. @sjrl
- (:pr:`248`) NWChem - fix HFexch specification bug.
- Psi4 -- QCFractal INCOMPLETE state bug https://github.com/MolSSI/QCEngine/issues/250 fixed by https://github.com/psi4/psi4/pull/1933 .
- (:pr:`253`) Make compatible with both py-cpuinfo 5 & 6, fixing issue 252.


v0.14.0 / 2020-02-06
--------------------

New Features
++++++++++++
- (:pr:`212`) NWChem - Adds CI for the NWChem harness.
- (:pr:`226`) OpenMM - Moves the OpenMM harness to a canonical forcefield based method/basis language combination.
- (:pr:`228`) RDKit - Adds MMFF94 force field capabilities.

Enhancements
++++++++++++
- (:pr:`201`) Psi4 - ``psi4 --version`` collection to only grab the last line.
- (:pr:`202`) Entos - Adds wavefunction parsing.
- (:pr:`203`) NWChem - Parses DFT empirical dispersion energy.
- (:pr:`204`) NWChem - Allows custom DFT functionals to be run.
- (:pr:`205`) NWChem - Improved gradient output and added Hessian support for NWChem.
- (:pr:`215`) Psi4 - if Psi4 location can be found by either PATH or PYTHONPATH, harness sets up both subprocesses and API execution.
- (:pr:`215`) ``get_program`` shows the helpful "install this" messages from ``found()`` rather than just saying "cannot be found".

Bug Fixes
+++++++++
- (:pr:`199`) Fix typo breaking NWChem property parsing.
- (:pr:`215`) NWChem complains *before* a calculation if the necessary ``networkx`` package not available.
- (:pr:`207`) NWChem - Minor bug fixes for NWChem when more than core per MPI rank is used.
- (:pr:`209`) NWChem - Fixed missing extras tags in NWChem harness.


v0.13.0 / 2019-12-10
--------------------

New Features
++++++++++++
- (:pr:`151`) Adds a OpenMM Harness for evaluation of SMIRNOFF force fields.
- (:pr:`189`) General MPI support and MPI CLI generator.

Enhancements
++++++++++++
- (:pr:`175`) Allows specifications for ``nnodes`` to begin MPI support.
- (:pr:`177`) NWChem - Parsing updates including Hessian abilities.
- (:pr:`180`) GAMESS - Output properties improvements.
- (:pr:`181`) NWChem - Output properties improvements.
- (:pr:`183`) Entos - Hessian and XTB support.
- (:pr:`185`) Entos - Improved subcommand support.
- (:pr:`187`) QChem - Support for raw log files without the binary file requirements and improved output properties support.
- (:pr:`188`) Automatic buffer reads to prevent deadlocking of process for very large outputs.
- (:pr:`194`) DFTD3 - Improved error message on failed evaluations.
- (:pr:`195`) Blackens the code base add GHA-based lint checks.

Bug Fixes
+++++++++
- (:pr:`179`) QChem - fixes print issue when driver is of an incorrect value.
- (:pr:`190`) Psi4 - fixes issues for methods without basis sets such as HF-3c.

v0.12.0 / 2019-11-13
--------------------

New Features
++++++++++++

- (:pr:`159`) Adds MolSSI Driver Interface support.
- (:pr:`160`) Adds Turbomole support.
- (:pr:`164`) Adds Q-Chem support.

Enhancements
++++++++++++

- (:pr:`155`) Support for Psi4 Wavefunctions using v1.4a2 or greater.
- (:pr:`162`) Adds test for geometry optimization with trajectory protocol truncation.
- (:pr:`167`) CFOUR and NWChem parsing improvements for CCSD(T) properties.
- (:pr:`168`) Standardizes on ``dispatch.out`` for the common output files.
- (:pr:`170`) Increases coverage and begins a common documentation page.
- (:pr:`171`) Add Molpro to the standard suite.
- (:pr:`172`) Models renamed according to https://github.com/MolSSI/QCElemental/issues/155, particularly ``ResultInput`` -> ``AtomicInput``, ``Result`` -> ``AtomicResult``, ``Optimization`` -> ``OptimizationResult``.

Bug Fixes
+++++++++


v0.11.0 / 2019-10-01
--------------------

New Features
++++++++++++

- (:pr:`162`) Adds a test to take advantage of Elemental's `Protocols <https://github.com/MolSSI/QCElemental/pull/140>`_.
  Although this PR does not technically change anything in Engine, bumping the minor version here allows
  upstream programs to note when this feature was available because the minimum version dependency on Elemental
  has been bumped as well.


Enhancements
++++++++++++

- (:pr:`143`) Updates to Entos and Molpro to allow Entos to execute functions from the Molpro Harness. Also helps
  the two drivers to conform to :pr:`86`.
- (:pr:`145`, :pr:`148`) Initial CLI tests have been added to help further ensure Engine is running proper.
- (:pr:`149`) The GAMESS Harness has been improved by adding testing.
- (:pr:`150`, :pr:`153`) TorchANI has been improved by adding a Hessian driver to it and additional information
  is returned in the ``extra`` field when ``energy`` is the driver.
  This also bumped the minimum version of TorchANI Engine supports from 0.5 to 0.9.
- (:pr:`154`) Molpro's harness has been improved to support ``callinfo_X`` properties, unrestricted HF and DFT
  calculations, and the initial support for parsing local correlation calculations.
- (:pr:`158`) Entos' output parsing has been improved to read the json dictionary produced by the program
  directly. Also updates the input file generation.
- (:pr:`161`) Updates MOPAC to have more sensible quantum-chemistry like keywords by default.

Bug Fixes
+++++++++
- (:pr:`156`) Fixed a compatibility bug in specific version of Intel-OpenMP by skipping version
  2019.5-281.
- (:pr:`161`) Improved error handling in MOPAC if the execution was incorrect.


v0.10.0 / 2019-08-25
--------------------

New Features
++++++++++++

- (:pr:`132`) Expands CLI for ``info``, ``run``, and ``run-procedure`` options.
- (:pr:`137`) A new CI pipeline through Azure has been developed which uses custom, private Docker images
  to house non-public code which will enable us to test Engine through integrated CI on these codes securely.
- (:pr:`140`) GAMESS, CFOUR, NWChem preliminary implementations.

Enhancements
++++++++++++

- (:pr:`138`) Documentation on Azure triggers.
- (:pr:`139`) Overhauls install documentation and clearly defines dev install vs production installs.



v0.9.0 / 2019-08-14
-------------------

New Features
++++++++++++

- (:pr:`120`) Engine now takes advantage of Elemental's new Msgpack serialization option for Models. Serialization
  defaults to msgpack when available (``conda install msgpack-python [-c conda-forge]``), falling back to JSON
  otherwise. This results in substantial speedups for both serialization and deserialization actions and should be a
  transparent replacement for users within Engine and Elemental themselves.

Enhancements
++++++++++++

- (:pr:`112`) The ``MolproHarness`` has been updated to handle DFT and CCSD(T) energies and gradients.
- (:pr:`116`) An environment context manager has been added to catch NumPy style parallelization with Python functions.
- (:pr:`117`) MOPAC and DFTD3 can now accept an ``extras`` field which can pass around additional
  data, conforming to the rest of the Harnesses.
- (:pr:`119`) Small visual improvements to the docs have been made.
- (:pr:`120`) Lists inside models are now generally converted to numpy arrays for internal storage to maximize the
  benefit of the new Msgpack feature from Elemental.
- (:pr:`133`) The GAMESS Harness now collects the CCSD as part of its output.

Bug Fixes
+++++++++

- (:pr:`127`) Removed unused imports from the NWChem Harvester module.
- (:pr:`129`) Missing type hints from the ``MolproHarness`` have been added.
- (:pr:`131`) A code formatting redundancy in the GAMESS input file parser has been removed.

v0.8.2 / 2019-07-25
-------------------

Bug Fixes
+++++++++

- (:pr:`114`) Make compute and compute_procedure not have required kwargs while debugging
  a Fractal serialization issue. This is intended to be a temporary change and likely reverted
  in a later release

v0.8.1 / 2019-07-22
-------------------

Enhancements
++++++++++++

- (:pr:`110`) Psi4's auto-retry exception handlers now catch more classes of random errors

Bug Fixes
+++++++++

- (:pr:`109`) Geometric auto-retry settings now correctly propagate through the base code.

v0.8.0 / 2019-07-19
-------------------

New Features
++++++++++++

- (:pr:`95`, :pr:`96`, :pr:`97`, and :pr:`98`) The NWChem interface from QCDB has been added.
  Thanks to @vivacebelles and @jygrace for this addition!
- (:pr:`100`) The MOPAC interface has now been added to QCEngine thanks help to from @godotalgorithm.

Enhancements
++++++++++++

- (:pr:`94`) The gradient and molecule parsed from a GAMESS calculation output file are now returned in ``parse_output``
- (:pr:`101`) Enabled extra files in TeraChem scratch folder to be requested by users, collected after program
  execution, and recorded in the ``Result`` object as extras.
- (:pr:`103`) Random errors can now be retried a finite, controllable number of times (current default is zero retries).
  Geometry optimizations automatically set retries to 2. This only impacts errors which are categorized as
  ``RandomError`` by QCEngine and all other errors are raised as normal.

Bug Fixes
+++++++++

- (:pr:`99`) QCEngine now manages an explicit folder for each Psi4 job to write into and passes the scratch directory
  via ``-s`` command line. This resolves a key mismatch which could cause an error.
- (:pr:`102`) DFTD3 errors are now correctly returned as a ``FailedOperation`` instead of a raw ``dict``.


v0.7.1 / 2019-06-18
-------------------

Bug Fixes
+++++++++

- (:pr:`92`) Added an ``__init__.py`` file to the ``programs/tests`` directory so they are correctly bundled with the
  package.


v0.7.0 / 2019-06-17
-------------------

Breaking Changes
++++++++++++++++

- (:pr:`85`) The resource file ``programs.dftd3.dashparam.py`` has relocated and renamed to
  ``programs.empirical_dispersion_resources.py``.
- (:pr:`89`) Function ``util.execute`` forgot str argument ``scratch_location`` and learned ``scratch_directory`` in
  the same role of existing directory within which temporary directories are created and cleaned up. Non-user-facing
  function ``util.scratch_directory`` renamed to ``util.temporary_directory``.

New Features
++++++++++++

- (:pr:`60`) WIP: QCEngine interface to GAMESS can run the program (after light editing of rungms)
  and parse selected output (HF, CC, FCI) into QCSchema.
- (:pr:`73`) WIP: QCEngine interface to CFOUR can run the program and parse a variety of output into QCSchema.
- (:pr:`59`, :pr:`71`, :pr:`75`, :pr:`76`, :pr:`78`, :pr:`88`) Molpro improvements: Molpro can be run by QCEngine; and
  the input generator and output parser now supports CCSD energy and gradient calculations. Large thanks to
  @sjrl for many of the improvements
- (:pr:`69`) Custom Exceptions have been added to QCEngine's returns which will make parsing and
  diagnosing them easier and more programmatic for codes which invoke QCEngine. Thanks to @dgasmith for implementation.
- (:pr:`82`) QCEngine interface to entos can create input files (dft energy and gradients), run the program,
  and parse the output.
- (:pr:`85`) MP2D interface switched to upstream repo (https://github.com/Chandemonium/MP2D v1.1) and now produces
  correct analytic gradients.

Enhancements
++++++++++++

- (:pr:`62`, :pr:`67`, :pr:`83`) A large block of TeraChem improvements thanks to @ffangliu contributions.
  Changed the input parser to call qcelemental to_string method with bohr unit, improved output of parser to turn stdout
  into Result, and modified how version is parsed.
- (:pr:`63`) QCEngine functions ``util.which``, ``util.which_version``, ``util.parse_version``, and
  ``util.safe_version`` removed after migrating to QCElemental.
- (:pr:`65`) Torchani can now handle the ANI1-x and ANI1-ccx models. Credit to @dgasmith for implementation
- (:pr:`74`) Removes caching and reduces pytorch overhead from Travis CI. Credit to @dgasmith for implementation
- (:pr:`77`) Rename ``ProgramExecutor`` to ``ProgramHarness`` and ``BaseProcedure`` to ``ProcedureHarness``.
- (:pr:`77`) Function ``util.execute(..., outfiles=[])`` learned to collect output files matching a globbed filename.
- (:pr:`81`) Function ``util.execute`` learned list argument ``as_binary`` to handle input or output
  files as binary rather than string.
- (:pr:`81`) Function ``util.execute`` learned bool argument ``scratch_exist_ok`` to run in a preexisting directory.
  This is handy for stringing together execute calls.
- (:pr:`84`) Function ``util.execute`` learned str argument ``scratch_suffix`` to identify temp dictionaries for debugging.
- (:pr:`90`) DFTD3 now supports preliminary parameters for zero and Becke-Johnson damping to use with SAPT0-D

Bug Fixes
+++++++++

- (:pr:`80`) Fix "psi4:qcvars" handling for older Psi4 versions.


v0.6.4 / 2019-03-21
-------------------

Bug Fixes
+++++++++

- (:pr:`54`) Psi4's Engine implementation now checks its key words in a case insensitive way to give the same value
  whether you called Psi4 or Engine to do the compute.
- (:pr:`55`) Fixed an error handling routine in Engine to match Psi4.
- (:pr:`56`) Complex inputs are now handled better through Psi4's wrapper which caused Engine to hang while trying
  to write to ``stdout``.


v0.6.3 / 2019-03-15
-------------------

New Features
++++++++++++

- (:pr:`28`) TeraChem is now a registered executor in Engine! Thanks to @ffangliu for implementing.
- (:pr:`46`) MP2D is now a registered executor in Engine! Thanks to @loriab for implementing.

Enhancements
++++++++++++

- (:pr:`46`) ``dftd3``'s workings received an overhaul. The ``mol`` keyword has been replaced with ``dtype=2``,
  full Psi4 support is now provided, and an MP2D interface has been added.

Bug Fixes
+++++++++

- (:pr:`50` and :pr:`51`) Executing Psi4 on a single node with multiprocessing is more stable because Psi4 temps are
  moved to scratch directories. This behavior is now better documented with an example as well.
- (:pr:`52`) Psi4 calls are now executed through the ``subprocess`` module to prevent possible multiprocessing issues
  and memory leak after thousands of runs. A trade off is this adds about 0.5 seconds to task start-up, but its safe.
  A future Psi4 release will correct this issue and the change can be reverted.


v0.6.2 / 2019-03-07
-------------------

Enhancements
++++++++++++

- (:pr:`38` and :pr:`39`) Documentation now pulls from the custom QC Archive Sphinx Theme, but can fall back to the standard
  RTD theme. This allows all docs across QCA to appear consistent with each other.
- (:pr:`43`) Added a base model for all ``Procedure`` objects to derive from. This allows
  procedures' interactions with compute programs to be more unified. This PR also ensured
  GeomeTRIC provides Provenance information.

Bug Fixes
+++++++++
- (:pr:`40`) This PR improved numerous back-end and testing quality of life aspects.
  Fixed ``setup.py`` to call ``pytest`` instead of ``unittest`` when running tests on install.
  Some conda packages for Travis-CI are cached to reduce the download time of the larger computation codes.
  Psi4 is now pinned to the 1.3 version to fix build-level pin of libint.
  Conda-build recipe removed to avoid possible confusion for everyone who isn't a Conda-Forge
  recipe maintainer. Tests now rely exclusively on the ``conda env`` setups.


v0.6.1 / 2019-02-20
-------------------

Bug Fixes
+++++++++

- (:pr:`37`) Fixed an issue where RDKit methods were not case agnostic.

v0.6.0 / 2019-02-28
-------------------

Breaking Changes
++++++++++++++++

- (:pr:`36`) **breaking change** Model objects are returned by default rather than a dictionary.

New Features
++++++++++++

- (:pr:`18`) Add the ``dftd3`` program to available computers.
- (:pr:`29`) Adds preliminary support for the ``Molpro`` compute engine.
- (:pr:`31`) Moves all computation to ``ProgramExecutor`` to allow for a more flexible input generation, execution, output parsing interface.
- (:pr:`32`) Adds a general ``execute`` process which safely runs subprocess jobs.

Enhancements
++++++++++++

- (:pr:`33`) Moves the ``dftd3`` executor to the new ``ProgramExecutor`` interface.
- (:pr:`34`) Updates models to the more strict QCElemental v0.3.0 model classes.
- (:pr:`35`) Updates CI to avoid pulling CUDA libraries for ``torchani``.
- (:pr:`36`) First pass at documentation.


v0.5.2 / 2019-02-13
-------------------

Enhancements
++++++++++++

- (:pr:`24`) Improves load times dramatically by delaying imports and cpuutils.
- (:pr:`25`) Code base linting.
- (:pr:`30`) Ensures Psi4 output is already returned and Pydantic v0.20+ changes.

v0.5.1 / 2019-01-29
-------------------

Enhancements
++++++++++++

- (:pr:`22`) Compute results are now returned as a dict of Python Primals which have
  been serialized-deserialized through Pydantic instead of returning un-processed Python objects
  or json-compatible string.

v0.5.0 / 2019-01-28
-------------------

New Features
++++++++++++

- (:pr:`8`) Adds the TorchANI program for ANI-1 like energies and potentials.
- (:pr:`16`) Adds QCElemental models based off QCSchema to QCEngine for both validation and object-based manipulation of input and output data.

Enhancements
++++++++++++

- (:pr:`14`) Migrates option to Pydantic objects for validation and creation.
- (:pr:`14`) Introduces NodeDescriptor (for individual node description) and JobConfig (individual job configuration) objects.
- (:pr:`17`) NodeDescriptor overhauled to work better with Parsl/Balsam/Dask/etc.
