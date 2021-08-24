Changelog
=========

.. vX.Y.0 / 2021-MM-DD
.. --------------------
..
.. New Features
.. ++++++++++++
..
.. Enhancements
.. ++++++++++++
..
.. Bug Fixes
.. +++++++++


v0.20.0 / 2021-MM-DD
--------------------

New Features
++++++++++++

Enhancements
++++++++++++
- (:pr:`307`) NWChem - learns to automatically increase the number of iterations when SCF, CC, etc. fails to converge.
- (:pr:`309`) ``qcengine info`` learned to print the location of found CMS programs, and geometric, OpenMM, and RDKit learned to return their versions.

Bug Fixes
+++++++++


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
