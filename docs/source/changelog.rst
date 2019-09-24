Changelog
=========

.. vX.Y.0 / 2019-MM-DD
.. -------------------
..
.. New Features
.. ++++++++++++
..
.. Enhancements
.. ++++++++++++
..
.. Bug Fixes
.. +++++++++

v0.10.0 / 2019-09-25
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
