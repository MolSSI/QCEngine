
Adding a New Program Harness
============================

Program harnesses are for community CMS codes that can independently
(or with a SCF bootstrap) compute single-point energies, derivatives,
or properties or components thereof (e.g., dispersion corrections).

A single CMS code generally has one program harness. However, if there
are drastically different ways of running a code (e.g., TeraChem text
input file and TeraChem PBS), separate harnesses may be created. Also,
if there are specialty capabilites not fitting into "single-point
energies ..." (e.g., GAMESS makefp task), an additional *procedure*
harness may be created.

This guide is a coarse path through adding a new program harness.

#. Open up communication with the QCEngine maintainers. Post an issue
   to GitHub and join the Slack channel (link off GH README) so you can
   get advice.

#. Copy a similar harness. Choose a program that your code behaves
   roughly like (mostly consider parsed vs. API access) and copy that
   harness, renaming it as your own and commenting out all but the
   structure. Search for the (copied) harness name to register it.

#. Fill in the ``_defaults`` section with program name and characteristics.

#. Fill in the ``def found`` function using ``which`` and
   ``which_import`` from QCElemental. See NWChem and OpenMM for examples
   of handling additional dependencies.

#. If your code's version can be extracted short of parsing an output
   file, fill in ``def get_version`` next. After this, ``> qcengine info``
   should show your code (provided it's in path).

#. To get a string output of QCSchema Molecule in your code's format,
   you may need to add a ``dtype`` at ``qcelemental/molparse/to_string.py``.

#. If your code's of the common translate-QCSchema-to-input, run,
   translate-output-to-QCSchema variety, next work on the ``def execute``
   function. This is fairly simple because it calls the powerful
   qcengine.util.execute to handle scratch, timeout, file writing and
   collection, etc. The harness function needs the names of input files
   (hard-code a string for now), the execution command, and the names
   of any scratch files to return for processing. Once ready, fill in
   ``def get_version`` if not done above.

#. Now fill in the short ``def compute`` entirely and ``def build_input``
   and ``def parse_output`` skeletally. Set up a simple molecule-and-model
   ``AtomicInput`` dictionary and run it with ``qcng.compute(atomicinput,
   "yourcode")`` to get something to iterate on.

#. Fill in ``def build_input`` to form your code's usual input format
   from the fields of ``AtomicInput``.

#. Fill in ``def parse_output`` to take results and put them
   into ``AtomicResult``. Most important is the ``return_result``
   field.  ``AtomicResultProperties`` can be populated when
   convenient. ``WavefunctionProperties`` is great but save for a later
   pass.

#. At this point your harness can correctly run one or more QCSchema
   inputs of your devising. Time to put it through paces. Register your
   code in ``_programs`` in ``testing.py``.  Most tests will then need a
   ``@using("yourcode")`` decorator so that they don't run (and fail the
   test suite) when your code isn't available.

#. Add basic tests to qcengine/tests/test_harness_canonical.py 
   * def test_compute_energy(program, model, keywords):
   * def test_compute_gradient(program, model, keywords):
   * def test_compute_energy_qcsk_basis(program, model, keywords):

#. Add basic failure tests to qcengine/tests/test_harness_canonical.py 
   * def test_compute_bad_models(program, model):

#. Add tests for the runtime config and to qcengine/programs/tests/test_canonical_config.py

#. For QM codes, consider adding lines to the
   qcengine/programs/tests/test_standard_suite.py to check energies,
   gradients, and Hessians against other codes.

#. For codes that can produce a Hartree--Fock, add lines to the
   qcengine/program/tests/test_alignment.py to check molecular and
   properties orientation handling.

#. If your code is available as a visible binary (e.g., pip, conda,
   docker download with no or trivial build), create a testing lane by
   adding to ``devtools/conda-envs/`` and ``.github/workflows/CI.yml``.
   This will check your code for every PR. We're looking into private
   testing for codes that aren't available.

#. Throughout, talk with the maintainers with questions. Error handling,
   especially, is intricate.

