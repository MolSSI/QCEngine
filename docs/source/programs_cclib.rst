cclib-backed Q-Chem and ORCA
============================

QCEngine provides two optional, non-production program selectors backed by the
same cclib parsing harness:

* ``cclib-qchem`` executes Q-Chem 5.1 or newer.
* ``cclib-orca`` executes ORCA 6.0 or newer.

The existing native ``qchem`` selector remains available and is not replaced by
``cclib-qchem``.  The cclib selectors generate native input from an
``AtomicInput``, execute the selected program in QCEngine-managed scratch space,
and parse the complete primary output exclusively through cclib.

Supported inputs
----------------

The first implementation accepts only the following atomic calculations:

+--------------+-----------------+---------------------+
| Driver       | Q-Chem job type | ORCA support        |
+==============+=================+=====================+
| ``energy``   | ``sp``          | supported           |
+--------------+-----------------+---------------------+
| ``gradient`` | ``force``       | not supported       |
+--------------+-----------------+---------------------+
| ``hessian``  | ``freq``        | not supported       |
+--------------+-----------------+---------------------+

Method strings pass through to the native program rather than being restricted
by a harness-maintained allowlist.  Methods unsupported by the selected native
program fail downstream.  ORCA gradient and Hessian calculations are disabled
because current cclib parsing does not support the verified output path.

The basis must be a non-empty string.  Structured QCSchema basis objects are
not supported.  Molecules must contain only real atoms; ghost atoms are
rejected.  Charge, multiplicity, atom order, and Cartesian geometry come from
the QCSchema molecule.  Geometry is written in bohr for Q-Chem and converted
to Angstrom for ORCA.

These selectors do not accept raw native input, multi-job or restart files,
node launchers, or unsupported drivers.  ORCA's node-parallel harness flag is
scheduler capability metadata; the generated input uses ``TaskConfig.ncores``
in ``%pal``.

Program keywords and resources
------------------------------

Q-Chem keywords are a flat QCSchema mapping.  Keys are converted to uppercase
``$rem`` keys and sorted; values may be strings, booleans, integers, or finite
floats.  For example:

.. code-block:: python

   "keywords": {
       "scf_algorithm": "diis",
       "thresh": 10,
   }

The harness derives ``JOBTYPE``, ``METHOD``, ``BASIS``, ``MEM_TOTAL``, and
``INPUT_BOHR`` and adds harvesting settings.  Consequently, the complete
case-insensitive reserved-key set is:

.. code-block:: text

   JOBTYPE METHOD BASIS MEM_TOTAL INPUT_BOHR
   SCF_FINAL_PRINT PRINT_GENERAL_BASIS PRINT_ORBITALS MOLDEN_FORMAT

Attempting to override any reserved key is an input error.  Q-Chem uses total
memory from ``TaskConfig`` and runs the resolved executable with ``-nt`` and the
configured core count.

ORCA keywords have exactly two optional top-level entries, ``simple`` and
``blocks``:

.. code-block:: python

   "keywords": {
       "simple": ["rks", "usesym"],
       "blocks": {
           "output": "PrintLevel Normal\nPrint[P_Basis] 2",
           "scf": "MaxIter 200",
       },
   }

``simple`` must be a list of non-empty, single-line strings.  ``blocks`` must
map valid ORCA block names to string bodies.  Bodies are limited to block-local
lines: ``end`` and ``$new_job`` cannot be first tokens, and ``%`` or ``*`` cannot
be the first non-whitespace character.  User output lines follow the harness
defaults, and other blocks are emitted deterministically.  The case-insensitive
block names ``pal`` and ``maxcore`` are reserved because core
and per-process memory values come from ``TaskConfig``.  The ``coords`` block
and all other coordinate injection are also reserved; geometry comes only from
the QCSchema molecule.  Unknown top-level keyword entries are rejected.

Installation and discovery
--------------------------

cclib remains optional: importing QCEngine and using other programs does not
import it.  A development cclib containing the required QCSchema-writer
behavior can be installed without adding cclib to QCEngine's core dependencies,
for example:

.. code-block:: console

   python -m pip install "cclib @ git+https://github.com/cclib/cclib.git"

At availability checking time, QCEngine lazily imports the required cclib
interfaces.  Each executable is discovered only through ``PATH`` and then
checked for the expected program identity and minimum version.  Cached
executable versions are keyed by resolved path, so an unrelated executable
named ``orca`` does not satisfy ``cclib-orca``.  Real writer output is validated
as QCSchema v1 during result conversion.

Configure ``PATH`` and all proprietary-program environment variables *before*
starting or importing in the Python process that will call QCEngine.  The
harness never sources setup scripts and does not accept a selector-specific
absolute executable path.  For Q-Chem, the inherited ``QC``, ``QCAUX``, and
``QCPROG`` resources must be readable (and the program driver executable).
QCEngine supplies a managed ``QCSCRATCH`` to the child.  Child process
environments otherwise begin with a copy of the Python process environment.

Results and native files
------------------------

cclib auto-detects the parser, and the harness rejects a parser that does not
match the selected program.  It also rejects parsed driver, method, or basis
mismatches rather than relabeling the result.  ``QCSchemaWriter`` first emits a
QCSchema v1 dictionary, which is validated as a v1 ``AtomicResult`` and then
converted to QCEngine's internal QCSchema v2 representation.  Public callers
may still request a v1 result with ``return_version=1``.

The result keeps cclib's flat extras, such as available atom charges and
coordinates, orbital data, SCF histories, correlated energies, and cclib unit
annotations.  It adds one namespaced object:

.. code-block:: python

   result.extras["cclib_harness"] == {
       "selector": "cclib-qchem",  # or cclib-orca
       "cclib_version": "...",
       "parser": "QChem",          # or ORCA
       "executable": "...",        # resolved executable
   }

The parsed program remains the provenance creator, and the complete primary
program output is returned as top-level ``stdout``.  The
``protocols.native_files`` choices produce these native-file keys:

+-----------+-----------------------------------+
| Protocol  | ``native_files`` keys             |
+===========+===================================+
| ``none``  | none                              |
+-----------+-----------------------------------+
| ``input`` | ``input``                         |
+-----------+-----------------------------------+
| ``all``   | ``input`` and ``dispatch.out``    |
+-----------+-----------------------------------+

Native-file keys are never named ``stdout`` or ``stderr``.  If an older active
QCSchema v1 model rejects ``native_files``, the generated input is retained as
``extras["cclib_harness"]["native_input"]`` when requested and native files
are otherwise omitted.

Adding another cclib-backed program
-----------------------------------

The package layout is:

.. code-block:: text

   qcengine/programs/cclib_programs/
   |-- __init__.py
   |-- base.py
   |-- cclib_orca.py
   |-- cclib_qchem.py
   `-- tests/
       |-- __init__.py
       `-- test_cclib.py

``base.py`` owns shared lazy loading, execution, parsing, and QCSchema
conversion.  A ``cclib_<program>.py`` module owns native input, preflight,
version probe, output selection, parser type, and its concrete subclass.
``ProgramDefinition`` is the immutable callback contract consumed by
``CCLibHarness``.

Program modules must not import cclib at module import time, which preserves
cclib as an optional dependency.  Every helper must have a
responsibility-focused docstring that describes the single boundary it owns.

Fixture and live verification
-----------------------------

Fixture integration requires a cclib source checkout with its ``data``
directory.  Set the checkout only for the test process:

.. code-block:: console

   CCLIB_SOURCE_ROOT=/path/to/cclib \
     python -m pytest qcengine/programs/cclib_programs/tests/test_cclib.py -q -k fixture

The approved fixture set has exactly seven passing cases.  It contains only
supported Q-Chem calculations and ORCA energy calculations.

Live tests run only when the corresponding optional selector is available:

.. code-block:: console

   python -m pytest qcengine/programs/cclib_programs/tests/test_cclib.py -q -m cclib_qchem
   python -m pytest qcengine/programs/cclib_programs/tests/test_cclib.py -q -m cclib_orca
