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

+------------+-----------------+---------------------+
| Driver     | Q-Chem job type | ORCA driver keyword |
+============+=================+=====================+
| ``energy`` | ``sp``          | none                |
+------------+-----------------+---------------------+
| ``gradient`` | ``force``     | ``engrad``          |
+------------+-----------------+---------------------+
| ``hessian`` | ``freq``       | ``freq``            |
+------------+-----------------+---------------------+

Methods are case-insensitive ``hf``, ``b3lyp``, ``bp86``, ``mp2``, and
``ccsd``.  The basis must be a non-empty string.  Structured QCSchema basis
objects are not supported.  Molecules must contain only real atoms; ghost atoms
are rejected.  Charge, multiplicity, atom order, and Cartesian geometry come
from the QCSchema molecule.  Geometry is written in bohr for Q-Chem and
converted to Angstrom for ORCA.

These selectors do not accept raw native input, multi-job or restart files,
node launchers, or drivers and methods outside the lists above.  ORCA's
node-parallel harness flag is scheduler capability metadata; the generated
input uses ``TaskConfig.ncores`` in ``%pal``.

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
map valid ORCA block names to string bodies.  User output lines follow the
harness defaults, and other blocks are emitted deterministically.  The
case-insensitive block names ``pal`` and ``maxcore`` are reserved because core
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

At availability checking time, QCEngine runs a small compatibility probe that
validates the cclib writer's QCSchema v1 output, geometry units, and flat
extras.  Each executable is discovered only through ``PATH`` and then checked
for the expected program identity and minimum version.  Cached executable
versions are keyed by resolved path, so an unrelated executable named ``orca``
does not satisfy ``cclib-orca``.

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

Fixture and live verification
-----------------------------

Fixture integration requires a cclib source checkout with its ``data``
directory.  A portable invocation is:

.. code-block:: console

   CCLIB_SOURCE_ROOT=/path/to/cclib \
     python -m pytest qcengine/programs/tests/test_cclib.py -q -k fixture

The approved fixture set has seven passing cases and one intentional skip:
ORCA ``dvb_ir`` Hessian is deferred because the successful cclib parse lacks
``metadata.functional`` and its QCSchema writer raises ``KeyError``.

Live tests run only when the corresponding optional selector is available:

.. code-block:: console

   python -m pytest qcengine/programs/tests/test_cclib.py -q -m "cclib-qchem"
   python -m pytest qcengine/programs/tests/test_cclib.py -q -m "cclib-orca"
   python qchem_water_mp2.py
   python orca_water_ccsd.py

Each demonstration writes a complete ``*.result.json`` file and reports one
``PASS`` or ``FAIL`` line per comparison.  The Q-Chem MP2 demonstration checks
all six SCF history rows, not only the iteration count:

.. code-block:: python

   [[[0.398], [0.0668], [0.00822], [0.0016], [2.83e-5], [8.23e-6]]]

The nested shape and every row value are compared with an absolute ``1e-6``
tolerance.  The demonstration also checks the six required numerical
properties, calculation dimensions, seven historical MO energies, Mulliken
charges, schema/molecule/provenance, and representative flat extras.  The ORCA
CCSD demonstration selects strict extra references from the exact
``major.minor.patch`` provenance version: ``6.0.x`` and ``6.1.x`` each have
version-aware atom-charge, coordinate, orbital, CCSD, and SCF references.
Floating-point extra values use an absolute ``1e-6`` tolerance and discrete
values are exact.  Missing, malformed, or unsupported minor versions (including
``6.2.x``) fail
the demonstration's reference check rather than silently using another
version's values.  Known anomalous MP2 fields in the historical ORCA CCSD
writer output must be present but are deliberately not numerically endorsed.

Non-portable verification example
---------------------------------

.. warning::

   The commands and paths in this section are machine-local examples only.
   They are not portable configuration, are not used by the harness or unit
   tests, and must be adapted for another installation.  Initialize this shell
   before starting Python.

.. code-block:: bash

   export PATH=/projects/cos-lab-cs207/common/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg:$PATH
   source ~/qchem_vars.sh
   which qchem
   # /projects/cos-lab-cs207/common/software/qchem5.1/bin/qchem
   which orca
   # /projects/cos-lab-cs207/common/software/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca

   python -m pip install -e /home/awallace43/gits/cclib
   python -m pip install -e '/home/awallace43/gits/qcengine[test]'

   CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib \
     python -m pytest qcengine/programs/tests/test_cclib.py -q -k fixture
   python -m pytest qcengine/programs/tests/test_cclib.py -q -m "cclib-qchem"
   python -m pytest qcengine/programs/tests/test_cclib.py -q -m "cclib-orca"
   python qchem_water_mp2.py
   python orca_water_ccsd.py
