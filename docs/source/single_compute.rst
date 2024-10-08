Single Compute
==============

QCEngine's primary purpose is to consume the MolSSI `QCSchema <https://github.com/MolSSI/QC_JSON_Schema>`_ and produce
QCSchema results for a variety of quantum chemistry, semiempirical, and molecular mechanics programs. Single QCSchema representation
comprises of a single ``energy``, ``gradient``, ``hessian``, or ``properties`` evaluation.

Input Description
-----------------

An input description has the following fields:

- ``molecule`` - A QCSchema compliant dictionary or Molecule model.
- ``driver`` - The ``energy``, ``gradient``, ``hessian``, or ``properties`` option.
- ``model`` - A description of the evaluation model. For quantum chemistry this is typically ``method`` and ``basis``. However,
  non-quantum chemistry models are often a simple ``method`` as in ``method = 'UFF'`` for forcefield evaluation.
- ``keywords`` - a dictionary of keywords to pass to the underlying program. These are program-specific keywords.

An example input is as follows:

.. code:: python

    >>> import qcengine as qcng
    >>> import qcelemental as qcel

    >>> mol = qcel.models.Molecule.from_data("""
    >>>     O  0.0  0.000  -0.129
    >>>     H  0.0 -1.494  1.027
    >>>     H  0.0  1.494  1.027
    >>> """)

    >>> inp = qcel.models.AtomicInput(
    >>>     molecule=mol,
    >>>     driver="energy",
    >>>     model={"method": "SCF", "basis": "sto-3g"},
    >>>     keywords={"scf_type": "df"}
    >>> )


Computation
-----------

A single computation can be evaluated with the ``compute`` function as follows:

.. code:: python

    >>> ret = qcng.compute(inp, "psi4")

By default the job is given resources relating to the compute environment it is in; however, these variables can be overridden:

.. code:: python

    >>> ret = qcng.compute(inp, "psi4", local_options={"memory": 2, "ncores": 3})

By default, the job returns a QCSchema Result of the same ``schema_version`` as the Input (v1 if Input version can't be determined). To request a specific version back, use the ``return_version`` keyword:

.. code:: python

    >>> ret = qcng.compute(inp, "psi4", return_version=2)


Results
-------

The results contain a complete record of the computation:

.. code:: python

    >>> ret.return_result
    -74.45994963230625

    >>> ret.properties.scf_dipole_moment
    [0.0, 0.0, 0.6635967188869244]

    >>> ret.provenance.cpu
    Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz


Input Fields
-------------

Note this is QCSchema v2 ``AtomicInput``. See the QCElemental docs for the
longstanding v1 model.

.. autopydantic_model:: qcelemental.models.v2.AtomicInput
   :noindex:

Returned Fields
---------------

Note this is QCSchema v2 ``AtomicResult``. See the QCElemental docs for the
longstanding v1 model.

.. autopydantic_model:: qcelemental.models.v2.AtomicResult
   :noindex:

FAQ
---

#. Where is scratch so I can access the CMS code's files?

   The QCArchive philosophy is that you shouldn't go looking in scratch for CMS-code-written files since the scratch directory is deleted automatically by QCEngine and even if preserved may be subject to autodeletion if run from a cluster. Instead, QCEngine brings back the primary input and output and any ancillary files from which it can harvest results. Whether these are returned to the user in ``AtomicResult`` can be controlled through protocols in the input like ``atomicinput.protocols.stdout = True`` and eventually (https://github.com/MolSSI/QCElemental/pull/275) ``atomicinput.protocols.native_files = "all"``.

   Nevertheless, you can, of course, access the scratch directory and CMS-code-written files. Pass an existing directory to the compute command (this directory will be parent) and tell it to not delete after the run: ``qcng.compute(..., local_options={"scratch_directory": "/existing/parent/dir", "scratch_messy": True})``.

#. sdfs
