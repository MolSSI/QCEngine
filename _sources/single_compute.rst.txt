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

.. autoclass:: qcelemental.models.AtomicInput

Returned Fields
---------------

.. autoclass:: qcelemental.models.AtomicResult

