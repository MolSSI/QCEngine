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
- ``model`` - A description of the evaluation model, for quantum chemistry this is typically ``method`` and ``basis``. However,
  non-quantum chemistry models are often a simple ``method`` as in ``method = 'UFF'`` for forcefield evaluation.
- ``keywords`` - a dictionary of keywords to pass to the underlying program, these are program-specific keywords.

An example input is as follows:

.. code:: python

    >>> import qcengine as qcng
    >>> import qcelemental as qcel

    >>> mol = qcel.models.Molecule.from_data("""
    O  0.0  0.000  -0.129
    H  0.0 -1.494  1.027
    H  0.0  1.494  1.027
    """)

    >>> inp = qcel.models.ResultInput(
        molecule=mol,
        driver="energy",
        model={"method": "SCF", "basis": "sto-3g"},
        keywords={"scf_type": "df"}
        )


Computation
-----------

A single computation can be evaluated with the ``compute`` function as follows:

.. code:: python

    >>> ret = qcng.compute(inp, "psi4")

By default the job is given resources relating to the compute environment it in; however, these variables can be overridden:

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

A short description of the fields is as follow:

- ``return_result`` - the direct return of the ``driver`` input. That is energy and gradient for a driver ``energy`` and ``gradient`` call, respectively.
- ``properties`` - Values associated with the ``return_result`` such as the ``scf_one_electron_energy``.
- ``stdout`` - The ``stdout`` or log of a programs run.
- ``provenance`` - A description of the calling program, version, wall time, etc.

A complete description of the input is also available in the output:

.. code:: python

    >>> ret.driver
    energy


Fields
------

A list of all fields is available through the ``fields`` property on the input and output:

.. code:: python

    >>> ret.driver
    ['molecule', 'driver', 'model', 'id', 'schema_name', 'schema_version', 'keywords',
     'extras', 'provenance', 'return_result', 'success', 'properties', 'stdout', 'stderr', 'error']

