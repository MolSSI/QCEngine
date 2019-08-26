.. QCEngine documentation master file, created by
   sphinx-quickstart on Fri Aug 17 09:45:43 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========
QCEngine
=========

*Quantum chemistry program executor and IO standardizer (QCSchema) for quantum chemistry.*

Program Execution
-----------------

A simple example of QCEngine's capabilities is as follows:

.. code:: python

    >>> import qcengine as qcng
    >>> import qcelemental as qcel

    >>> mol = qcel.models.Molecule.from_data("""
    O  0.0  0.000  -0.129
    H  0.0 -1.494  1.027
    H  0.0  1.494  1.027
    """)

    >>> input = qcel.models.ResultInput(
        molecule=mol,
        driver="energy",
        model={"method": "SCF", "basis": "sto-3g"},
        keywords={"scf_type": "df"}
        )

These input specifications can be executed with the ``compute`` syntax along with a program specifier:

.. code:: python

    >>> ret = qcng.compute(input, "psi4")

The results contain a complete record of the computation:

.. code:: python

    >>> ret.return_result
    -74.45994963230625

    >>> ret.properties.scf_dipole_moment
    [0.0, 0.0, 0.6635967188869244]

    >>> ret.provenance.cpu
    Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz

Backends
--------

Currently available compute backends for single results are as follow:

- Quantum Chemistry:

  - `Terachem <http://www.petachem.com>`_
  - `Molpro <https://www.molpro.net>`_
  - `Psi4 <http://www.psicode.org>`_

- AI Evaluation:

  - `TorchANI <https://github.com/aiqm/torchani>`_

- Molecular Mechanics:

  - `RDKit <http://rdkit.org>`_

In addition, several procedures are available:

- Geometry Optimization:

  - `geomeTRIC <https://github.com/leeping/geomeTRIC>`_

Information about detected compute backends and procedures may be printed using the :doc:`CLI <cli>`::

    qcengine info --programs
    qcengine info --procedures

Configuration Determination
---------------------------

In addition, QCEngine can automatically determine the following quantites:

- The number of physical cores on the system and to use.
- The amount of physical memory on the system and the amount to use.
- The provenance of a computation (hardware, software versions, and compute resources).
- Location of scratch disk space.
- Location of quantum chemistry programs binaries or Python modules.

Each of these options can be specified by the user as well.

.. code:: python

    >>> qcng.get_config()
    <JobConfig ncores=2 memory=2.506 scratch_directory=None>

    >>> qcng.get_config(local_options={"scratch_directory": "/tmp"})
    <JobConfig ncores=2 memory=2.506 scratch_directory='/tmp'>

    >>> os.environ["SCRATCH"] = "/my_scratch"
    >>> qcng.get_config(local_options={"scratch_directory": "$SCRATCH"})
    <JobConfig ncores=2 memory=2.506 scratch_directory='/my_scratch'>

Configuration information may be printed using the :doc:`CLI <cli>`::

   qcengine info --config

========

.. toctree::
    :maxdepth: 2
    :caption: Contents:



========

Index
-----

**Getting Started**

* :doc:`install`

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Getting Started

    install


**User Interface**

* :doc:`single_compute`
* :doc:`environment`
* :doc:`cli`

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: User Interface

    single_compute
    environment
    cli


**Developer Documentation**

* :doc:`api`
* :doc:`changelog`

.. toctree::
    :maxdepth: 1
    :caption: Developer Documentation

    api
    changelog
