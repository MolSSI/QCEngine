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

    >>> import qcengine

    >>> job = {
      "schema_name": "qc_schema_input",
      "schema_version": 1,
      "molecule": {
        "geometry": [
          0.0,  0.0,  -0.129,
          0.0, -1.494, 1.027,
          0.0,  1.494, 1.027
        ],
        "symbols": [ "O", "H", "H"],
      },
      "driver": "energy",
      "model": {
        "method": "SCF",
        "basis": "sto-3g"
      },
      "keywords": {
        "scf_type": "df"
      }
    }

    >>> qcengine.compute(job, "psi4")
    {
      ...,
      "provenance": {
        "creator": "Psi4",
        "version": "1.2",
        "routine": "psi4.json.run_json",
        "memory": 2.235,
        "nthreads": 2,
        "cpu": "Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz",
        "hostname": "xxx.xxx.com",
        "username": "username",
        "wall_time": 0.8911292552947998
      },
      "return_result": -74.96475169985112,
      "properties": {
        "calcinfo_nbasis": 7,
        "calcinfo_nmo": 7,
        "calcinfo_nalpha": 5,
        "calcinfo_nbeta": 5,
        "calcinfo_natom": 3,
        "scf_one_electron_energy": -121.67648822935482,
        "scf_two_electron_energy": 37.91027446887428,
        "nuclear_repulsion_energy": 8.80146206062943,
        "scf_dipole_moment": [0.0, 0.0, 1.668684476563345],
        "scf_iterations": 6,
        "scf_total_energy": -74.96475169985112,
        "return_energy": -74.96475169985112
      }
    }


The QCEngine middleware can automatically determine:

- The number of physical cores on the system and to use.
- The amount of physical memory on the system and the amount to use.
- The provenance of a computation (hardware, software versions, and compute resources).
- Locations of scratch disk space.
- Locations of quantum chemistry programs.


QCArchive
---------

This module is part of the QCArchive project which sets out to answer the
the fundamental question of "How do we compile, aggregate, query, and share quantum
chemistry data to accelerate the understanding of new method performance,
fitting of novel force fields, and supporting the incredible data needs of
machine learning for computational molecular science?"

The QCArchive project is made up of three primary modules:

- `QCSchema <https://github.com/MolSSI/QC_JSON_Schema>`_ - A key/value schema for quantum chemistry.
- `QCEngine <https://github.com/MolSSI/QCEngine>`_ - A computational middleware to provide IO to a variety of quantum chemistry programs.
- `QCFractal <https://github.com/MolSSI/QCFractal>`_ - A distributed compute and database platform powered by QCEngine and QCSchema.

The QCArchive project's primary support comes from `The Molecular Sciences Software Institute <https://molssi.org>`_.

========

.. toctree::
   :maxdepth: 2
   :caption: Contents:



========

Index
-----

**Getting Started**

* :doc:`install`
* :doc:`community`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting Started

   install
   community


**Developer Documentationd**

* :doc:`api`
* :doc:`changelog`
* :doc:`dev_guidelines`

.. toctree::
   :maxdepth: 2
   :caption: Developer Documentation

   api
   changelog
   dev_guidelines
