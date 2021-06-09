Molecular Mechanics
===================

For Molecular Mechanics (MM) engines to fit the AtomicInput/Result schema
the following convention is used:

- Method: The force field used such as MMFF94, GAFF, OpenFF-1.0.0.
- Basis: The typing engine used to find the required paramters.

For all MM computations the input Molecule object must have connectivity and this will not be automatically assigned for you.

Example
-------

.. code:: python

    >>> mol = qcel.models.Molecule(
    >>>     symbols=["O", "H", "H"],
    >>>     geometry=[[0, 0, 0], [0, 0, 2], [0, 2, 0]],
    >>>     connectivity=[[0, 1, 1], [0, 2, 1]],
    >>> )

    >>> model = qcel.models.AtomicInput(
    >>>     molecule=mol,
    >>>     driver="energy",
    >>>     model={"method": "openff-1.0.0", "basis": "smirnoff"},
    >>> )
    >>> ret = qcng.compute(model, "openmm")
    >>> ret.return_result
    0.011185654397410195



OpenMM
------

Currently OpenMM only supports the smirnoff typing engine from the
``openff-toolkit``. Currently available force fields are the following:

+----------------------------+------------+
| Method                     | Basis      |
+============================+============+
| smirnoff99Frosst-1.1.0     | smirnoff   |
+----------------------------+------------+
| openff-1.0.0               | smirnoff   |
+----------------------------+------------+
| openff_unconstrained-1.0.0 | smirnoff   |
+----------------------------+------------+

Other forcefields may be available depending on your version of the ``openff-toolkit``, see their `docs <https://open-forcefield-toolkit.readthedocs.io>`_ for more information.

RDKit
-----

RDKit force fields currently do not require a typing engine and the basis is omitted in all computations. Currently available force fields are the following:

+----------------------------+------------+
| Method                     | Basis      |
+============================+============+
| UFF                        | None       |
+----------------------------+------------+
| MMFF94                     | None       |
+----------------------------+------------+
| MMFF94s                    | None       |
+----------------------------+------------+


xtb
---

Experimental access to force fields are available with the ``xtb`` engine.
Note that the ``xtb`` engine will not require nor use a topology information provided in the input schema.

=========== ======== ==============================
 Method      Basis    Reference
=========== ======== ==============================
 GFN-FF      None     10.1002/anie.202004239
=========== ======== ==============================
