Semiempirical Quantum Mechanics
===============================

For semiempirical quantum mechanics (SQM) engines to fit the AtomicInput/Result schema the following convention is used:

- Method: The unique method name (PM7 or GFN2-xTB) including the parametrisation information is provided, no basis is needed

As for quantum mechanical methods a minimal Molecule object is sufficient as input.

.. note::

   Semiemprical engines might not handle the concept of ghost atoms correctly, check carefully how the used engine handles ghost atoms.
   To be sure remove ghost atoms from input to semiempirical engines beforehand.


Example
-------

For example, running a calculation with the GFN2-xTB method using the ``xtb`` engine would work like any other QM engine with

.. code:: python

   >>> import qcelemental as qcel
   >>> mol = qcel.models.Molecule(
   ...     symbols=["O", "H", "H"],
   ...     geometry=[
   ...         [ 0.00000000000000, 0.00000000000000,-0.73578586109551],
   ...         [ 1.44183152868459, 0.00000000000000, 0.36789293054775],
   ...         [-1.44183152868459, 0.00000000000000, 0.36789293054775],
   ...     ],
   ... )
   ...
   >>> model = qcel.models.AtomicInput(
   ...     molecule=mol,
   ...     driver="energy",
   ...     model={"method": "GFN2-xTB"},
   ... )
   ...
   >>> import qcengine as qcng
   >>> ret = qcng.compute(model, "xtb")
   >>> ret.return_result
   -5.070451354836705


MOPAC
-----

The following semiempirical Hamiltonians are supported with the MOPAC engine.

============= ===========
 Method        Basis
============= ===========
 mndo          None
 am1           None
 pm3           None
 rm1           None
 mndod         None
 pm6           None
 pm6-d3        None
 pm6-dh+       None
 pm6-dh2       None
 pm6-dh2x      None
 pm6-d3h4      None
 pm6-3dh4x     None
 pm7           None
 pm7-ts        None
============= ===========


xtb
---

The following extended tight binding Hamiltonians are available with the ``xtb`` engine.

=========== ======== ==============================
 Method      Basis    Reference
=========== ======== ==============================
 GFN2-xTB    None     10.1021/acs.jctc.8b01176
 GFN1-xTB    None     10.1021/acs.jctc.7b00118
 GFN0-xTB    None     10.26434/chemrxiv.8326202.v1
=========== ======== ==============================
