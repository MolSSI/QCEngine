Installing QCEngine
===================

You can install qcengine with ``conda``, with ``pip``, or by installing from source.

Conda
-----

You can update qcengine using `conda <https://www.anaconda.com/download/>`_::

    conda install qcengine -c conda-forge

This installs qcengine and the NumPy dependancy.

The qcengine package is maintained on the
`conda-forge channel <https://conda-forge.github.io/>`_.


Pip
---

To install qcengine with ``pip`` there are a few options, depending on which
dependencies you would like to keep up to date:

*   ``pip install qcengine``

Install from Source
-------------------

To install qcengine from source, clone the repository from `github
<https://github.com/molssi/qcengine>`_::

    git clone https://github.com/molssi/qcengine.git
    cd qcengine
    python setup.py install

or use ``pip`` for a local install::

    pip install -e .


Test
----

Test qcengine with ``py.test``::

    cd qcengine
    py.test
