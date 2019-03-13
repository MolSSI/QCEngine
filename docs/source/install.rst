Installing QCEngine
===================

You can install ``qcengine`` with ``conda``, with ``pip``, or by installing from source.

Conda
-----

You can update ``qcengine`` using `conda <https://www.anaconda.com/download/>`_::

    conda install qcengine -c conda-forge

This command installs qcengine and its dependencies. The ``qcengine`` package is maintained on the
`conda-forge channel <https://conda-forge.github.io/>`_.


Pip
---

To install ``qcengine`` with ``pip`` there are a few options, depending on which
dependencies you would like to keep up to date::

    pip install qcengine

Install from Source
-------------------

To install qcengine from source, clone the repository from `github
<https://github.com/molssi/qcengine>`_::

    git clone https://github.com/MolSSI/QCEngine.git
    cd qcengine
    python setup.py install

or use ``pip`` for a local install::

    pip install -e .

We recommend building a development environment with the following lines::

    cd qcengine
    python devtools/scripts/conda_env.py -n=qcngdev -p=3.6 devtools/conda-envs/psi.yaml
    conda activate qcngdev

This will build out a new environment with several compute backends for
``qcengine`` which provides a platform to test and develop the code.


Test
----

Test a ``qcengine`` local install with ``pytest``::

    cd qcengine
    pytest -v
