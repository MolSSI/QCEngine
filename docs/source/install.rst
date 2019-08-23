Install QCEngine
================

You can install qcengine with ``conda`` or with ``pip``.

Conda
-----

You can install qcengine using `conda <https://www.anaconda.com/download/>`_:

.. code-block:: console

    >>> conda install qcengine -c conda-forge

This installs QCEngine and its dependencies. The qcengine package is maintained on the
`conda-forge channel <https://conda-forge.github.io/>`_.


Pip
---

You can also install QCEngine using ``pip``:

.. code-block:: console

   >>> pip install qcengine


Test the Installation
---------------------

.. note::

   QCEngine is a wrapper for other quantum chemistry codes. The tests for QCEngine will only test the wrapper for a
   given code if its detected in the ``$PATH`` or current Python Environment, otherwise the tests for that package are
   skipped. Keep this in mind if you see many ``skip`` or ``s`` codes output from PyTest.

You can test to make sure that Engine is installed correctly by first installing ``pytest``.

From ``conda``:

.. code-block:: console

   >>> conda install pytest -c conda-forge

From ``pip``:

.. code-block:: console

   >>> pip install pytest

Then, run the following command:

.. code-block::

   >>> pytest --pyargs qcengine


Developing from Source
----------------------

If you are a developer and want to make contributions Engine, you can access the source code from
`github <https://github.com/molssi/qcengine>`_.
