Command Line Interface
======================

QCEngine provides a command line interface with three commands:

* ``qcengine info`` displays information about the environment detected by QCEngine.
* ``qcengine run`` runs a program.
* ``qcengine run-procedure`` runs a procedure.

Info Command
------------

Command Invocation
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    qcengine info <options>

Command Description
~~~~~~~~~~~~~~~~~~~

This command prints information about the QCEngine environment.

Arguments
~~~~~~~~~

``category``
    The information categories to show. Choices include:

    * ``version``: Print version of QCEngine and QCElemental.
    * ``programs``: Print detected and supported programs.
    * ``procedures``: Print detected and supported procedures.
    * ``config``: Print host, compute, and job configuration
    * ``all``: Print all available information.

    By default, all available information is printed.

Run Command
-----------

Command Invocation
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    qcengine run <program> <data>

Command Description
~~~~~~~~~~~~~~~~~~~

This command runs a program on a given task and outputs the result as a JSON blob.

Arguments
~~~~~~~~~

``program``
    The program to run.

``data``
    Data describing the task. One of:

    * A JSON blob.
    * A file name.
    * '-', indicating data will be read from STDIN.


Run-Procedure Command
---------------------

Command Invocation
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    qcengine run-procedure <program> <data>

Command Description
~~~~~~~~~~~~~~~~~~~

This command runs a procedure on a given task and outputs the result as a JSON blob.

Arguments
~~~~~~~~~

``procedure``
    The procedure to run.

``data``
    Data describing the task. One of:

    * A JSON blob.
    * A file name.
    * '-', indicating data will be read from STDIN.