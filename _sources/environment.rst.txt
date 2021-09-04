Environment Detection
======================

QCEngine can inspect the current compute environment to determine the resources available to it.

Node Description
----------------

QCEngine can detect node descriptions to obtain general information about the current node.

.. code:: python

    >>> qcng.config.get_node_descriptor()
    <NodeDescriptor hostname_pattern='*' name='default' scratch_directory=None
                    memory=5.568 memory_safety_factor=10 ncores=4 jobs_per_node=2>

Config
------

The configuration file operated based on the current node descriptor and can be overridden:

.. code:: python

    >>> qcng.get_config()
    <JobConfig ncores=2 memory=2.506 scratch_directory=None>

    >>> qcng.get_config(local_options={"scratch_directory": "/tmp"})
    <JobConfig ncores=2 memory=2.506 scratch_directory='/tmp'>

    >>> os.environ["SCRATCH"] = "/my_scratch"
    >>> qcng.get_config(local_options={"scratch_directory": "$SCRATCH"})
    <JobConfig ncores=2 memory=2.506 scratch_directory='/my_scratch'>

Global Environment
-------------------

The global environment can also be inspected directly.

.. code:: python

    >>> qcng.config.get_global()
    {
        'hostname': 'qcarchive.molssi.org',
        'memory': 5.568,
        'username': 'user',
        'ncores': 4,
        'cpuinfo': {
            'python_version': '3.6.7.final.0 (64 bit)',
            'brand': 'Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz',
            'hz_advertised': '2.9000 GHz',
            ...
        },
        'cpu_brand': 'Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz'
    }

Configuration Files
-------------------

The computational environment defaults can be overridden by configuration files.

Configuration files must be named ``qcengine.yaml`` and stored either in the directory
from which you run QCEngine, a folder named ``.qcarchive`` in your home directory,
or in a folder specified by the ``DQM_CONFIG_PATH`` environmental variable.
Only one configuration file will be used if multiple are available.
The ``DQM_CONFIG_PATH`` configuration file takes precedence over the current directory,
which takes precedence over the ``.qcarchive`` folder.

The configuration file is a YAML file that contains a dictionary of different node configurations.
The keys in the YAML file are human-friendly names for the configurations.
The values are dictionaries that define configurations for different nodes,
following the ``NodeDescription`` schema:

.. autoclass:: qcengine.config.NodeDescriptor

When running QCEngine, the proper configuration for a node is determined based on the hostname of the node
and matching the ``hostname_pattern`` to each of the configurations defined in ``qcengine.yaml``.

An example ``qcengine.yaml`` file that sets the scratch directory for all nodes is as follows:

.. code:: yaml

    all:
      hostname_pattern: "*"
      scratch_directory: ./scratch

Cluster Configuration
---------------------

A node configuration file is required when using node-parallel tasks on a compute cluster.
The configuration file must contain a description of the command used to launch MPI tasks and,
in some cases, the designation that a certain node is a compute node.
See the descriptions for ``mpiexec_command`` and ``is_batch_node`` in the ``NodeDescriptor``
documentation for further details.
