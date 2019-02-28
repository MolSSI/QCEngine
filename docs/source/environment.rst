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
