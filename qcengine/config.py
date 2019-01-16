"""
Creates globals for the qcengine module
"""

import copy
import fnmatch
import getpass
import logging
import os
import socket

from typing import Dict, Union

import cpuinfo
import psutil
import pydantic
import yaml

__all__ = ["get_config", "get_provenance", "global_repr", "NodeDescriptor"]

# Start a globals dictionary with small starting values
CPUINFO = cpuinfo.get_cpu_info()

# We want physical cores
if hasattr(psutil.Process(), "cpu_affinity"):
    cpu_cnt = len(psutil.Process().cpu_affinity())
else:
    cpu_cnt = psutil.cpu_count(logical=False)
    if cpu_cnt is None:
        cpu_cnt = psutil.cpu_count(logical=True)
CPUINFO["count"] = cpu_cnt

# Generic globals
GLOBALS = {}
GLOBALS["hostname"] = socket.gethostname()
GLOBALS["available_memory"] = round(psutil.virtual_memory().available / (1024**3), 3)
GLOBALS["cpu"] = CPUINFO["brand"]
GLOBALS["username"] = getpass.getuser()

# Handle logger
LOGGER = logging.getLogger("QCEngine")
LOGGER.setLevel(logging.CRITICAL)

NODE_DESCRIPTORS = {}


class NodeDescriptor(pydantic.BaseModel):
    """
    Description of an individual node
    """

    # Host data
    hostname_pattern: str
    name: str
    available_memory: int = None
    total_cores: int = None

    # Specifications
    jobs_per_node: int = 2
    memory_safety_factor: int = 10  # Percentage of memory as a safety factor
    memory_per_job: Union[int, str] = "auto"  # Amount of memory in Gb per node
    nthreads_per_job: Union[int, str] = "auto"  # Number of nthreads per job
    scratch_directory: str = None  # What location to use as scratch

    def __init__(self, **kwargs):
        """
        Initalize the pydantic class after processing options
        """

        super().__init__(**kwargs)

        for k, var in self.__values__.items():
            if isinstance(var, str) and var.startswith("$"):
                var = var.lstrip("$")
                if var in os.environ:
                    var = os.environ[var]
                else:
                    var = None
                self.__values__[k] = var

        # Pull from defaults
        if self.available_memory is None:
            self.available_memory = GLOBALS["available_memory"]

        if self.total_cores is None:
            self.total_cores = CPUINFO["count"]

        # Handle any auto keywords
        if self.nthreads_per_job == "auto":
            # Figure out number of threads per job
            self.nthreads_per_job = int(self.total_cores / self.jobs_per_node)
            if self.nthreads_per_job < 1:
                raise KeyError("Number of jobs per node exceeds the number of available cores.")

        if self.memory_per_job == "auto":
            memory_coeff = (1 - self.memory_safety_factor / 100)
            self.memory_per_job = round(self.available_memory * memory_coeff / self.jobs_per_node, 3)


class JobConfig(pydantic.BaseModel):

    # Specifications
    nthreads: int  # Number of nthreads per job
    memory: str  # Amount of memory in Gb per node
    scratch_directory: str  # What location to use as scratch


def _load_defaults():
    """
    Pulls the defaults from the QCA folder
    """

    # Find the config
    load_path = None
    test_paths = [os.getcwd(), os.path.join(os.path.expanduser('~'), ".qcarchive")]

    if "DQM_CONFIG_PATH" in os.environ:
        test_paths.insert(0, os.environ["DQM_CONFIG_PATH"])

    for path in test_paths:
        path = os.path.join(path, "qcengine.yaml")
        if os.path.exists(path):
            load_path = path
            break

    found_default = False
    if load_path is None:
        LOGGER.info("Could not find 'qcengine.yaml'. Searched the following paths: {}".format(", ".join(test_paths)))
        LOGGER.info("Using default options...")

    else:
        LOGGER.info("Found 'qcengine.yaml' at path: {}".format(load_path))
        with open(load_path, "r") as stream:
            user_config = yaml.load(stream)

        for k, v in user_config.items():
            NODE_DESCRIPTORS[k] = NodeDescriptor(name=k, **v)

    # Make sure we have a default
    if "default" not in NODE_DESCRIPTORS:
        NODE_DESCRIPTORS["default"] = NodeDescriptor(hostname_pattern="*", name="default")


# Pull in the local variables
_load_defaults()


def global_repr():
    """
    A representation of the current global configuration.
    """

    ret = ""
    ret += "Host information:\n"
    ret += '-' * 80 + "\n"
    prov = get_provenance()
    for k in ["username", "hostname", "cpu"]:
        ret += "{:<30} {:<30}\n".format(k, prov)
    ret += "{:<30} {:<30}\n".format("cores", CPUINFO["count"])

    ret += "\nJob information:\n"
    ret += '-' * 80 + "\n"
    for k, v in get_config().items():
        if (v is None) or isinstance(v, (str, int, float)):
            ret += "  {:<28} {}\n".format(k, v)
        else:
            raise TypeError("Global printing error, type '{}' not understood for key '{}'.".format(type(v), v))

    ret += '-' * 80 + "\n"

    return ret


def get_node_descriptor(hostname=None):
    """
    Find the correct NodeDescriptor based off current hostname
    """
    if hostname is None:
        hostname = GLOBALS["hostname"]

    # Find a match
    for name, node in NODE_DESCRIPTORS.items():
        if name == "default": continue

        if fnmatch.fnmatch(hostname, node.hostname_pattern):
            config = host_config
            break
    else:
        config = NODE_DESCRIPTORS["default"]

    return config


def get_config(hostname=None, local_options=None):
    """
    Returns the configuration key for qcengine.
    """

    node = get_node_descriptor(hostname)
    config = {
        "nthreads": node.nthreads_per_job,
        "memory": node.memory_per_job,
        "scratch_directory": node.scratch_directory
    }
    if local_options is not None:
        config.update(local_options)

    return JobConfig(**config)


def get_provenance():
    from qcengine import __version__
    ret = {
        "cpu": GLOBAL["cpu"],
        "hostname": GLOBAL["hostname"],
        "username": GLOBAL["username"],
        "qcengine_version": __version__
    }

    return ret


def get_logger():
    return LOGGER
