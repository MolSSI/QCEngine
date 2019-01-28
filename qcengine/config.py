"""
Creates globals for the qcengine module
"""

import fnmatch
import getpass
import logging
import os
import socket
from typing import Optional

import cpuinfo
import psutil
import pydantic
import yaml

__all__ = ["get_config", "get_provenance_augments", "global_repr", "NodeDescriptor"]

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
GLOBALS["memory"] = round(psutil.virtual_memory().available / (1024**3), 3)
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
    scratch_directory: Optional[str] = None  # What location to use as scratch

    memory: Optional[float] = None
    memory_safety_factor: int = 10  # Percentage of memory as a safety factor

    # Specifications
    ncores: Optional[int] = None
    jobs_per_node: int = 2

    def __init__(self, **data):

        data = parse_environment(data)
        super().__init__(**data)

    class Config:
        ignore_extra = False


class JobConfig(pydantic.BaseModel):

    # Specifications
    ncores: int  # Number of ncores per job
    memory: float  # Amount of memory in GiB per node
    scratch_directory: Optional[str]  # What location to use as scratch

    class Config:
        ignore_extra = False


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

    if load_path is None:
        LOGGER.info("Could not find 'qcengine.yaml'. Searched the following paths: {}".format(", ".join(test_paths)))
        LOGGER.info("Using default options...")

    else:
        LOGGER.info("Found 'qcengine.yaml' at path: {}".format(load_path))
        with open(load_path, "r") as stream:
            user_config = yaml.load(stream)

        for k, v in user_config.items():
            NODE_DESCRIPTORS[k] = NodeDescriptor(name=k, **v)


# Pull in the local variables
_load_defaults()


def global_repr():
    """
    A representation of the current global configuration.
    """

    ret = ""
    ret += "Host information:\n"
    ret += "-" * 80 + "\n"

    prov = get_provenance_augments()
    for k in ["username", "hostname", "cpu"]:
        ret += "{:<30} {:<30}\n".format(k, prov[k])

    ret += "\nNode information:\n"
    ret += "-" * 80 + "\n"
    for k, v in get_node_descriptor():
        ret += "  {:<28} {}\n".format(k, v)

        if k in ["scratch_directory", "memory_per_job"]:
            ret += "\n"

    ret += "\nJob information:\n"
    ret += "-" * 80 + "\n"
    for k, v in get_config():
        ret += "  {:<28} {}\n".format(k, v)

    ret += "-" * 80 + "\n"

    return ret


def get_node_descriptor(hostname=None):
    """
    Find the correct NodeDescriptor based off current hostname
    """
    if isinstance(hostname, NodeDescriptor):
        return hostname

    if hostname is None:
        hostname = GLOBALS["hostname"]

    # Find a match
    for name, node in NODE_DESCRIPTORS.items():

        if fnmatch.fnmatch(hostname, node.hostname_pattern):
            config = node
            break
    else:
        config = NodeDescriptor(
            name="default", hostname_pattern="*", memory=GLOBALS["memory"], ncores=CPUINFO["count"])

    return config


def parse_environment(data):
    """
    Parses a dictionary looking for environmental variables
    """
    ret = {}
    for k, var in data.items():
        if isinstance(var, str) and var.startswith("$"):
            var = var.replace("$", "", 1)
            if var in os.environ:
                var = os.environ[var]
            else:
                var = None

        ret[k] = var

    return ret


def get_config(*, hostname=None, local_options=None):
    """
    Returns the configuration key for qcengine.
    """

    if local_options is None:
        local_options = {}

    local_options = parse_environment(local_options)

    # Node data
    node = get_node_descriptor(hostname)
    ncores = node.ncores or CPUINFO["count"]
    scratch_directory = local_options.get("scratch_directory", None) or node.scratch_directory

    # Jobs per node
    jobs_per_node = local_options.pop("jobs_per_node", None) or node.jobs_per_node

    # Handle memory
    memory = local_options.pop("memory", None)
    if memory is None:
        memory = node.memory or GLOBALS["memory"]
        memory_coeff = (1 - node.memory_safety_factor / 100)
        memory = round(memory * memory_coeff / jobs_per_node, 3)

    # Handle ncores
    ncores = local_options.pop("ncores", None) or int(ncores / jobs_per_node)
    if ncores < 1:
        raise KeyError("Number of jobs per node exceeds the number of available cores.")

    config = {"ncores": ncores, "memory": memory, "scratch_directory": scratch_directory}
    if local_options is not None:
        config.update(local_options)

    return JobConfig(**config)


def get_provenance_augments():
    from qcengine import __version__
    return dict(
        cpu=GLOBALS["cpu"],
        hostname=GLOBALS["hostname"],
        username=GLOBALS["username"],
        qcengine_version=__version__,
    )


def get_logger():
    return LOGGER
