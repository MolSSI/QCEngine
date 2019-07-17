"""
Creates globals for the qcengine module
"""

import fnmatch
import getpass
import logging
import os
import socket
from typing import Any, Dict, Optional, Union

import pydantic

from .extras import get_information

__all__ = ["get_config", "get_provenance_augments", "global_repr", "NodeDescriptor"]

# Start a globals dictionary with small starting values
_global_values = None
NODE_DESCRIPTORS = {}
LOGGER = logging.getLogger("QCEngine")
LOGGER.setLevel(logging.CRITICAL)


# Generic globals
def get_global(key: Optional[str] = None) -> Union[str, Dict[str, Any]]:
    import cpuinfo
    import psutil
    global _global_values
    if _global_values is None:
        _global_values = {}
        _global_values["hostname"] = socket.gethostname()
        _global_values["memory"] = round(psutil.virtual_memory().available / (1024**3), 3)
        _global_values["username"] = getpass.getuser()

        # Work through VMs and logical cores.
        if hasattr(psutil.Process(), "cpu_affinity"):
            cpu_cnt = len(psutil.Process().cpu_affinity())
        else:
            cpu_cnt = psutil.cpu_count(logical=False)
            if cpu_cnt is None:
                cpu_cnt = psutil.cpu_count(logical=True)

        _global_values["ncores"] = cpu_cnt

        _global_values["cpuinfo"] = cpuinfo.get_cpu_info()
        _global_values["cpu_brand"] = _global_values["cpuinfo"]["brand"]

    if key is None:
        return _global_values.copy()
    else:
        return _global_values[key]


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
    retries: int = 0

    def __init__(self, **data: Dict[str, Any]) -> 'BaseModel':

        data = parse_environment(data)
        super().__init__(**data)

    class Config:
        extra = "forbid"


class JobConfig(pydantic.BaseModel):

    # Specifications
    ncores: int  # Number of ncores per job
    memory: float  # Amount of memory in GiB per node
    scratch_directory: Optional[str]  # What location to use as scratch
    retries: int # Number of retries on random failures

    class Config:
        extra = "forbid"


def _load_defaults() -> None:
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
        import yaml
        LOGGER.info("Found 'qcengine.yaml' at path: {}".format(load_path))
        with open(load_path, "r") as stream:
            user_config = yaml.load(stream)

        for k, v in user_config.items():
            NODE_DESCRIPTORS[k] = NodeDescriptor(name=k, **v)


# Pull in the local variables
_load_defaults()


def global_repr() -> str:
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


def get_node_descriptor(hostname: Optional[str] = None) -> NodeDescriptor:
    """
    Find the correct NodeDescriptor based off current hostname
    """
    if isinstance(hostname, NodeDescriptor):
        return hostname

    if hostname is None:
        hostname = get_global("hostname")

    # Find a match
    for name, node in NODE_DESCRIPTORS.items():

        if fnmatch.fnmatch(hostname, node.hostname_pattern):
            config = node
            break
    else:
        config = NodeDescriptor(name="default",
                                hostname_pattern="*",
                                memory=get_global("memory"),
                                ncores=get_global("ncores"))

    return config


def parse_environment(data: Dict[str, Any]) -> Dict[str, Any]:
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


def get_config(*, hostname: Optional[str] = None, local_options: Dict[str, Any] = None) -> JobConfig:
    """
    Returns the configuration key for qcengine.
    """

    if local_options is None:
        local_options = {}

    local_options = parse_environment(local_options)
    config = {}

    # Node data
    node = get_node_descriptor(hostname)
    ncores = node.ncores or get_global("ncores")
    config["scratch_directory"] = local_options.pop("scratch_directory", node.scratch_directory)
    config["retries"] = local_options.pop("retries", node.retries)

    # Jobs per node
    jobs_per_node = local_options.pop("jobs_per_node", None) or node.jobs_per_node

    # Handle memory
    memory = local_options.pop("memory", None)
    if memory is None:
        memory = node.memory or get_global("memory")
        memory_coeff = (1 - node.memory_safety_factor / 100)
        memory = round(memory * memory_coeff / jobs_per_node, 3)

    config["memory"] = memory

    # Handle ncores
    ncores = local_options.pop("ncores", int(ncores / jobs_per_node))
    if ncores < 1:
        raise KeyError("Number of jobs per node exceeds the number of available cores.")

    config["ncores"] = ncores

    if local_options is not None:
        config.update(local_options)

    return JobConfig(**config)


def get_provenance_augments() -> Dict[str, str]:
    return {
        "cpu": get_global("cpu_brand"),
        "hostname": get_global("hostname"),
        "username": get_global("username"),
        "qcengine_version": get_information("version")
    }


def get_logger() -> 'Logger':
    return LOGGER
