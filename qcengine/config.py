"""
Creates globals for the qcengine module
"""

import copy
import fnmatch
import getpass
import logging
import os
import socket

from typing import Dict

import cpuinfo
import psutil
import pydantic
import yaml

__all__ = ["get_global", "get_config", "get_provenance", "global_repr"]

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
GLOBALS["available_memory"] = psutil.virtual_memory().available
GLOBALS["cpu"] = CPUINFO["brand"]
GLOBALS["username"] = getpass.getuser()

# User file configuration
USER_CONFIG = {}

# Handle logger
LOGGER = logging.getLogger("QCEngine")
LOGGER.setLevel(logging.CRITICAL)


def _process_options(data):
    """
    Expands environmental variables and sets automatic keywords
    """
    for k, var in data.items():
        if isinstance(var, str) and var.startswith("$"):
            var = var.lstrip("$")
            if var in os.environ:
                var = os.environ[var]
            else:
                var = None
            data[k] = var

    if data.get("nthreads_per_job", False) == "auto":
        # Figure out number of threads per job
        data["nthreads"] = int(CPUINFO["count"] / data["jobs_per_node"])
        if data["nthreads"] < 1:
            raise KeyError("Number of jobs per node exceeds the number of available cores.")

    if data.get("memory", False) == "auto":
        data["memory"] = round(data["available_memory"] * 0.9 / data["jobs_per_node"] / (1024**3), 3)


class Options(pydantic.BaseModel):

    # Host data
    hostname: str
    available_memory: str
    cpu: str
    username: str

    # Program paths
    psi_path: None
    rdkit_path: None

    # Specifications
    nthreads: int = None  # Number of nthreads per job
    memory: str = None  # Amount of memory in Gb per node
    scratch_directory: None  # What location to use as scratch

    def __init__(self, **kwargs):
        """
        Initalize the pydantic class after processing options
        """

        _process_options(kwargs)
        kwargs.update(GLOBALS)

        super().__init__(**kwargs)


def _load_defaults():
    """
    Pulls the defaults from the QCA folder
    """

    # Find the config
    load_path = None
    test_paths = [os.getcwd(), os.path.join(os.path.expanduser('~'), ".qca")]

    if "DQM_CONFIG_PATH" in os.environ:
        test_paths.insert(0, os.environ["DQM_CONFIG_PATH"])

    for path in test_paths:
        path = os.path.join(path, "qcengine_config.yaml")
        if os.path.exists(path):
            load_path = path
            break

    if load_path is None:
        LOGGER.info("Could not find 'qcengine_config.yaml'. Searched the following paths: %s" % ", ".join(test_paths))
        LOGGER.info("Using default options...")

    else:
        if isinstance(load_path, str):
            with open(load_path, "r") as stream:
                USER_CONFIG = yaml.load(stream)
        elif isinstance(load_path, dict):
            USER_CONFIG = load_path
        else:
            raise TypeError("Unknown options load")


# Pull in the local variables
_load_defaults()


def get_hostname():
    """
    Returns the global current hostname
    """

    return GLOBALS["hostname"]


def global_repr():

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


def get_config(hostname=None, local_options=None):
    """
    Returns the configuration key for qcengine.
    """
    config = None
    hostname_match = None
    if hostname is None:
        hostname = GLOBALS["hostname"]

    # Find a match
    for host, host_config in GLOBALS["other_compute"].items():
        if fnmatch.fnmatch(hostname, host_config["hostname"]):
            hostname_match = host
            config = host_config
            break

    local_config = {}

    # Use default
    if hostname_match is None:
        config = GLOBALS["default_compute"]

    return Options(**config)

def get_provenance():
    ret = {"cpu": GLOBAL["cpu"], "hostname": GLOBAL["hostname"], "username": GLOBAL["username"]}

    return ret


def get_logger():
    return LOGGER
