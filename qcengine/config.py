"""
Creates globals for the qcengine module
"""

import copy
import fnmatch
import getpass
import os
import socket

import cpuinfo
import psutil
import yaml

__all__ = ["get_global", "get_config", "get_provenance", "global_repr"]

# Start a globals dictionary with small starting values
_globals = {}
_cpuinfo = cpuinfo.get_cpu_info()

# We want physical cores
_cpuinfo["count"] = psutil.cpu_count(logical=False)

_globals["hostname"] = socket.gethostname()
_globals["cpu"] = _cpuinfo["brand"]
_globals["username"] = getpass.getuser()
_globals["default_compute"] = {

    # Program paths
    "psi_path": None,  # Path for the Psi4 API

    # Specifications
    "jobs_per_node": 2,  # Number of jobs per node
    "nthreads_per_job": 'auto',  # Number of nthreads per job
    "memory_per_job": 'auto',  # Amount of memory in Gb per node
    "scratch_directory": None,  # What location to use as scratch
}
_globals["other_compute"] = {}


def _process_variables(var):
    # Environmental var
    if isinstance(var, str) and var.startswith("$"):
        var = var.lstrip("$")
        if var in os.environ:
            return os.environ[var]
        else:
            return None

    # Normal var
    else:
        return var


def _process_autos(data):
    if data.get("nthreads_per_job", False) == "auto":
        nthreads = _cpuinfo["count"]

        # Figure out number of threads per job 
        data["nthreads_per_job"] = int(nthreads / data["jobs_per_node"])
        if data["nthreads_per_job"] < 1:
            raise KeyError("Number of jobs per node exceeds the number of available cores.")

    if data.get("memory_per_job", False) == "auto":
        data["memory_per_job"] = round(psutil.virtual_memory().available * 0.9 / data["jobs_per_node"] / (1024**3), 3)


def load_options(load_path):
    """
    Options can be loaded from a specific path
    """

    # Load the library
    with open(load_path) as stream:
        user_config = yaml.load(stream)

    _globals["config_path"] = load_path

    # Override default keys
    default_keys = list(_globals["default_compute"].keys())

    if "default_compute" in user_config:
        for k, v in user_config["default_compute"].items():
            if k not in default_keys:
                raise KeyError("Key %s not accepted for default_compute" % k)
            _globals["default_compute"][k] = _process_variables(v)

    default_keys.append("hostname")

    if "other_compute" in user_config:
        for host, config in user_config["other_compute"].items():
            _globals["other_compute"][host] = _globals["default_compute"].copy()

            if "hostname" not in config:
                raise KeyError("Other_compute must have a hostname to help identify the server")
            for k, v in config.items():
                if k not in default_keys:
                    raise KeyError("Key %s not accepted for default_compute" % k)
                _globals["other_compute"][host][k] = _process_variables(v)

    # Process autos
    _process_autos(_globals["default_compute"])
    for k, v in _globals["other_compute"].items():
        _process_autos(v)


def _load_locals():

    # Find the config
    load_path = None
    test_paths = [os.getcwd(), os.path.join(os.path.expanduser('~'), ".qc")]

    if "DQM_CONFIG_PATH" in os.environ:
        test_paths.insert(0, os.environ["DQM_CONFIG_PATH"])

    for path in test_paths:
        path = os.path.join(path, "qcengine_config.yaml")
        if os.path.exists(path):
            load_path = path
            break

    if load_path is None:
        print("Could not find 'qcengine_config.yaml'. Searched the following paths: %s" % ", ".join(test_paths))
        print("Using default options...")

        # Process autos
        _process_autos(_globals["default_compute"])
        for k, v in _globals["other_compute"].items():
            _process_autos(v)
    else:
        load_options(load_path)



# Pull in the local variables
_load_locals()


def get_hostname():
    """
    Returns the global current hostname
    """

    return _globals["hostname"]


def global_repr():

    ret = ""
    ret += "Host information:\n"
    ret += '-' * 80 + "\n"
    for k in ["username", "hostname", "cpu"]:
        ret += "{:<30} {:30}\n".format(k, get_global(k))

    ret += "\nJob information:\n"
    ret += '-' * 80 + "\n"
    for k, v in get_config().items():
        if (v is None) or isinstance(v, (str, int, float)):
            ret += "  {:<28} {}\n".format(k, v)
        else:
            raise TypeError("Global printing error, type '{}' not understood for key '{}'.".format(type(v), v))

    ret += '-' * 80 + "\n"

    return ret


def get_global(name):
    return copy.deepcopy(_globals[name])


def get_config(key=None, hostname=None):
    """
    Returns the configuration key for qcengine.
    """
    config = None
    hostname_match = None
    if hostname is None:
        hostname = _globals["hostname"]

    # Find a match
    for host, config in _globals["other_compute"].items():
        if fnmatch.fnmatch(hostname, config["hostname"]):
            hostname_match = host
            config = config
            break

    # Use default
    if hostname_match is None:
        config = _globals["default_compute"]

    if key is None:
        return config.copy()
    else:
        if key not in config:
            raise Exception("Key '{}' asked for, but not in local data".format(key))
        return config[key]


def get_provenance():
    ret = {"cpu": get_global("cpu"), "hostname": get_global("hostname"), "username": get_global("username")}

    return ret
