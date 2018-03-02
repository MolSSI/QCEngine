"""
Creates globals for the dqm_compute module
"""

import os
import socket
import fnmatch
import yaml
import json

__all__ = ["get_config"]

# Start a globals dictionary with small starting values
_globals = {}

_globals["hostname"] = socket.gethostname()
_globals["default_compute"] = {
    "psi_path": None,  # Path for the Psi4 API
    "jobs_per_node": 1,  # Number of jobs per node
    "cores_per_job": 1,  # Number of cores per job
    "memory_per_job": 2,  # Amount of memory in Gb per node
    "scratch_directory": None, # What location to use as scratch
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

def _load_locals():

    # Find the dqm_config
    load_path = None
    test_paths = [os.getcwd(), os.path.join(os.path.expanduser('~'), ".dqm")]

    if "FW_CONFIG_FILE" in os.environ:
        test_paths.insert(0, os.environ["FW_CONFIG_FILE"])

    for path in test_paths:
        path = os.path.join(path, "dqm_config.yaml")
        if os.path.exists(path):
            load_path = path
            break

    if load_path is None:
        raise OSError("Could not find 'dqm_config.yaml'. Search the following paths: %s" % ", ".join(test_paths))

    # Load the library
    with open(load_path) as stream:
        user_config = yaml.load(stream)

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


# Pull in the local variables
_load_locals()

def get_hostname():
    """
    Returns the global current hostname
    """

    return _globals["hostname"]

def get_config(key=None, hostname=None):
    """
    Returns the configuration key for dqm_compute.
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
            raise Exception("Key %s asked for, but not in local data")
        return config[key]
