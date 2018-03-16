"""
Calls the Psi4 executable.
"""

from .dqm_config import get_config
import time
import sys


def test_pass():
    """
    Runs a quick Psi4 test
    """

    json_data = {"something": "stuff"}

    return run_pass(json_data)


def run_pass(json):
    """
    Runs Psi4 in API mode
    """

    t = time.time()
    json["memory"] = int(get_config("memory_per_job") * 1024 * 1024 * 1024 * 0.9)

    json["success"] = True
    json["wall_time"] = time.time() - t
    return json
