"""
Calls the Psi4 executable.
"""

from .dqm_config import get_config
import time
import sys


def test_psi4():
    """
    Runs a quick Psi4 test
    """

    json_data = {}
    json_data["molecule"] = """He 0 0 0\n--\nHe 0 0 1"""
    json_data["driver"] = "energy"
    json_data["method"] = 'SCF'
    json_data["options"] = {"BASIS": "STO-3G"}
    json_data["return_output"] = False

    return run_psi4(json_data)


def run_psi4(json):
    """
    Runs Psi4 in API mode
    """
    t = time.time()

    psiapi = get_config("psi_path")
    if psiapi is not None:
        sys.path.insert(1, psiapi)

    import psi4

    # Setup the job
    psi4.set_num_threads(get_config("cores_per_job"))
    json["memory"] = int(get_config("memory_per_job") * 1024 * 1024 * 1024 * 0.9)

    scratch = get_config("scratch_directory")
    if scratch is not None:
        json["scratch_location"] = scratch

    # Compute!
    json = psi4.json_wrapper.run_json(json)

    json["wall_time"] = time.time() - t
    return json
