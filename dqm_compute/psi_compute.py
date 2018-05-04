"""
Calls the Psi4 executable.
"""

import time
import sys

from . import dqm_config


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

    # Insert API path if needed
    psiapi = dqm_config.get_config("psi_path")
    if (psiapi is not None) and (psiapi != sys.path[1]):
        sys.path.insert(1, psiapi)

    import psi4

    # Setup the job
    psi4.set_num_threads(dqm_config.get_config("cores_per_job"))
    json["memory"] = int(dqm_config.get_config("memory_per_job") * 1024 * 1024 * 1024 * 0.9)

    scratch = dqm_config.get_config("scratch_directory")
    if scratch is not None:
        json["scratch_location"] = scratch

    # Check if RHF/UHF
    mol = psi4.geometry(json["molecule"])
    wfn = psi4.core.Wavefunction.build(mol, "def2-SVP")
    if wfn.molecule().multiplicity() != 1:
        json["options"]["reference"] = "uks"

    # Compute!
    json = psi4.json_wrapper.run_json(json)
    psi4.core.clean()

    if json["success"] is False:

        # Dispatch errors, PSIO Errors are not recoverable for future runs
        if "PSIO Error" in json["error"]:
            raise ValueError(json["error"])

    return json
