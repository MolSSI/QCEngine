"""
Calls the Psi4 executable.
"""

import time
import sys

from . import config


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
    psiapi = config.get_config("psi_path")
    if (psiapi is not None) and (psiapi not in sys.path):
        sys.path.insert(1, psiapi)

    import psi4

    # Setup the job
    psi4.set_num_threads(config.get_config("cores_per_job"))
    json["memory"] = int(config.get_config("memory_per_job") * 1024 * 1024 * 1024 * 0.9)

    scratch = config.get_config("scratch_directory")
    if scratch is not None:
        json["scratch_location"] = scratch

    pversion = psi4.__version__ 

    if pversion == "1.1":
        # Check if RHF/UHF
        mol = psi4.geometry(json["molecule"])
        wfn = psi4.core.Wavefunction.build(mol, "def2-SVP")
        if wfn.molecule().multiplicity() != 1:
            json["options"]["reference"] = "uks"

        json["args"] = (json["method"], )

        # Compute!
        rjson = psi4.json_wrapper.run_json(json)
        psi4.core.clean()
        if rjson is False:
            json["success"] = False
            if "error" not in json:
                json["error"] = "Unspecified error occured"

    else:
        raise TypeError("Psi4 version '{}' not understood".format(pversion))

    if json["success"] is False:

        # Dispatch errors, PSIO Errors are not recoverable for future runs
        if "PSIO Error" in json["error"]:
            raise ValueError(json["error"])

    return json
