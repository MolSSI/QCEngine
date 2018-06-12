"""
Calls the Psi4 executable.
"""

import sys

from pkg_resources import parse_version

from . import config


def _build_v11_mol(jm):
    """
    Converts a QC JSON molecule to the Psi4 expected string for Psi4 version 1.1.
    """

    # Cannot handle fragments
    if "fragments" in jm:
        return False, "Psi4 v1.1 cannot support fragments"

    if "real" in jm:
        real = jm["real"]
    else:
        real = [True for x in range(len(jm["symbols"]))]

    # Charge and multiplicity, skip otherwise
    psimol = ""
    if ("molecular_charge" in jm) or ("molecular_multiplicity" in jm):
        psimol += "\n    "
        if "molecular_charge" in jm:
            psimol += "%d  " % jm["molecular_charge"]
        else:
            psimol += "0 "

        if "molecular_multiplicity" in jm:
            psimol += "%d" % jm["molecular_multiplicity"]
        else:
            psimol += "1"

    psimol += "\n"

    for x in range(len(jm["symbols"])):
        shift = x * 3
        if real[x]:
            psimol += "   %-5s" % jm["symbols"][x]
        else:
            psimol += "Gh(%-5s)" % jm["symbols"][x]

        psimol += "    % 14.10f % 14.10f % 14.10f\n" % (jm["geometry"][shift], jm["geometry"][shift + 1],
                                                        jm["geometry"][shift + 2])

    psimol += "   units bohr"
    if "fix_com" in jm and jm["fix_com"]:
        psimol += "\n   no_com"

    if "fix_orientation" in jm and jm["fix_orientation"]:
        psimol += "\n   no_reorient"

    return psimol


def _parse_psi_version(version):
    if "undef" in version:
        raise TypeError(
            "Using custom build Psi4 without tags. Please `git pull origin master --tags` and recompile Psi4.")

    return parse_version(version)


def run_psi4(input_data):
    """
    Runs Psi4 in API mode
    """

    # Insert API path if needed
    psiapi = config.get_config("psi_path")
    if (psiapi is not None) and (psiapi not in sys.path):
        sys.path.insert(1, psiapi)

    try:
        import psi4
    except ImportError:
        raise ImportError("Could not find Psi4 in the Python path.")

    # Setup the job
    input_data["nthreads"] = config.get_config("nthreads_per_job")
    input_data["memory"] = int(config.get_config("memory_per_job") * 1024 * 1024 * 1024 * 0.9)
    input_data["success"] = False

    scratch = config.get_config("scratch_directory")
    if scratch is not None:
        input_data["scratch_location"] = scratch

    psi_version = _parse_psi_version(psi4.__version__)

    if psi_version == parse_version("1.1"):

        json_mol = input_data["molecule"]
        mol_str = _build_v11_mol(json_mol)

        input_data["options"] = input_data["keywords"]
        input_data["options"]["BASIS"] = input_data["model"]["basis"]
        psi4.set_num_threads(input_data["nthreads"], quiet=True)

        # Check if RHF/UHF
        mol = psi4.geometry(mol_str)
        wfn = psi4.core.Wavefunction.build(mol, "def2-SVP")
        if wfn.molecule().multiplicity() != 1:
            input_data["options"]["reference"] = "uks"

        input_data["args"] = (input_data["model"]["method"], )

        # v1.1 wanted an actual string
        input_data["molecule"] = mol_str

        # Compute!
        ret_data = psi4.json_wrapper.run_json(input_data)
        psi4.core.clean()
        if ret_data is False:
            ret_data["success"] = False
            if "error" not in rjson:
                ret_data["error"] = "Unspecified error occured."

        ret_data["molecule"] = json_mol

        # Manuall add in some fields
        ret_data["properties"] = {"return_energy": ret_data["variables"]["CURRENT ENERGY"]}
        del ret_data["variables"]

        if isinstance(ret_data["return_value"], float):
            ret_data["return_result"]
        else:
            import numpy as np  # Will have this if using Psi4
            arr = np.fromstring(ret_data["return_value"]["data"][0], dtype=np.double)
            ret_data["return_result"] = arr.ravel().tolist()

        del ret_data["return_value"]

        ret_data["provenance"] = {
            "version": "1.1",
            "routine": "psi4.json.run_json",
            "creator": "Psi4",
        }

    elif psi_version > parse_version("1.2rc2.dev500"):

        # Handle slight RC2 weirdness
        psi_12rc2_tweak = (psi_version == parse_version("1.2rc2"))
        psi_12rc2_tweak &= psi4.metadata.version_formatter("") != '(inplace)'

        if psi_12rc2_tweak:
            input_data["schema_name"] = "QC_JSON"
            input_data["schema_version"] = 0
            psi4.set_num_threads(data["nthreads"], quiet=True)

        mol = psi4.core.Molecule.from_schema(input_data)
        if mol.multiplicity() != 1:
            input_data["keywords"]["reference"] = "uks"

        ret_data = psi4.json_wrapper.run_json(input_data)

        # Handle slight RC2 weirdness once more
        if psi_12rc2_tweak:
            ret_data["schema_name"] = "qc_schema_output"
            ret_data["schema_version"] = 1

    else:
        raise TypeError("Psi4 version '{}' not understood".format(psi_version))

    # Dispatch errors, PSIO Errors are not recoverable for future runs
    if ret_data["success"] is False:

        if "PSIO Error" in ret_data["error"]:
            raise ValueError(ret_data["error"])

    # Move several pieces up a level
    if ret_data["success"]:
        ret_data["provenance"]["memory"] = input_data["memory"]
        ret_data["provenance"]["nthreads"] = input_data["nthreads"]
        del ret_data["memory"], input_data["nthreads"]

    return ret_data
