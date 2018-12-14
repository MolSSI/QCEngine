"""
Calls the Psi4 executable.
"""

import sys

from pkg_resources import parse_version

from qcengine import config


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
    input_data["memory"] = int(config.get_config("memory_per_job") * 1024 * 1024 * 1024 * 0.95)
    input_data["success"] = False

    scratch = config.get_config("scratch_directory")
    if scratch is not None:
        input_data["scratch_location"] = scratch

    psi_version = _parse_psi_version(psi4.__version__)

    if psi_version > parse_version("1.2rc2.dev500"):

        # Handle slight RC2 weirdness
        psi_12rc2_tweak = (psi_version == parse_version("1.2rc2"))
        psi_12rc2_tweak &= psi4.metadata.version_formatter("") != '(inplace)'

        if psi_12rc2_tweak:
            input_data["schema_name"] = "QC_JSON"
            input_data["schema_version"] = 0
            psi4.set_num_threads(input_data["nthreads"], quiet=True)

        mol = psi4.core.Molecule.from_schema(input_data)
        if mol.multiplicity() != 1:
            input_data["keywords"]["reference"] = "uks"

        output_data = psi4.json_wrapper.run_json(input_data)

        # Handle slight RC2 weirdness once more
        if psi_12rc2_tweak:
            output_data["schema_name"] = "qc_schema_output"
            output_data["schema_version"] = 1

    else:
        raise TypeError("Psi4 version '{}' not understood".format(psi_version))

    # Dispatch errors, PSIO Errors are not recoverable for future runs
    if output_data["success"] is False:

        if "PSIO Error" in output_data["error"]:
            raise ValueError(output_data["error"])

    # Move several pieces up a level
    if output_data["success"]:
        output_data["provenance"]["memory"] = round(input_data["memory"] / (1024**3), 3)
        output_data["provenance"]["nthreads"] = input_data["nthreads"]
        del output_data["memory"], input_data["nthreads"]

    return output_data
