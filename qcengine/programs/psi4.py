"""
Calls the Psi4 executable.
"""

from pkg_resources import parse_version
from qcelemental.models import Result, FailedOperation


def _parse_psi_version(version):
    if "undef" in version:
        raise TypeError(
            "Using custom build Psi4 without tags. Please `git pull origin master --tags` and recompile Psi4.")

    return parse_version(version)


def psi4(input_model, config):
    """
    Runs Psi4 in API mode
    """

    try:
        import psi4
    except ImportError:
        raise ImportError("Could not find Psi4 in the Python path.")

    # Setup the job
    input_data = input_model.copy().dict()
    input_data["nthreads"] = config.ncores
    input_data["memory"] = int(config.memory * 1024 * 1024 * 1024 * 0.95)  # Memory in bytes
    input_data["success"] = False
    reset_schema = False
    if input_data["schema_name"] == "qcschema_input":
        input_data["schema_name"] = "qc_schema_input"
        reset_schema = True

    scratch = config.scratch_directory
    if scratch is not None:
        input_data["scratch_location"] = scratch

    psi_version = _parse_psi_version(psi4.__version__)

    if psi_version > parse_version("1.2"):

        mol = psi4.core.Molecule.from_schema(input_data)
        if mol.multiplicity() != 1:
            input_data["keywords"]["reference"] = "uks"

        output_data = psi4.json_wrapper.run_json(input_data)

    else:
        raise TypeError("Psi4 version '{}' not understood.".format(psi_version))

    if reset_schema:
        output_data["schema_name"] = "qcschema_input"
    # Dispatch errors, PSIO Errors are not recoverable for future runs
    if output_data["success"] is False:

        if "PSIO Error" in output_data["error"]:
            raise ValueError(output_data["error"])

    # Move several pieces up a level
    if output_data["success"]:
        output_data["provenance"]["memory"] = round(input_data["memory"] / (1024**3), 3)
        output_data["provenance"]["nthreads"] = input_data["nthreads"]
        del output_data["memory"], input_data["nthreads"]
        return Result(**output_data)
    return FailedOperation(
        success=output_data.pop("success", False), error=output_data.pop("error"), input_data=output_data)
