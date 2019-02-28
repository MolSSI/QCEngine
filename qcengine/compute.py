"""
Integrates the computes together
"""
from typing import Any, Dict, Optional, Union

from qcelemental.models import ComputeError, FailedOperation, Optimization, OptimizationInput, ResultInput

from .config import get_config
from .programs import get_program
from .util import compute_wrapper, get_module_function, handle_output_metadata, model_wrapper

__all__ = ["compute", "compute_procedure"]


def compute(input_data: Union[Dict[str, Any], 'ResultInput'],
            program: str,
            raise_error: bool = False,
            local_options: Optional[Dict[str, str]] = None,
            return_dict: bool = False) -> 'Result':
    """Executes a single quantum chemistry program given a QC Schema input.

    The full specification can be found at:
        http://molssi-qc-schema.readthedocs.io/en/latest/index.html#

    Parameters
    ----------
    input_data:
        A QC Schema input specification in dictionary or model from QCElemental.models
    program:
        The program to execute the input with
    raise_error:
        Determines if compute should raise an error or not.
    capture_output:
        Determines if stdout/stderr should be captured.
    local_options:
        A dictionary of local configuration options
    return_dict:
        Returns a dict instead of qcelemental.models.ResultInput

    Returns
    -------
    : Result
        A QC Schema output or type depending on return_dict key

    """

    # Build the model and validate
    input_data = model_wrapper(input_data, ResultInput, raise_error)
    if isinstance(input_data, FailedOperation):
        if return_dict:
            return input_data.dict()
        return input_data

    # Build out local options
    if local_options is None:
        local_options = {}

    input_engine_options = input_data.extras.pop("_qcengine_local_config", {})

    local_options = {**local_options, **input_engine_options}
    config = get_config(local_options=local_options)

    # Run the program
    with compute_wrapper(capture_output=False) as metadata:

        output_data = input_data.copy()  # Initial in case of error handling
        try:
            output_data = get_program(program).compute(input_data, config)
        except KeyError as e:
            output_data = FailedOperation(
                input_data=output_data.dict(),
                success=False,
                error=ComputeError(
                    error_type='program_error',
                    error_message="QCEngine Call Error:\nProgram {} not understood."
                    "\nError Message: {}".format(program, str(e))))

    return handle_output_metadata(output_data, metadata, raise_error=raise_error, return_dict=return_dict)


def compute_procedure(input_data: Dict[str, Any],
                      procedure: str,
                      raise_error: bool = False,
                      local_options: Optional[Dict[str, str]] = None,
                      return_dict: bool = False) -> 'BaseModel':
    """Runs a procedure (a collection of the quantum chemistry executions)

    Parameters
    ----------
    input_data : dict or qcelemental.models.OptimizationInput
        A JSON input specific to the procedure executed in dictionary or model from QCElemental.models
    procedure : {"geometric"}
        The name of the procedure to run
    raise_error : bool, option
        Determines if compute should raise an error or not.
    local_options : dict, optional
        A dictionary of local configuration options
    return_dict : bool, optional, default True
        Returns a dict instead of qcelemental.models.ResultInput

    Returns
    ------
    dict, Optimization, FailedOperation
        A QC Schema representation of the requested output, type depends on return_dict key.
    """

    # Build the model and validate
    input_data = model_wrapper(input_data, OptimizationInput, raise_error)
    if isinstance(input_data, FailedOperation):
        if return_dict:
            return input_data.dict()
        return input_data

    config = get_config(local_options=local_options)

    # Run the procedure
    with compute_wrapper(capture_output=False) as metadata:
        # Create a base output data in case of errors
        output_data = input_data.copy()  # lgtm [py/multiple-definition]
        if procedure == "geometric":
            # Augment the input
            geometric_input = input_data.dict()

            # Older QCElemental compat, can be removed in v0.6
            if "extras" not in geometric_input["input_specification"]:
                geometric_input["input_specification"]["extras"] = {}

            geometric_input["input_specification"]["extras"]["_qcengine_local_config"] = config.dict()

            # Run the program
            output_data = get_module_function(procedure, "run_json.geometric_run_json")(geometric_input)

            output_data["schema_name"] = "qcschema_optimization_output"
            output_data["input_specification"]["extras"].pop("_qcengine_local_config", None)
            if output_data["success"]:
                output_data = Optimization(**output_data)

        else:
            output_data = FailedOperation(
                input_data=input_data.dict(),
                success=False,
                error=ComputeError(
                    error_type="program_error",
                    error_message="QCEngine Call Error:"
                    "\nProcedure {} not understood".format(procedure)))

    return handle_output_metadata(output_data, metadata, raise_error=raise_error, return_dict=return_dict)
