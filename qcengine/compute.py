"""
Integrates the computes together
"""
from typing import Any, Dict, Optional, Union

from qcelemental.models import ComputeError, FailedOperation, Optimization, OptimizationInput, ResultInput

from .config import get_config
from .procedures import get_procedure, list_all_procedures, list_available_procedures
from .programs import get_program, list_all_programs, list_available_programs
from .util import compute_wrapper, get_module_function, handle_output_metadata, model_wrapper

__all__ = ["compute", "compute_procedure"]


def _process_failure_and_return(model, return_dict, raise_error):
    if isinstance(model, FailedOperation):
        if raise_error:
            raise ValueError(model.error.error_message)
        elif return_dict:
            return model.dict()
        else:
            return model
    else:
        return False


def compute(input_data: Union[Dict[str, Any], 'ResultInput'],
            program: str,
            raise_error: bool=False,
            local_options: Optional[Dict[str, str]]=None,
            return_dict: bool=False) -> 'Result':
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

    program = program.lower()
    if program not in list_all_programs():
        input_data = FailedOperation(
            input_data=input_data,
            error=ComputeError(
                error_type="not_registered",
                error_message="QCEngine Call Error:\n"
                "Program {} is not registered with QCEngine".format(program)))
    elif program not in list_available_programs():
        input_data = FailedOperation(
            input_data=input_data,
            error=ComputeError(
                error_type="not_available",
                error_message="QCEngine Call Error:\n"
                "Program {} is registered with QCEngine, but cannot be found".format(program)))
    error = _process_failure_and_return(input_data, return_dict, raise_error)
    if error:
        return error

    # Build the model and validate
    input_data = model_wrapper(input_data, ResultInput)
    error = _process_failure_and_return(input_data, return_dict, raise_error)
    if error:
        return error

    # Grab the executor and build the input model
    executor = get_program(program)

    # Build out local options
    if local_options is None:
        local_options = {}

    input_engine_options = input_data.extras.pop("_qcengine_local_config", {})

    local_options = {**local_options, **input_engine_options}
    config = get_config(local_options=local_options)

    # Run the program
    with compute_wrapper(capture_output=False) as metadata:

        output_data = input_data.copy()  # lgtm [py/multiple-definition]
        output_data = executor.compute(input_data, config)

    return handle_output_metadata(output_data, metadata, raise_error=raise_error, return_dict=return_dict)


def compute_procedure(input_data: Union[Dict[str, Any], 'BaseModel'],
                      procedure: str,
                      raise_error: bool=False,
                      local_options: Optional[Dict[str, str]]=None,
                      return_dict: bool=False) -> 'BaseModel':
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

    procedure = procedure.lower()
    if procedure not in list_all_procedures():
        input_data = FailedOperation(
            input_data=input_data,
            error=ComputeError(
                error_type="not_registered",
                error_message="QCEngine Call Error:\n"
                "Procedure {} is not registered with QCEngine".format(procedure)))
    elif procedure not in list_available_procedures():
        input_data = FailedOperation(
            input_data=input_data,
            error=ComputeError(
                error_type="not_available",
                error_message="QCEngine Call Error:\n"
                "Procedure {} is registered with QCEngine, but cannot be found".format(procedure)))
    error = _process_failure_and_return(input_data, return_dict, raise_error)
    if error:
        return error

    # Grab the executor and build the input model
    executor = get_procedure(procedure)

    config = get_config(local_options=local_options)
    input_data = executor.build_input_model(input_data)
    error = _process_failure_and_return(input_data, return_dict, raise_error)
    if error:
        return error

    # Run the procedure
    with compute_wrapper(capture_output=False) as metadata:

        # Create a base output data in case of errors
        output_data = input_data.copy()  # lgtm [py/multiple-definition]
        output_data = executor.compute(input_data, config)

    return handle_output_metadata(output_data, metadata, raise_error=raise_error, return_dict=return_dict)
