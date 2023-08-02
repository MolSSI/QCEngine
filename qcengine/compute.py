"""
Integrates the computes together
"""
import warnings
from typing import TYPE_CHECKING, Any, Dict, Optional, Union

from qcelemental.models import AtomicInput, AtomicResult, FailedOperation, OptimizationResult

from .config import get_config
from .exceptions import InputError, RandomError
from .procedures import get_procedure
from .programs import get_program
from .util import compute_wrapper, environ_context, handle_output_metadata, model_wrapper

if TYPE_CHECKING:
    try:
        from pydantic.v1.main import BaseModel
    except ImportError:
        from pydantic.main import BaseModel
    from qcelemental.models import AtomicResult


__all__ = ["compute", "compute_procedure"]


def _process_failure_and_return(model, return_dict, raise_error):
    if isinstance(model, FailedOperation):
        if raise_error:
            raise InputError(model.error.error_message)
        elif return_dict:
            return model.dict()
        else:
            return model
    else:
        return False


def compute(
    input_data: Union[Dict[str, Any], "AtomicInput"],
    program: str,
    raise_error: bool = False,
    task_config: Optional[Dict[str, Any]] = None,
    local_options: Optional[Dict[str, Any]] = None,
    return_dict: bool = False,
) -> Union["AtomicResult", "FailedOperation", Dict[str, Any]]:
    """Executes a single CMS program given a QCSchema input.

    The full specification can be found at:
        http://molssi-qc-schema.readthedocs.io/en/latest/index.html#

    Parameters
    ----------
    input_data
        A QCSchema input specification in dictionary or model from QCElemental.models
    program
        The CMS program with which to execute the input.
    raise_error
        Determines if compute should raise an error or not.
    retries : int, optional
        The number of random tries to retry for.
    task_config
        A dictionary of local configuration options corresponding to a TaskConfig object.
    local_options
        Deprecated parameter, renamed to ``task_config``
    return_dict
        Returns a dict instead of qcelemental.models.AtomicResult

    Returns
    -------
    result
        AtomicResult, FailedOperation, or Dict representation of either object type
        A QCSchema representation of the requested output, type depends on return_dict key.
    """

    output_data = input_data.copy()  # lgtm [py/multiple-definition]
    with compute_wrapper(capture_output=False, raise_error=raise_error) as metadata:

        # Grab the executor and build the input model
        executor = get_program(program)

        # Build the model and validate
        input_data = model_wrapper(input_data, AtomicInput)

        # Build out task_config
        if task_config is None:
            task_config = {}

        if local_options:
            warnings.warn(
                "Using the `local_options` keyword argument is deprecated in favor of using `task_config`, "
                "in version 0.30.0 it will stop working.",
                category=FutureWarning,
                stacklevel=2,
            )
            task_config = {**local_options, **task_config}

        input_engine_options = input_data.extras.pop("_qcengine_local_config", {})

        task_config = {**task_config, **input_engine_options}
        config = get_config(task_config=task_config)

        # Set environment parameters and execute
        with environ_context(config=config):

            # Handle optional retries
            for x in range(config.retries + 1):
                try:
                    output_data = executor.compute(input_data, config)
                    break
                except RandomError as e:

                    if x == config.retries:
                        raise e
                    else:
                        metadata["retries"] += 1
                except:
                    raise

    return handle_output_metadata(output_data, metadata, raise_error=raise_error, return_dict=return_dict)


def compute_procedure(
    input_data: Union[Dict[str, Any], "BaseModel"],
    procedure: str,
    raise_error: bool = False,
    task_config: Optional[Dict[str, str]] = None,
    local_options: Optional[Dict[str, str]] = None,
    return_dict: bool = False,
) -> Union["OptimizationResult", "FailedOperation", Dict[str, Any]]:
    """Runs a procedure (a collection of the quantum chemistry executions)

    Parameters
    ----------
    input_data : dict or qcelemental.models.OptimizationInput
        A JSON input specific to the procedure executed in dictionary or model from QCElemental.models
    procedure : {"geometric", "berny"}
        The name of the procedure to run
    raise_error : bool, option
        Determines if compute should raise an error or not.
    task_config
        A dictionary of local configuration options corresponding to a TaskConfig object.
    local_options
        Deprecated parameter, renamed to ``task_config``
    return_dict : bool, optional, default True
        Returns a dict instead of qcelemental.models.AtomicInput

    Returns
    ------
    dict, OptimizationResult, FailedOperation
        A QC Schema representation of the requested output, type depends on return_dict key.
    """
    # Build out task_config
    if task_config is None:
        task_config = {}

    if local_options:
        warnings.warn(
            "Using the `local_options` keyword argument is depreciated in favor of using `task_config`, "
            "in version 0.30.0 it will stop working.",
            category=FutureWarning,
            stacklevel=2,
        )
        task_config = {**local_options, **task_config}

    output_data = input_data.copy()  # lgtm [py/multiple-definition]
    with compute_wrapper(capture_output=False, raise_error=raise_error) as metadata:

        # Grab the executor and build the input model
        executor = get_procedure(procedure)

        config = get_config(task_config=task_config)
        input_data = executor.build_input_model(input_data)

        # Create a base output data in case of errors
        output_data = input_data.copy()  # lgtm [py/multiple-definition]

        # Set environment parameters and execute
        with environ_context(config=config):
            output_data = executor.compute(input_data, config)

    return handle_output_metadata(output_data, metadata, raise_error=raise_error, return_dict=return_dict)
