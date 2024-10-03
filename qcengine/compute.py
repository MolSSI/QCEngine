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
    input_data: Union[Dict[str, Any], "BaseModel"],  # TODO Input base class
    program: str,
    raise_error: bool = False,
    task_config: Optional[Dict[str, Any]] = None,
    local_options: Optional[Dict[str, Any]] = None,
    return_dict: bool = False,
    schema_version: int = -1,
) -> Union[
    "BaseModel", "FailedOperation", Dict[str, Any]
]:  # TODO Output base class, was AtomicResult OptimizationResult
    """Executes a single CMS program given a QCSchema input.

    The full specification can be found at:
        http://molssi-qc-schema.readthedocs.io/en/latest/index.html#

    Parameters
    ----------
    input_data
        A QCSchema input specification in dictionary or model from QCElemental.models
    program
        The CMS program or procedure with which to execute the input.
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
    schema_version
        The schema version to return. If -1, the input schema_version is used.

    Returns
    -------
    result
        AtomicResult, OptimizationResult, FailedOperation, etc., or Dict representation of any object type
        A QCSchema representation of the requested output, type depends on return_dict key.
    """

    try:
        output_data = input_data.model_copy()
    except AttributeError:
        output_data = input_data.copy()  # lgtm [py/multiple-definition]

    with compute_wrapper(capture_output=False, raise_error=raise_error) as metadata:

        # Grab the executor and build the input model
        try:
            executor = get_procedure(program)
        except InputError:
            executor = get_program(program)

        # Build the model and validate
        input_data = executor.build_input_model(input_data)  # calls model_wrapper

        # Build out task_config
        if task_config is None:
            task_config = {}
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


def compute_procedure(*args, **kwargs):
    warnings.warn(
        "Using the `compute_procedure` function is deprecated in favor of using `compute`, "
        "and as soon as version 0.70.0 it will stop working.",
        category=FutureWarning,
        stacklevel=2,
    )
    if "procedure" in kwargs:
        kwargs["program"] = kwargs.pop("procedure")
    return compute(*args, **kwargs)
