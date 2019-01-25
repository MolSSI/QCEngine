"""
Integrates the computes together
"""

from qcelemental.models import ResultInput, ComputeError, OptimizationInput, Optimization, FailedOperation
from pydantic import ValidationError

from .programs import get_program
from .util import compute_wrapper, handle_output_metadata, get_module_function
from .config import get_config


def compute(input_data, program, raise_error=False, capture_output=True, local_options=None,
            return_dict=True):
    """Executes a single quantum chemistry program given a QC Schema input.

    The full specification can be found at:
        http://molssi-qc-schema.readthedocs.io/en/latest/index.html#

    Parameters
    ----------
    input_data : dict or qcelemental.models.ResultInput
        A QC Schema input specification in dictionary or model from QCElemental.models
    program : {"psi4", "rdkit"}
        The program to run the input under
    raise_error : bool, optional
        Determines if compute should raise an error or not.
    capture_output : bool, optional
        Determines if stdout/stderr should be captured.
    local_options : dict, optional
        A dictionary of local configuration options
    return_dict : bool, optional, default True
        Returns a dict instead of qcelemental.models.ResultInput

    Returns
    -------
    ret : dict, Result, FailedOperation
        A QC Schema output, type depends on return_dict key
        A FailedOperation returns

    """
    try:
        if isinstance(input_data, dict):
            input_data = ResultInput(**input_data)
    except ValidationError as val_error:
        failure = FailedOperation(input_data=input_data,
                                  success=False,
                                  error=ComputeError(error_type="input_validation_error",
                                                     error_message="Input could not be validated for the following "
                                                                   "reasons:\n{}".format(val_error.json())))
        if return_dict:
            return failure.dict()

        return failure

    except Exception as unk_err:
        failure = FailedOperation(input_data=input_data,
                                  success=False,
                                  error=ComputeError(
                                      error_type="unknown_input_error",
                                      error_message="Input could not be validated for unknown reasons, "
                                                    "likely due to invalid input data types:\n{}".format(unk_err)))
        if return_dict:
            return failure.dict()

        return failure

    if local_options is None:
        local_options = {}

    try:
        input_engine_options = input_data._qcengine_local_config
        input_data = input_data.copy(exclude={'_qcengine_local_config'})
    except AttributeError:
        input_engine_options = {}

    local_options = {**local_options, **input_engine_options}
    config = get_config(local_options=local_options)

    # Run the program
    with compute_wrapper(capture_output=capture_output) as metadata:
        try:
            output_data = get_program(program)(input_data, config)

        except KeyError as e:
            output_data = FailedOperation(input_data=input_data.dict(),
                                          success=False,
                                          error=ComputeError(
                                              error_type='program_error',
                                              error_message="QCEngine Call Error:\nProgram {} not understood."
                                                            "\nError Message: {}".format(program, str(e))))

        except ValidationError as e:
            output_data = FailedOperation(input_data=input_data.dict(),
                                          success=False,
                                          error=ComputeError(
                                             error_type="validation_error",
                                             error_message=e.json()))

        except Exception as e:
            output_data = FailedOperation(input_data=input_data.dict(),
                                          success=False,
                                          error=ComputeError(
                                              error_type="runtime_error",
                                              error_message="QCEngine Error:\nError Message: {}".format(str(e))
                                          ))

    return handle_output_metadata(output_data, metadata, raise_error=raise_error, return_dict=return_dict)


def compute_procedure(input_data, procedure, raise_error=False, capture_output=True, local_options=None,
                      return_dict=True):
    """Runs a procedure (a collection of the quantum chemistry executions)

    Parameters
    ----------
    input_data : dict or qcelemental.models.OptimizationInput
        A JSON input specific to the procedure executed in dictionary or model from QCElemental.models
    procedure : {"geometric"}
        The name of the procedure to run
    raise_error : bool, option
        Determines if compute should raise an error or not.
    capture_output : bool, optional
        Determines if stdout/stderr should be captured.
    local_options : dict, optional
        A dictionary of local configuration options
    return_dict : bool, optional, default True
        Returns a dict instead of qcelemental.models.ResultInput

    Returns
    ------
    dict, Optimization, FailedOperation
        A QC Schema representation of the requested output, type depends on return_dict key.
    """

    try:
        if isinstance(input_data, dict):
            input_data = OptimizationInput(**input_data)

    except ValidationError as val_error:
        # Fix this more procedure-centric
        failure = FailedOperation(input_data=input_data,
                                  success=False,
                                  error=ComputeError(
                                      error_type="input_validation_error",
                                      error_message="Input could not be validated for the following "
                                                    "reasons:\n{}".format(val_error.json())))
        if return_dict:
            return failure.dict()
        return failure

    except Exception as unk_err:
        failure = FailedOperation(input_data=input_data,
                                  success=False,
                                  error=ComputeError(
                                      error_type="unknown_input_error",
                                      error_message="Input could not be validated for unknown reasons, "
                                                    "likely due to invalid input data types:\n{}".format(unk_err)))
        if return_dict:
            return failure.dict()

        return failure

    config = get_config(local_options=local_options)

    # Run the procedure
    with compute_wrapper(capture_output=capture_output) as metadata:
        if procedure == "geometric":
            try:
                # Augment the input
                geometric_input = input_data.dict()
                geometric_input["input_specification"]["_qcengine_local_config"] = config.dict()
                output_data = get_module_function(procedure, "run_json.geometric_run_json")(geometric_input)
                if output_data["schema_name"] == "qc_schema_optimization_output":
                    output_data["schema_name"] = "qcschema_optimization_output"
                output_data["input_specification"].pop("_qcengine_local_config", None)
                output_data = Optimization(**output_data)

            except ModuleNotFoundError:
                output_data = FailedOperation(input_data=input_data.dict(),
                                                 success=False,
                                                 error=ComputeError(
                                                     error_type="import_error",
                                                     error_message="Could not import {}".format(procedure)
                                                 ))

            except ValidationError as e:
                output_data = FailedOperation(input_data=input_data.dict(),
                                              success=False,
                                              error=ComputeError(
                                                  error_type="validation_error",
                                                  error_message="Could not form output into Optimization, could be "
                                                                "bad program outputs, or some other error\n"
                                                                "Error Message: {}".format(e.json())))

            except Exception as e:
                output_data = FailedOperation(input_data=input_data.dict(),
                                              success=False,
                                              error=ComputeError(
                                                  error_type="runtime_error",
                                                  error_message="QCEngine RuntimeError:\n"
                                                                "Something went wrong in procedure execution, see "
                                                                "message for details:\n{}".format(str(e))))

        else:
            output_data = FailedOperation(input_data=input_data.dict(),
                                          success=False,
                                          error=ComputeError(
                                              error_type="program_error",
                                              error_message="QCEngine Call Error:"
                                                            "\nProcedure {} not understood".format(procedure))
)

    return handle_output_metadata(output_data, metadata, raise_error=raise_error, return_dict=return_dict)
