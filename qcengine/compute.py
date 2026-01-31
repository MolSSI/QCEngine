"""
Integrates the computes together
"""
import os
import sys
import warnings
from typing import TYPE_CHECKING, Any, Dict, Optional, Union

import qcelemental

from .config import get_config
from .exceptions import InputError, RandomError
from .procedures import get_procedure
from .programs import get_program
from .util import QCEL_V1V2_SHIM_CODE, compute_wrapper, environ_context, handle_output_metadata, model_wrapper

if TYPE_CHECKING:
    from pydantic.main import BaseModel


__all__ = ["compute", "compute_procedure"]


def _process_failure_and_return(model, return_dict, raise_error):
    if isinstance(model, (qcelemental.models.v1.FailedOperation, qcelemental.models.v2.FailedOperation)):
        if raise_error:
            raise InputError(model.error.error_message)
        elif return_dict:
            return model.model_dump()
        else:
            return model
    else:
        return False


def compute(
    input_data: Union[Dict[str, Any], "BaseModel"],  # TODO Input base class
    program: str,
    raise_error: bool = False,
    task_config: Optional[Dict[str, Any]] = None,
    return_dict: bool = False,
    return_version: int = -1,
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
        The CMS program or procedure with which to execute the input. E.g., "psi4", "rdkit", "geometric".
    raise_error
        If computation doesn't succeed, should compute raise an error (`True`) or encode in a
        `FailedOperation` (`False`). See also: `return_dict`.
    retries : int, optional
        The number of random tries to retry for.
    task_config
        A dictionary of local configuration options corresponding to a TaskConfig object.
        Formerly local_options.
    return_dict
        Returns a dictionary serialization instead of a Result or FailedOperation model from
        QCElemental.models. Note that `raise_error` take precedence.
    return_version
        The schema version to return. If -1, the input schema_version is used.

    Returns
    -------
    result
        AtomicResult, OptimizationResult, FailedOperation, etc., or Dict representation of any object type
        A QCSchema representation of the requested output, type depends on return_dict key.

        .. _`table:compute_result`:

        +-----------+-------------+-------------+--------------------------------------------------+
        | good_calc | raise_error | return_dict | output                                           |
        +===========+=============+=============+==================================================+
        | T         | F (def)     | F (def)     | ``AtomicResult`` object                          |
        +-----------+-------------+-------------+--------------------------------------------------+
        | T         | T           | F (def)     | ``AtomicResult`` object                          |
        +-----------+-------------+-------------+--------------------------------------------------+
        | T         | T           | T           | dict of ``AtomicResult``                         |
        +-----------+-------------+-------------+--------------------------------------------------+
        | T         | F (def)     | T           | dict of ``AtomicResult``                         |
        +-----------+-------------+-------------+--------------------------------------------------+
        | F         | F (def)     | F (def)     | ``FailedOperation`` object                       |
        +-----------+-------------+-------------+--------------------------------------------------+
        | F         | T           | F (def)     | raises ``InputError`` (type encoded in FailedOp) |
        +-----------+-------------+-------------+--------------------------------------------------+
        | F         | T           | T           | raises ``InputError`` (type encoded in FailedOp) |
        +-----------+-------------+-------------+--------------------------------------------------+
        | F         | F (def)     | T           | dict of ``FailedOperation``                      |
        +-----------+-------------+-------------+--------------------------------------------------+

        .. _`table:compute_result_schver`:

        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | good_calc? | input_data ver | return_version | return_dict | output                                                                 |
        +============+================+================+=============+========================================================================+
        | T          | 1              | -1 (def) or 1  | F (def)     | ``v1.AtomicResult`` object (not avail. Py 3.14+ ^^)                    |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | T          | 1              |  2             | F (def)     | ``v2.AtomicResult`` object                                             |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | T          | 1              | -1 (def) or 1  | T           | dict of ``v1.AtomicResult`` @@                                         |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | T          | 1              |  2             | T           | dict of ``v2.AtomicResult``                                            |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | T          | 2              | -1 (def) or 2  | F (def)     | ``v2.AtomicResult`` object                                             |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | T          | 2              |  1             | F (def)     | ``v1.AtomicResult`` object (not avail. Py 3.14+ ^^)                    |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | T          | 2              | -1 (def) or 2  | T           | dict of ``v2.AtomicResult``                                            |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | T          | 2              |  1             | T           | dict of ``v1.AtomicResult`` @@                                         |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | F          | 1              | -1 (def) or 1  | F (def)     | ``v1.FailedOperation`` object (not avail. Py 3.14+ ^^) **              |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | F          | 1              |  2             | F (def)     | ``v2.FailedOperation`` object                                          |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | F          | 1              | -1 (def) or 1  | T           | dict of ``v1.FailedOperation`` @@ **                                   |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | F          | 1              |  2             | T           | dict of ``v2.FailedOperation``                                         |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | F          | 2              | -1 (def) or 2  | F (def)     | ``v2.FailedOperation`` object **                                       |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | F          | 2              |  1             | F (def)     | ``v1.FailedOperation`` object (not avail. Py 3.14+ ^^)                 |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | F          | 2              | -1 (def) or 2  | T           | dict of ``v2.FailedOperation`` **                                      |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | F          | 2              |  1             | T           | dict of ``v1.FailedOperation`` @@                                      |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        | ** F       |                | ** -1          |             | for errors *before* input ver detected, returns v2 if Py 3.14+ else v1 |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        |            |                | ^^ 1           | ^^ F        | for Py 3.14+, returns v2.FailedOp for minimal API surprise             |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+
        |            |                | @@ 1           | @@ T        | for Py 3.14+, returns dict of v1 (like all other Py) with brittle      |
        +            +                +                +             +                                                                        +
        |            |                |                |             |   :envvar:`QCNG_USE_V1V2_SHIM=1`, otherwise dict of v2.FailedOp        |
        +------------+----------------+----------------+-------------+------------------------------------------------------------------------+

    .. versionadded:: 0.50.0
       input_data can newly be QCSchema v2 as well as longstanding v1.
       Note that the QCSchema v2 layout will not be official until v0.60.0 .
       The compute_procedure is newly incorporated into compute.
       The *return_version* parameter was added.

    """

    try:
        # models, v1 or v2
        output_data = input_data.model_copy()
    except AttributeError:
        # dicts
        output_data = input_data.copy()  # lgtm [py/multiple-definition]

    with compute_wrapper(capture_output=False, raise_error=raise_error) as metadata:
        # Grab the executor harness
        try:
            executor = get_procedure(program)
        except InputError:
            executor = get_program(program)

        # Build the model and validate
        # * calls model_wrapper with the (Atomic|Optimization|etc)Input for which the harness was designed
        # * upon return, input_data is a model of the type (e.g., Atomic) and version (e.g., 1 or 2) the harness prefers: all v2.
        input_data, input_schema_version = executor.build_input_model(input_data, return_input_schema_version=True)
        return_version = input_schema_version if return_version == -1 else return_version

        # V1V2TEST if return_version == 1:
        if sys.version_info >= (3, 14) and return_version == 1:  # forbidden by pydantic ...

            _MSG314 = (
                f"QCSchema v1 models (this: {input_data.schema_name}) cannot be instantiated in this environment. "
                + "Reason: pydantic.v1 is unavailable on Python 3.14+. You can: "
                + "(a) use Python <3.14; or "
                + "(b) use QCSchema v2 (either input or ask for output with `return_version=2`; "
                + " note that QCSchema v2 not finalized until QCElemental v0.60); or "
                + "(c) ask for a QCSchema v1 dictionary back rather than the model with `return_dict=True` "
                + " while setting envvar `QCNG_USE_V1V2_SHIM=1` to acknowledge this is brittle."
            )

            # V1V2TEST if input_data.schema_name == "qcschema_atomic_input":
            if (
                return_dict is True and input_data.schema_name == "qcschema_atomic_input"
            ):  # ... but we have a workaround ...

                if bool(os.environ.get("QCNG_USE_V1V2_SHIM", False)):  # ... if choose to use it
                    # return_version = -12 = QCEL_V1V2_SHIM_CODE signals to use the shim classes that represent certain
                    #   QCSchema v1 layouts (QCSk v1 only exist in pydantic.v1 API) in pydantic v2 API.
                    #   We never want to release these into the wild so only available if returning
                    #   dict and only for Atomic models (i.e., shims not avail for Opt, TD, MBE).
                    return_version = QCEL_V1V2_SHIM_CODE
                else:
                    # reset retver so that FailedOp w/_MSG314 can be constructed
                    return_version = 2
                    raise RuntimeError(_MSG314)
            else:
                # reset retver so that FailedOp w/_MSG314 can be constructed
                return_version = 2
                raise RuntimeError(_MSG314)

        # Build out task_config
        if task_config is None:
            task_config = {}
        input_engine_options = input_data.specification.extras.pop("_qcengine_local_config", {})
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
                    if return_version >= 2:
                        output_data = input_data

                    if x == config.retries:
                        raise e
                    else:
                        metadata["retries"] += 1
                except:
                    if return_version >= 2:
                        output_data = input_data
                    raise

    return handle_output_metadata(
        output_data, metadata, raise_error=raise_error, return_dict=return_dict, convert_version=return_version
    )


def compute_procedure(*args, **kwargs):
    from qcelemental.models.common_models import _qcsk_v2_default_v1_importpathschange

    warnings.warn(
        "Using the `compute_procedure` function is deprecated in favor of using `compute`, "
        f"and as soon as version {_qcsk_v2_default_v1_importpathschange} it will stop working.",
        category=FutureWarning,
        stacklevel=2,
    )
    if "procedure" in kwargs:
        kwargs["program"] = kwargs.pop("procedure")
    return compute(*args, **kwargs)
