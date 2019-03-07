"""
Several import utilities
"""

import importlib
import io
import json
import operator
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import traceback
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union

from qcelemental.models import ComputeError, FailedOperation

from .config import LOGGER, get_provenance_augments

__all__ = ["compute_wrapper", "get_module_function", "model_wrapper", "handle_output_metadata"]


def model_wrapper(input_data: Dict[str, Any], model: 'BaseModel') -> 'BaseModel':
    """
    Wrap input data in the given model, or return a controlled error
    """

    try:
        if isinstance(input_data, dict):
            input_data = model(**input_data)
        elif isinstance(input_data, model):
            input_data = input_data.copy()
        else:
            raise KeyError("Input type of {} not understood.".format(type(model)))

        # Older QCElemental compat
        try:
            input_data.extras
        except AttributeError:
            input_data = input_data.copy(update={"extras": {}})

    except Exception:
        input_data = FailedOperation(
            input_data=input_data,
            success=False,
            error=ComputeError(
                error_type="input_error",
                error_message=("Input data could not be processed correctly:\n" + traceback.format_exc())))

    return input_data


@contextmanager
def compute_wrapper(capture_output: bool = True) -> Dict[str, Any]:
    """Wraps compute for timing, output capturing, and raise protection
    """

    metadata = {"stdout": None, "stderr": None}

    # Start timer
    comp_time = time.time()

    # Capture stdout/err
    new_stdout = io.StringIO()
    new_stderr = io.StringIO()
    if capture_output:

        old_stdout, sys.stdout = sys.stdout, new_stdout
        old_stderr, sys.stderr = sys.stderr, new_stderr

    try:
        yield metadata
        metadata["success"] = True
    except Exception as e:
        metadata["error_message"] = "QCEngine Call Error:\n" + traceback.format_exc()
        metadata["success"] = False

    # Place data
    metadata["wall_time"] = time.time() - comp_time
    if capture_output:
        sys.stdout = old_stdout
        sys.stderr = old_stderr
        # Pull over values
        metadata["stdout"] = new_stdout.getvalue() or None
        metadata["stderr"] = new_stderr.getvalue() or None


def get_module_function(module: str, func_name: str, subpackage=None) -> Callable[..., Any]:
    """Obtains a function from a given string

    Parameters
    ----------
    module : str
        The module to pull the function from
    func_name : str
        The name of the function to acquire, can be in a subpackage
    subpackage : None, optional
        Explicitly import a subpackage if required

    Returns
    -------
    ret : function
        The requested functions

    Example
    -------

    # Import numpy.linalg.eigh
    f = get_module_function("numpy", "linalg.eigh")
    f(np.ones((2, 2)))

    """
    # Will throw import error if we fail
    pkg = importlib.import_module(module, subpackage)

    return operator.attrgetter(func_name)(pkg)


def handle_output_metadata(output_data: Union[Dict[str, Any], 'BaseModel'],
                           metadata: Dict[str, Any],
                           raise_error: bool = False,
                           return_dict: bool = True) -> Union[Dict[str, Any], 'BaseModel']:
    """
    Fuses general metadata and output together.

    Returns
    -------
    result : dict or pydantic.models.Result
        Output type depends on return_dict or a dict if an error was generated in model construction
    """

    if isinstance(output_data, dict):
        output_fusion = output_data  # Error handling
    else:
        output_fusion = output_data.dict()

    # Do not override if computer generates
    output_fusion["stdout"] = output_fusion.get("stdout", None) or metadata["stdout"]
    output_fusion["stderr"] = output_fusion.get("stderr", None) or metadata["stderr"]

    if metadata["success"] is not True:
        output_fusion["success"] = False
        output_fusion["error"] = {"error_type": "meta_error", "error_message": metadata["error_message"]}

    # Raise an error if one exists and a user requested a raise
    if raise_error and (output_fusion["success"] is not True):
        msg = "stdout:\n{}".format(output_fusion["stdout"])
        msg += "\nstderr:\n{}".format(output_fusion["stderr"])
        LOGGER.info(msg)
        raise ValueError(output_fusion["error"]["error_message"])

    # Fill out provenance datadata
    provenance_augments = get_provenance_augments()
    provenance_augments["wall_time"] = metadata["wall_time"]
    if "provenance" in output_fusion:
        output_fusion["provenance"].update(provenance_augments)
    else:
        # Add onto the augments with some missing info
        provenance_augments["creator"] = "QCEngine"
        provenance_augments["version"] = provenance_augments["qcengine_version"]
        output_fusion["provenance"] = provenance_augments

    # Make sure pydantic sparsity is upheld
    for val in ["stdout", "stderr"]:
        if output_fusion[val] is None:
            output_fusion.pop(val)

    # We need to return the correct objects; e.g. Results, Procedures
    if output_fusion["success"]:
        # This will only execute if everything went well
        ret = output_data.__class__(**output_fusion)
    else:
        # Should only be reachable on failures
        ret = FailedOperation(
            success=output_fusion.pop("success", False), error=output_fusion.pop("error"), input_data=output_fusion)

    if return_dict:
        return json.loads(ret.json())  # Use Pydantic to serialize, then reconstruct as Python dict of Python Primals
    else:
        return ret


def terminate_process(proc: Any, timeout: int = 15) -> None:
    if proc.poll() is None:

        # Sigint (keyboard interupt)
        if sys.platform.startswith('win'):
            proc.send_signal(signal.CTRL_BREAK_EVENT)
        else:
            proc.send_signal(signal.SIGINT)

        try:
            start = time.time()
            while (proc.poll() is None) and (time.time() < (start + timeout)):
                time.sleep(0.02)

        # Flat kill
        finally:
            proc.kill()


@contextmanager
def popen(args: List[str], append_prefix: bool = False,
          popen_kwargs: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Opens a background task

    Code and idea from dask.distributed's testing suite
    https://github.com/dask/distributed
    """
    args = list(args)
    if popen_kwargs is None:
        popen_kwargs = {}
    else:
        popen_kwargs = popen_kwargs.copy()

    # Bin prefix
    if sys.platform.startswith('win'):
        bin_prefix = os.path.join(sys.prefix, 'Scripts')
    else:
        bin_prefix = os.path.join(sys.prefix, 'bin')

    # Do we prefix with Python?
    if append_prefix:
        args[0] = os.path.join(bin_prefix, args[0])

    if sys.platform.startswith('win'):
        # Allow using CTRL_C_EVENT / CTRL_BREAK_EVENT
        popen_kwargs['creationflags'] = subprocess.CREATE_NEW_PROCESS_GROUP

    popen_kwargs['stdout'] = subprocess.PIPE
    popen_kwargs['stderr'] = subprocess.PIPE
    LOGGER.info('Popen', args, popen_kwargs)
    ret = {"proc": subprocess.Popen(args, **popen_kwargs)}
    try:
        yield ret
    except Exception:
        raise

    finally:
        try:
            terminate_process(ret["proc"])
        finally:
            output, error = ret["proc"].communicate()
            ret["stdout"] = output.decode()
            ret["stderr"] = error.decode()


def execute(command: List[str],
            infiles: Optional[Dict[str, str]] = None,
            outfiles: Optional[List[str]] = None,
            *,
            scratch_name: Union[bool, str] = True,
            scratch_location: Optional[str] = None,
            scratch_messy: bool = False,
            blocking_files: Optional[List[str]] = None,
            timeout: Optional[int] = None,
            interupt_after: Optional[int] = None,
            environment: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """
    Runs a process in the background until complete.

    Returns True if exit code zero

    Parameters
    ----------
    command : list of str
    infiles : Dict[str] = str
        Input file names (names, not full paths) and contents.
        to be written in scratch dir. May be {}.
    outfiles : List[str] = None
        Output file names to be collected after execution into
        values. May be {}.
    scratch_name : str, optional
        Passed to scratch_directory
    scratch_location : str, optional
        Passed to scratch_directory
    scratch_messy : bool, optional
        Passed to scratch_directory
    blocking_files : list, optional
        Files which should stop execution if present beforehand.
    timeout : int, optional
        Stop the process after n seconds.
    interupt_after : int, optional
        Interupt the process (not hard kill) after n seconds.
    environment : dict, optional
        The environment to run in

    Raises
    ------
    FileExistsError
        If any file in `blocking` is present

    """

    # Format inputs
    if infiles is None:
        infiles = {}

    if outfiles is None:
        outfiles = []
    outfiles = {k: None for k in outfiles}

    # Check for blocking files
    if blocking_files is not None:
        for fl in blocking_files:
            if os.path.isfile(fl):
                raise FileExistsError('Existing file can interfere with execute operation.', fl)

    # Format popen
    popen_kwargs = {}
    if environment is not None:
        popen_kwargs["env"] = {k: v for k, v in environment.items() if v is not None}

    # Execute
    with scratch_directory(scratch_name, scratch_location, scratch_messy) as scrdir:
        popen_kwargs["cwd"] = scrdir
        with disk_files(infiles, outfiles, scrdir) as extrafiles:
            with popen(command, popen_kwargs=popen_kwargs) as proc:

                if interupt_after is None:
                    proc["proc"].wait(timeout=timeout)
                else:
                    time.sleep(interupt_after)
                    terminate_process(proc["proc"])

            retcode = proc["proc"].poll()
        proc['outfiles'] = extrafiles
    proc['scratch_directory'] = scrdir

    return retcode == 0, proc


@contextmanager
def scratch_directory(child: bool = True, parent: bool = None, messy: bool = False) -> str:
    """

    Parameters
    ----------
    child : str, optional
        By default, `True`, quarantine directory generated through
        `tempfile.mdktemp` so guaranteed unique and safe. When specified,
        quarantine directory has exactly `name`.
    parent : str, optional
        Create dirctory `child` elsewhere than TMP default.
    messy : bool, optional
        Leave scratch directory and contents on disk after completion.

    Yields
    ------
    str
        Full path of scratch directory.

    Raises
    ------
    FileExistsError
        If `child` specified and directory already exists (perhaps from a
        previous `messy=True` run).

    """
    if child is True:
        tmpdir = tempfile.mkdtemp(dir=parent)
    else:
        if parent is None:
            parent = Path(tempfile.gettempdir())
        else:
            parent = Path(parent)
        tmpdir = parent / child
        os.mkdir(tmpdir)
    try:
        yield tmpdir

    finally:
        if not messy:
            shutil.rmtree(tmpdir)
            LOGGER.info(f'... Removing {tmpdir}')


@contextmanager
def disk_files(infiles: Dict[str, str], outfiles: Dict[str, None], cwd: Optional[str] = None) -> Dict[str, str]:
    """

    Parameters
    ----------
    infiles : Dict[str] = str
        Input file names (names, not full paths) and contents.
        to be written in scratch dir. May be {}.
    outfiles : Dict[str] = None
        Output file names to be collected after execution into
        values. May be {}.
    cwd : str, optional
        Directory to which to write and read files.

    Yields
    ------
    Dict[str] = str
        outfiles with RHS filled in.

    """
    if cwd is None:
        lwd = Path.cwd()
    else:
        lwd = Path(cwd)
    try:
        # need an extra flag for text/binary
        for fl, content in infiles.items():
            filename = lwd / fl
            with open(filename, 'w') as fp:
                fp.write(content)
                LOGGER.info(f'... Writing: {filename}')

        yield outfiles

    finally:
        for fl in outfiles.keys():
            try:
                filename = lwd / fl
                with open(filename, 'r') as fp:
                    outfiles[fl] = fp.read()
                    LOGGER.info(f'... Writing: {filename}')
            except (OSError, FileNotFoundError) as err:
                outfiles[fl] = None


def which(command, return_bool=False):
    # environment is $PATH, less any None values
    lenv = {'PATH': ':' + os.environ.get('PATH')}
    lenv = {k: v for k, v in lenv.items() if v is not None}

    ans = shutil.which(command, mode=os.F_OK | os.X_OK, path=lenv['PATH'])

    if return_bool:
        return bool(ans)
    else:
        return ans
