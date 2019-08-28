"""
Utilities for the testing suite.
"""

import os
import signal
import sys
import time
import subprocess
from contextlib import contextmanager
from typing import List

import numpy as np
import pytest
from pkg_resources import parse_version

import qcelemental as qcel
import qcengine as qcng
from qcelemental.util import which, which_import


def is_program_new_enough(program, version_feature_introduced):
    """Returns True if `program` registered in QCEngine, locatable in
    environment, has parseable version, and that version in normalized
    form is equal to or later than `version_feature_introduced`.

    """
    if program not in qcng.list_available_programs():
        return False
    candidate_version = qcng.get_program(program).get_version()

    return parse_version(candidate_version) >= parse_version(version_feature_introduced)


@pytest.fixture(scope="function")
def failure_engine():
    unique_name = "testing_random_name"

    class FailEngine(qcng.programs.ProgramHarness):
        iter_modes: List[str] = []
        ncalls: int = 0
        start_distance: float = 5
        equilibrium_distance: float = 4

        _defaults = {
            "name": unique_name,
            "scratch": False,
            "thread_safe": True,
            "thread_parallel": False,
            "node_parallel": False,
            "managed_memory": False,
        }

        class Config(qcng.programs.ProgramHarness.Config):
            allow_mutation: True

        @staticmethod
        def found(raise_error: bool = False) -> bool:
            return True

        def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
            self.ncalls += 1
            mode = self.iter_modes.pop(0)

            geom = input_data.molecule.geometry
            if geom.shape[0] != 2:
                raise ValueError("Failure Test must have an input size of two.")

            grad_value = np.abs(np.linalg.norm(geom[0] - geom[1]) - self.equilibrium_distance)
            grad = [0, 0, -grad_value, 0, 0, grad_value]

            if mode == "pass":
                return qcel.models.Result(
                    **{
                        **input_data.dict(),
                        **{
                            "properties": {
                                "return_energy": grad_value
                            },
                            "return_result": grad,
                            "success": True,
                            "extras": {
                                "ncalls": self.ncalls
                            },
                            "provenance": {
                                "creator": "failure_engine",
                                "ncores": config.ncores
                            }
                        }
                    })
            elif mode == "random_error":
                raise qcng.exceptions.RandomError("Whoops!")
            elif mode == "input_error":
                raise qcng.exceptions.InputError("Whoops!")
            else:
                raise KeyError("Testing error, should not arrive here.")

        def get_job(self):
            json_data = {
                "molecule": {
                    "symbols": ["He", "He"],
                    "geometry": [0, 0, 0, 0, 0, self.start_distance]
                },
                "driver": "gradient",
                "model": {
                    "method": "something"
                }
            }

            return json_data

    engine = FailEngine()
    qcng.register_program(engine)

    yield engine

    qcng.unregister_program(engine.name)


# Figure out what is imported
_programs = {
    "dftd3": which('dftd3', return_bool=True),
    "geometric": which_import("geometric", return_bool=True),
    "psi4": is_program_new_enough("psi4", "1.2"),
    "rdkit": which_import("rdkit", return_bool=True),
    "qcdb": which_import("qcdb", return_bool=True),
    "torchani": which_import("torchani", return_bool=True),
    "mp2d": which('mp2d', return_bool=True),
    "terachem": which("terachem", return_bool=True),
    "molpro": is_program_new_enough("molpro", "2018.1"),
    "mopac": is_program_new_enough("mopac", "2016"),
    "entos": is_program_new_enough("entos", "0.5"),
    "cfour": which('xcfour', return_bool=True),
    "gamess": which('rungms', return_bool=True),
    "nwchem": which('nwchem', return_bool=True),
}


def has_program(name):
    return _programs[name]


def _build_pytest_skip(program):
    import_message = "Not detecting module {}. Install package if necessary to enable tests."
    return pytest.mark.skipif(has_program(program) is False, reason=import_message.format(program))


def terminate_process(proc):
    if proc.poll() is None:

        # Sigint (keyboard interrupt)
        if sys.platform.startswith('win'):
            proc.send_signal(signal.CTRL_BREAK_EVENT)
        else:
            proc.send_signal(signal.SIGINT)

        try:
            start = time.time()
            while (proc.poll() is None) and (time.time() < (start + 15)):
                time.sleep(0.02)
        # Flat kill
        finally:
            proc.kill()


@contextmanager
def popen(args, **kwargs):
    """
    Opens a background task.

    Code and idea from dask.distributed's testing suite
    https://github.com/dask/distributed
    """
    args = list(args)

    # Bin prefix
    if sys.platform.startswith('win'):
        bin_prefix = os.path.join(sys.prefix, 'Scripts')
    else:
        bin_prefix = os.path.join(sys.prefix, 'bin')

    # Do we prefix with Python?
    if kwargs.pop("append_prefix", True):
        args[0] = os.path.join(bin_prefix, args[0])

    # Add coverage testing
    if kwargs.pop("coverage", False):
        coverage_dir = os.path.join(bin_prefix, "coverage")
        if not os.path.exists(coverage_dir):
            print("Could not find Python coverage, skipping cov.")

        else:
            src_dir = os.path.dirname(os.path.abspath(__file__))
            coverage_flags = [coverage_dir, "run", "--append", "--source=" + src_dir]

            # If python script, skip the python bin
            if args[0].endswith("python"):
                args.pop(0)
            args = coverage_flags + args

    # Do we optionally dumpstdout?
    dump_stdout = kwargs.pop("dump_stdout", False)

    if sys.platform.startswith('win'):
        # Allow using CTRL_C_EVENT / CTRL_BREAK_EVENT
        kwargs['creationflags'] = subprocess.CREATE_NEW_PROCESS_GROUP

    kwargs['stdout'] = subprocess.PIPE
    kwargs['stderr'] = subprocess.PIPE
    proc = subprocess.Popen(args, **kwargs)
    try:
        yield proc
    except Exception:
        dump_stdout = True
        raise

    finally:
        try:
            terminate_process(proc)
        finally:
            output, error = proc.communicate()
            if dump_stdout:
                print('\n' + '-' * 30)
                print("\n|| Process command: {}".format(" ".join(args)))
                print('\n|| Process stderr: \n{}'.format(error.decode()))
                print('-' * 30)
                print('\n|| Process stdout: \n{}'.format(output.decode()))
                print('-' * 30)


def run_process(args, **kwargs):
    """
    Runs a process in the background until complete.

    Returns True if exit code zero.
    """

    timeout = kwargs.pop("timeout", 30)
    terminate_after = kwargs.pop("interupt_after", None)
    with popen(args, **kwargs) as proc:
        if terminate_after is None:
            proc.wait(timeout=timeout)
        else:
            time.sleep(terminate_after)
            terminate_process(proc)

        retcode = proc.poll()

    return retcode == 0


# Add flags
using_dftd3 = _build_pytest_skip("dftd3")
using_entos = _build_pytest_skip("entos")
using_geometric = _build_pytest_skip("geometric")
using_mopac = _build_pytest_skip("mopac")
using_molpro = _build_pytest_skip("molpro")
using_mp2d = _build_pytest_skip("mp2d")
using_psi4 = _build_pytest_skip("psi4")
using_qcdb = _build_pytest_skip("qcdb")
using_rdkit = _build_pytest_skip("rdkit")
using_torchani = _build_pytest_skip("torchani")
using_terachem = _build_pytest_skip("terachem")
using_cfour = _build_pytest_skip("cfour")
using_gamess = _build_pytest_skip("gamess")
using_nwchem = _build_pytest_skip("nwchem")

using_dftd3_321 = pytest.mark.skipif(is_program_new_enough("dftd3", "3.2.1") is False,
                                     reason='DFTD3 does not include 3.2.1 features. Update package and add to PATH')
