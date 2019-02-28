"""
Utilities for the testing suite.
"""

import os
import shutil
import subprocess
from contextlib import contextmanager

import pytest

from .util import which


@contextmanager
def environ_context(env):
    """Temporarily set environment variables inside the context manager and
    fully restore previous environment afterwards
    """
    original_env = {key: os.getenv(key) for key in env}
    os.environ.update(env)
    try:
        yield
    finally:
        for key, value in original_env.items():
            if value is None:
                del os.environ[key]
            else:
                os.environ[key] = value


def _plugin_import(plug):
    """
    Tests to see if a module is available
    """
    import sys
    if sys.version_info >= (3, 4):
        from importlib import util
        plug_spec = util.find_spec(plug)
    else:
        import pkgutil
        plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True


def is_psi4_new_enough(version_feature_introduced):
    if not _plugin_import('psi4'):
        return False
    import psi4
    from pkg_resources import parse_version
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


def is_dftd3_new_enough(version_feature_introduced):
    if not which('dftd3', return_bool=True):
        return False
    # Note: anything below v3.2.1 will return the help menu here. but that's fine as version compare evals to False.
    command = [which('dftd3'), '-version']
    proc = subprocess.run(command, stdout=subprocess.PIPE)
    candidate_version = proc.stdout.decode('utf-8').strip()

    from pkg_resources import parse_version
    return parse_version(candidate_version) >= parse_version(version_feature_introduced)


# Figure out what is imported
_programs = {
    "dftd3": which('dftd3', return_bool=True),
    "geometric": _plugin_import("geometric"),
    "psi4": is_psi4_new_enough("1.2"),
    "rdkit": _plugin_import("rdkit"),
    "qcdb": _plugin_import("qcdb"),
    "torchani": _plugin_import("torchani"),
}

def has_program(name):
    return _programs[name]

def _build_pytest_skip(program):
    import_message = "Not detecting module {}. Install package if necessary to enable tests."
    return pytest.mark.skipif(has_program(program) is False, reason=import_message.format(program))

# Add flags
using_dftd3 = _build_pytest_skip("dftd3")
using_geometric = _build_pytest_skip("geometric")
using_psi4 = _build_pytest_skip("psi4")
using_rdkit = _build_pytest_skip("rdkit")
using_qcdb = _build_pytest_skip("qcdb")
using_torchani = _build_pytest_skip("torchani")

using_dftd3_321 = pytest.mark.skipif(
    is_dftd3_new_enough("3.2.1") is False,
    reason='DFTD3 does not include 3.2.1 features. Update package and add to PATH')
