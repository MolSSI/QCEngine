"""
Utilities for the testing suite.
"""

import os
import shutil
import subprocess

import pytest
from contextlib import contextmanager
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


# Add flags
using_psi4 = pytest.mark.skipif(
    is_psi4_new_enough("1.2") is False,
    reason="Could not find psi4 or version too old. Please install the package to enable tests")

using_rdkit = pytest.mark.skipif(
    _plugin_import("rdkit") is False, reason="Could not find rdkit. Please install the package to enable tests")

using_geometric = pytest.mark.skipif(
    _plugin_import("geometric") is False,
    reason="could not find geomeTRIC. Please install the package to enable tests")

using_torchani = pytest.mark.skipif(
    _plugin_import("torchani") is False, reason="Could not find TorchAni. Please install the package to enable tests")

using_qcdb = pytest.mark.skipif(
    _plugin_import("qcdb") is False,
    reason='Not detecting common driver. Install package if necessary and add to envvar PYTHONPATH')

using_dftd3 = pytest.mark.skipif(
    which('dftd3', return_bool=True) is False,
    reason='Not detecting executable dftd3. Install package if necessary and add to envvar PATH')

using_dftd3_321 = pytest.mark.skipif(
    is_dftd3_new_enough("3.2.1") is False,
    reason='DFTD3 does not include 3.2.1 features. Update package and add to PATH')
