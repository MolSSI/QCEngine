"""
Utilities for the testing suite.
"""

import os
import shutil
import subprocess
from contextlib import contextmanager
from pkg_resources import parse_version

import pytest

from qcelemental.util import which, which_import

from .programs import list_available_programs, get_program


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


def is_program_new_enough(program, version_feature_introduced):
    """Returns True if `program` registered in QCEngine, locatable in
    environment, has parseable version, and that version in normalized
    form is equal to or later than `version_feature_introduced`.

    """
    if program not in list_available_programs():
        return False
    candidate_version = get_program(program).get_version()

    return parse_version(candidate_version) >= parse_version(version_feature_introduced)


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
    "entos": is_program_new_enough("entos", "0.5")
}

def has_program(name):
    return _programs[name]

def _build_pytest_skip(program):
    import_message = "Not detecting module {}. Install package if necessary to enable tests."
    return pytest.mark.skipif(has_program(program) is False, reason=import_message.format(program))

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

using_dftd3_321 = pytest.mark.skipif(
    is_program_new_enough("dftd3", "3.2.1") is False,
    reason='DFTD3 does not include 3.2.1 features. Update package and add to PATH')