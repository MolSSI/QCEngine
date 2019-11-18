"""
Utilities for the testing suite.
"""

from typing import List

import numpy as np
import pytest
from pkg_resources import parse_version

import qcelemental as qcel
import qcengine as qcng
from qcelemental.util import which, which_import

QCENGINE_RECORDS_COMMIT = "fe14d77"


def _check_qcenginerecords(return_data=False):

    skip = True
    try:
        import qcenginerecords

        qcer_hash = qcenginerecords.__git_revision__[:7]
        if qcer_hash != QCENGINE_RECORDS_COMMIT[:7]:
            msg = f"Incorrect QCEngineRecord Git Revsion, found {qcer_hash} need {QCENGINE_RECORDS_COMMIT[:7]}."
        else:
            skip = False
            msg = "Works!"

    except ModuleNotFoundError:
        msg = "Could not find QCEngineRecords in PYTHONPATH"

    if return_data:
        return skip, msg
    else:
        return pytest.mark.skipif(skip, reason=msg)


using_qcenginerecords = _check_qcenginerecords()


def qcengine_records(program):

    skip, msg = _check_qcenginerecords(return_data=True)
    if skip:
        pytest.skip(msg, allow_module_level=True)

    import qcenginerecords

    return qcenginerecords.get_info(program)


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

        def compute(self, input_data: "AtomicInput", config: "JobConfig") -> "AtomicResult":
            self.ncalls += 1
            mode = self.iter_modes.pop(0)

            geom = input_data.molecule.geometry
            if geom.shape[0] != 2:
                raise ValueError("Failure Test must have an input size of two.")

            grad_value = np.abs(np.linalg.norm(geom[0] - geom[1]) - self.equilibrium_distance)
            grad = [0, 0, -grad_value, 0, 0, grad_value]

            if mode == "pass":
                return qcel.models.AtomicResult(
                    **{
                        **input_data.dict(),
                        **{
                            "properties": {"return_energy": grad_value},
                            "return_result": grad,
                            "success": True,
                            "extras": {"ncalls": self.ncalls},
                            "provenance": {"creator": "failure_engine", "ncores": config.ncores},
                        },
                    }
                )
            elif mode == "random_error":
                raise qcng.exceptions.RandomError("Whoops!")
            elif mode == "input_error":
                raise qcng.exceptions.InputError("Whoops!")
            else:
                raise KeyError("Testing error, should not arrive here.")

        def get_job(self):
            json_data = {
                "molecule": {"symbols": ["He", "He"], "geometry": [0, 0, 0, 0, 0, self.start_distance]},
                "driver": "gradient",
                "model": {"method": "something"},
            }

            return json_data

    engine = FailEngine()
    qcng.register_program(engine)

    yield engine

    qcng.unregister_program(engine.name)


# Figure out what is imported
_programs = {
    "dftd3": which("dftd3", return_bool=True),
    "geometric": which_import("geometric", return_bool=True),
    "psi4": is_program_new_enough("psi4", "1.2"),
    "psi4_14": is_program_new_enough("psi4", "1.4a2.dev250"),
    "qchem": is_program_new_enough("qchem", "5.2"),
    "rdkit": which_import("rdkit", return_bool=True),
    "qcdb": which_import("qcdb", return_bool=True),
    "torchani": is_program_new_enough("torchani", "0.9"),
    "mp2d": which("mp2d", return_bool=True),
    "terachem": which("terachem", return_bool=True),
    "molpro": is_program_new_enough("molpro", "2018.1"),
    "mopac": is_program_new_enough("mopac", "2016"),
    "entos": is_program_new_enough("entos", "0.6"),
    "cfour": which("xcfour", return_bool=True),
    "gamess": which("rungms", return_bool=True),
    "nwchem": which("nwchem", return_bool=True),
    "turbomole": which("define", return_bool=True),
    "mdi": which_import("mdi", return_bool=True),
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
using_psi4_14 = _build_pytest_skip("psi4_14")
using_qcdb = _build_pytest_skip("qcdb")
using_qchem = _build_pytest_skip("qchem")
using_rdkit = _build_pytest_skip("rdkit")
using_torchani = _build_pytest_skip("torchani")
using_terachem = _build_pytest_skip("terachem")
using_cfour = _build_pytest_skip("cfour")
using_gamess = _build_pytest_skip("gamess")
using_nwchem = _build_pytest_skip("nwchem")
using_turbomole = _build_pytest_skip("turbomole")
using_mdi = _build_pytest_skip("mdi")

using_dftd3_321 = pytest.mark.skipif(
    is_program_new_enough("dftd3", "3.2.1") is False,
    reason="DFTD3 does not include 3.2.1 features. Update package and add to PATH",
)
