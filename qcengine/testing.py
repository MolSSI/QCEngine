"""
Utilities for the testing suite.
"""

import sys
from typing import List

import numpy as np
import pytest
import qcelemental as qcel
from packaging.version import parse
from pydantic import ConfigDict
from qcelemental.util import which, which_import

import qcengine as qcng

QCENGINE_RECORDS_COMMIT = "19b843b"


def _check_qcenginerecords(return_data=False):

    skip = True
    try:
        import qcenginerecords

        qcer_hash = qcenginerecords.__git_revision__[:7]
        if qcer_hash != QCENGINE_RECORDS_COMMIT[:7]:
            msg = f"Incorrect QCEngineRecord Git Revision, found {qcer_hash} need {QCENGINE_RECORDS_COMMIT[:7]}."
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
    if program in qcng.list_all_procedures():
        if program not in qcng.list_available_procedures():
            return False
        candidate_version = qcng.get_procedure(program).get_version()
    else:
        if program not in qcng.list_available_programs():
            return False
        candidate_version = qcng.get_program(program).get_version()

    return parse(candidate_version) >= parse(version_feature_introduced)


def is_mdi_new_enough(version_feature_introduced):
    if which_import("mdi", return_bool=True):
        import mdi

        candidate_version = ".".join(
            [str(mdi.MDI_MAJOR_VERSION), str(mdi.MDI_MINOR_VERSION), str(mdi.MDI_PATCH_VERSION)]
        )
        return parse(candidate_version) >= parse(version_feature_introduced)
    else:
        return False


@pytest.fixture(scope="function")
def failure_engine(schema_versions, request):
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

        model_config = ConfigDict(
            frozen=False,
        )

        @staticmethod
        def found(raise_error: bool = False) -> bool:
            return True

        def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
            self.ncalls += 1
            mode = self.iter_modes.pop(0)

            geom = input_data.molecule.geometry
            if geom.shape[0] != 2:
                raise ValueError("Failure Test must have an input size of two.")

            grad_value = np.abs(np.linalg.norm(geom[0] - geom[1]) - self.equilibrium_distance)
            grad = [0, 0, -grad_value, 0, 0, grad_value]

            if mode == "pass":
                return qcel.models.v2.AtomicResult(
                    **{
                        "input_data": input_data,
                        "molecule": input_data.molecule,
                        "properties": {"return_energy": grad_value},
                        "return_result": grad,
                        "success": True,
                        "extras": {"ncalls": self.ncalls},
                        "provenance": {"creator": "failure_engine", "ncores": config.ncores},
                    }
                )
            elif mode == "random_error":
                raise qcng.exceptions.RandomError("Whoops!")
            elif mode == "input_error":
                raise qcng.exceptions.InputError("Whoops!")
            else:
                raise KeyError("Testing error, should not arrive here.")

        def get_job(self):
            if from_v2(request.node.name):
                json_data = {
                    "molecule": {"symbols": ["He", "He"], "geometry": [0, 0, 0, 0, 0, self.start_distance]},
                    "specification": {
                        "driver": "gradient",
                        "model": {"method": "something"},
                    },
                }
            else:
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
    "adcc": is_program_new_enough("adcc", "0.15.7"),
    "cfour": which("xcfour", return_bool=True),
    "dftd3": which("dftd3", return_bool=True),
    "dftd3_321": is_program_new_enough("dftd3", "3.2.1"),
    "dftd4": which_import("dftd4", return_bool=True),
    "dftd4_350": is_program_new_enough("dftd4", "3.5.0"),
    "s-dftd3": which_import("dftd3", return_bool=True),
    "qcore": is_program_new_enough("qcore", "0.8.9"),
    "gamess": which("rungms", return_bool=True),
    "mctc-gcp": is_program_new_enough("mctc-gcp", "2.3.0"),
    "gcp": which("gcp", return_bool=True),
    "geometric": which_import("geometric", return_bool=True),
    "berny": which_import("berny", return_bool=True),
    "mdi": is_mdi_new_enough("1.2"),
    "molpro": is_program_new_enough("molpro", "2018.1"),
    "mopac": is_program_new_enough("mopac", "2016"),
    "mp2d": which("mp2d", return_bool=True),
    "nwchem": which("nwchem", return_bool=True),
    "optking": which_import("optking", return_bool=True),
    "psi4": is_program_new_enough("psi4", "1.2"),
    "psi4_runqcsk": is_program_new_enough("psi4", "1.4a2.dev160"),
    "psi4_mp2qcsk": is_program_new_enough("psi4", "1.4a2.dev580"),
    "psi4_derqcsk": is_program_new_enough("psi4", "1.5a1.dev117"),
    "qcdb": which_import("qcdb", return_bool=True),
    "qchem": is_program_new_enough("qchem", "5.1"),
    "qcmanybody": which_import("qcmanybody", return_bool=True),
    "rdkit": which_import("rdkit", return_bool=True),
    "terachem": which("terachem", return_bool=True),
    "terachem_pbs": is_program_new_enough("terachem_pbs", "0.7.2"),
    "torchani": is_program_new_enough("torchani", "0.9"),
    "torsiondrive": which_import("torsiondrive", return_bool=True),
    "turbomole": which("define", return_bool=True),
    "xtb": which_import("xtb", return_bool=True),
    "mrchem": is_program_new_enough("mrchem", "1.0.0"),
    "mace": is_program_new_enough("mace", "0.3.2"),
    "aimnet2": which_import("pyaimnet2", return_bool=True),
}
_programs["openmm"] = _programs["rdkit"] and which_import(".openmm", package="simtk", return_bool=True)


def has_program(name):
    if name in _programs:
        return _programs[name]
    else:
        raise KeyError(f"Program {name} not registered with QCEngine testing.")


_using_cache = {}


def using(program):

    if program not in _using_cache:
        import_message = f"Not detecting module {program}. Install package if necessary to enable tests."
        skip = pytest.mark.skipif(has_program(program) is False, reason=import_message)
        _using_cache[program] = skip

    return _using_cache[program]


@pytest.fixture(scope="function", params=[None, "as_v1", "as_v2", "to_v1", "to_v2"])
def schema_versions(request):
    # V1V2TEST uncomment below
    if sys.version_info >= (3, 14) and request.param != "as_v2":
        pytest.skip("Only QCSchema v2-to-v2 available for Py >=3.14")

    if request.param == "as_v1":
        return qcel.models.v1, -1, qcel.models.v1
    elif request.param == "to_v2":
        return qcel.models.v1, 2, qcel.models.v2
    elif request.param == "as_v2":
        return (qcel.models.v2, -1, qcel.models.v2)
    elif request.param == "to_v1":
        return qcel.models.v2, 1, qcel.models.v1
    else:
        return qcel.models, -1, qcel.models


@pytest.fixture(scope="function", params=["as_v1", "as_v2"])
def schema_versions2(request):
    if sys.version_info >= (3, 14) and request.param != "as_v2":
        pytest.skip("Only QCSchema v2-to-v2 available for Py >=3.14")

    if request.param == "as_v1":
        return qcel.models.v1, -1, qcel.models.v1
    elif request.param == "as_v2":
        # TODO with dict-in and dict-out and models indiscriminable and defaulting to v1
        #   the as_v2 is often not reliable, so paper over it with 2 for now. return to -1 when fixed.
        return (qcel.models.v2, 2, qcel.models.v2)


@pytest.fixture(scope="function", params=[None])
def schema_versions0(request):
    if sys.version_info >= (3, 14):
        pytest.skip("Only QCSchema v2-to-v2 available for Py >=3.14")

    return qcel.models, -1, qcel.models


@pytest.fixture(scope="function", params=[None, "as_v1", "as_v2", "to_v1", "to_v2"])
def schema_versions5(request):
    # like schema_versions except not avoiding QCSchema v1 for Py 3.14

    if request.param == "as_v1":
        return qcel.models.v1, -1, qcel.models.v1
    elif request.param == "to_v2":
        return qcel.models.v1, 2, qcel.models.v2
    elif request.param == "as_v2":
        return (qcel.models.v2, -1, qcel.models.v2)
    elif request.param == "to_v1":
        return qcel.models.v2, 1, qcel.models.v1
    else:
        return qcel.models, -1, qcel.models


def checkver_and_convert(mdl, tnm, prepost, vercheck: bool = True, cast_dict_as=None):
    import json

    import pydantic

    def check_model_v1(m):
        assert isinstance(m, pydantic.v1.BaseModel), f"type({m.__class__.__name__}) = {type(m)} ⊄ v1.BaseModel (Pyd v1)"
        assert isinstance(
            m, qcel.models.v1.basemodels.ProtoModel
        ), f"type({m.__class__.__name__}) = {type(m)} ⊄ v1.ProtoModel"
        if vercheck:
            assert m.schema_version == 1, f"{m.__class__.__name__}.schema_version = {m.schema_version} != 1"

    def check_model__v1v2(m):
        assert isinstance(
            m, (pydantic.v1.BaseModel, pydantic.BaseModel)
        ), f"type({m.__class__.__name__}) = {type(m)} ⊄ v1.BaseModel (Pyd v1)"
        assert isinstance(
            m, (qcel.models.v1.basemodels.ProtoModel, qcel.models.v2.basemodels.ProtoModel)
        ), f"type({m.__class__.__name__}) = {type(m)} ⊄ v1.ProtoModel"
        if vercheck:
            assert m.schema_version == 1, f"{m.__class__.__name__}.schema_version = {m.schema_version} != 1"

    def check_model_v2(m):
        assert isinstance(m, pydantic.BaseModel), f"type({m.__class__.__name__}) = {type(m)} ⊄ BaseModel (Pyd v2)"
        assert isinstance(
            m, qcel.models.v2.basemodels.ProtoModel
        ), f"type({m.__class__.__name__}) = {type(m)} ⊄ v2.ProtoModel"
        if vercheck:
            assert m.schema_version == 2, f"{m.__class__.__name__}.schema_version = {m.schema_version} != 2"

    if prepost == "pre":
        dict_in = isinstance(mdl, dict)
        cast_smodel = cast_dict_as or "AtomicInput"
        if "as_v1" in tnm or "to_v2" in tnm or "None" in tnm:
            if dict_in:
                mdl = getattr(qcel.models.v1, cast_smodel)(**mdl)
            # V1V2TEST check_model__v1v2(mdl)
            check_model_v1(mdl)
        elif "as_v2" in tnm or "to_v1" in tnm:
            if dict_in:
                mdl = getattr(qcel.models.v2, cast_smodel)(**mdl)
            check_model_v2(mdl)
            # NOW IN COMPUTE mdl = mdl.convert_v(1)

        if dict_in:
            mdl = mdl.model_dump()

    elif prepost == "post":
        # excuse_as_v2 is only for "post" when the compute input was a dict and while dicts are not discriminable.
        #   for now these always go to v1 in programs/model.py so as_v2 returns wrongly as v1
        # follow-up: there are too many ways this can happen, so now it's forestalled by the schema_versions fixture passing 2 to as_v2
        dict_in = isinstance(mdl, dict)
        cast_smodel = cast_dict_as or "AtomicResult"
        if "as_v1" in tnm or "to_v1" in tnm or "None" in tnm:
            if dict_in:
                mdl = getattr(qcel.models.v1, cast_smodel)(**mdl)
            # V1V2TEST check_model__v1v2(mdl)
            check_model_v1(mdl)
        elif "as_v2" in tnm or "to_v2" in tnm:
            if dict_in:
                mdl = getattr(qcel.models.v2, cast_smodel)(**mdl)
            # NOW IN COMPUTE mdl = mdl.convert_v(2)
            check_model_v2(mdl)

        if dict_in:
            # imitates compute(..., return_dict=True)
            mdl = json.loads(mdl.model_dump_json())

    return mdl


def from_v2(tnm: str) -> bool:
    """Convenience test for partitioning tests expecting QCSchema v2 input."""
    return ("to_v1" in tnm) or ("as_v2" in tnm)
