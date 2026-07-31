import pytest
from qcelemental.models.v2 import AtomicInput

from qcengine.config import TaskConfig
from qcengine.programs.cclib_programs.cclib_qchem import build_input


HE_INPUT = AtomicInput(
    molecule={"symbols": ["He"], "geometry": [0.0, 0.0, 0.0]},
    specification={
        "driver": "energy",
        "model": {"method": "hf", "basis": "sto-3g"},
        "keywords": {},
    },
)
TASK_CONFIG = TaskConfig(
    ncores=1,
    nnodes=1,
    memory=1.0,
    scratch_directory=None,
    retries=0,
    mpiexec_command=None,
)


@pytest.mark.parametrize(
    "version,input_model,expected_input",
    [
        (
            "minimum",
            HE_INPUT,
            "$comment\n"
            "QCEngine CCLibHarness\n"
            "$end\n\n"
            "$molecule\n"
            "0 1\n"
            "He 0.0 0.0 0.0\n"
            "$end\n\n"
            "$rem\n"
            "JOBTYPE sp\n"
            "METHOD hf\n"
            "BASIS sto-3g\n"
            "MEM_TOTAL 1024\n"
            "INPUT_BOHR TRUE\n"
            "SCF_FINAL_PRINT 2\n"
            "PRINT_GENERAL_BASIS TRUE\n"
            "PRINT_ORBITALS TRUE\n"
            "MOLDEN_FORMAT FALSE\n"
            "$end\n",
        )
    ],
)
def test_generated_inputs(version, input_model, expected_input):
    assert build_input(input_model, TASK_CONFIG, "/resolved/program").input_text == expected_input
