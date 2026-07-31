import pytest
from qcelemental.models.v2 import AtomicInput

from qcengine.config import TaskConfig
from qcengine.programs.cclib_programs.cclib_orca import build_input


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
            "! hf sto-3g\n"
            "%output\n"
            "PrintLevel Normal\n"
            "Print[P_Basis] 2\n"
            "Print[P_MOs] 1\n"
            "Print[P_Overlap] 1\n"
            "Print[P_Hirshfeld] 1\n"
            "end\n"
            "%pal\n"
            "nprocs 1\n"
            "end\n"
            "%MaxCore 1024\n"
            "* xyz 0 1\n"
            "He 0.0 0.0 0.0\n"
            "*\n",
        )
    ],
)
def test_generated_inputs(version, input_model, expected_input):
    assert build_input(input_model, TASK_CONFIG, "/resolved/program").input_text == expected_input
