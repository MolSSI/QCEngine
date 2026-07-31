import os
import re
from pathlib import Path

import numpy as np
import pytest
import qcengine as qcng
from qcelemental.models.v2 import AtomicInput

from qcengine.config import TaskConfig
from qcengine.programs.cclib_programs import cclib_orca
from qcengine.programs.cclib_programs.base import (
    ExecutionResult,
    _load_cclib_api,
    _parse_and_convert,
    _validate_v1_atomic_result,
)

TASK_CONFIG = TaskConfig(
    ncores=2,
    nnodes=1,
    memory=1.5,
    scratch_directory=None,
    retries=0,
    mpiexec_command=None,
)


def _input(method, basis, keywords=None, geometry=None):
    return AtomicInput(
        molecule={
            "symbols": ["H", "H"],
            "geometry": geometry or [0.0, 0.0, -0.7, 0.0, 0.0, 0.7],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        },
        specification={
            "driver": "energy",
            "model": {"method": method, "basis": basis},
            "keywords": keywords or {},
        },
    )


# One compact request per supported fixture version. These fragments cover ORCA's
# simple keywords, user blocks, TaskConfig resources, and Bohr-to-Angstrom coordinates.
INPUT_ROWS = [
    pytest.param(
        "5.0",
        _input("hf", "sto-3g", {"simple": ["TightSCF"], "blocks": {"scf": "MaxIter 150"}}),
        (r"! hf sto-3g TightSCF", r"%scf\nMaxIter 150\nend", r"nprocs 2", r"%MaxCore 768"),
        id="5.0-simple-and-block",
    ),
    pytest.param(
        "6.0",
        _input(
            "b3lyp",
            "def2-svp",
            {"simple": ["RIJCOSX"], "blocks": {"cpcm": "epsilon 78.4"}},
            geometry=[0.0, 0.0, 0.0, 0.0, 0.0, 1.4],
        ),
        (r"! b3lyp def2-svp RIJCOSX", r"%cpcm\nepsilon 78.4\nend", r"H 0.0 0.0 0.740848", r"nprocs 2"),
        id="6.0-coordinates-and-block",
    ),
    pytest.param(
        "6.1",
        _input("mp2", "def2-tzvp", {"blocks": {"scf": "ConvForced true"}}),
        (r"! mp2 def2-tzvp", r"%scf\nConvForced true\nend", r"%MaxCore 768"),
        id="6.1-method-basis-resources",
    ),
]


# ORCA 5.0 has no collector-eligible DFT output; 6.0's DFT fixtures likewise
# fail cclib QCSchema conversion. The supported HF, MP2, and CCSD fixture classes
# are represented for every version with eligible output data.
OUTPUT_ROWS = [
    pytest.param("5.0", "ORCA/basicORCA5.0/Trp_polar.out", "properties.return_energy", -673.59057112, id="5.0-hf"),
    pytest.param("5.0", "ORCA/basicORCA5.0/water_mp2.out", "properties.return_energy", -74.999373815, id="5.0-mp2"),
    pytest.param("5.0", "ORCA/basicORCA5.0/water_ccsd.out", "properties.return_energy", -75.013487814, id="5.0-ccsd"),
    pytest.param("6.0", "ORCA/basicORCA6.0/dvb_sp_hf.out", "properties.return_energy", -379.7689629142, id="6.0-hf"),
    pytest.param("6.0", "ORCA/basicORCA6.0/water_mp2.out", "properties.return_energy", -74.999374598, id="6.0-mp2"),
    pytest.param("6.0", "ORCA/basicORCA6.0/water_ccsd.out", "properties.return_energy", -75.013487814, id="6.0-ccsd"),
]


def _fixture_path(relative_output):
    source_root = os.environ.get("CCLIB_SOURCE_ROOT")
    if source_root is None:
        pytest.skip("CCLIB_SOURCE_ROOT is not set")
    fixture = Path(source_root) / "data" / relative_output
    if not fixture.is_file():
        pytest.skip(f"cclib fixture is absent: {fixture}")
    return fixture


def _fixture_input(fixture):
    api = _load_cclib_api()
    parser = api.ccopen(str(fixture))
    try:
        writer_output = api.QCSchemaWriter(parser.parse()).as_dict(validate=False)
    finally:
        parser.inputfile.close()
    initial_v1 = _validate_v1_atomic_result(writer_output)
    return AtomicInput(
        molecule=initial_v1.molecule.convert_v(2),
        specification={"driver": writer_output["driver"], "model": writer_output["model"], "keywords": {}},
    )


def _dotted_value(value, path):
    for part in path.split("."):
        value = value[part] if isinstance(value, dict) else getattr(value, part)
    return value


@pytest.mark.parametrize("version,input_model,fragments", INPUT_ROWS)
def test_generated_inputs(version, input_model, fragments):
    input_text = cclib_orca.build_input(input_model, TASK_CONFIG, "/resolved/orca").input_text
    for fragment in fragments:
        assert re.search(fragment, input_text), f"ORCA {version} input lacks {fragment!r}"


@pytest.mark.parametrize("version,relative_output,result_path,expected", OUTPUT_ROWS)
def test_parsed_outputs(version, relative_output, result_path, expected):
    fixture = _fixture_path(relative_output)
    output_text = fixture.read_text(encoding="utf-8", errors="replace")
    definition = cclib_orca.ORCA_DEFINITION
    input_model = _fixture_input(fixture)
    result = _parse_and_convert(
        definition,
        ExecutionResult(
            process_success=True,
            executable="/fixture/orca",
            input_filename=definition.input_filename,
            output_filename=definition.output_filename,
            input_text=cclib_orca.build_input(input_model, TASK_CONFIG, "/fixture/orca").input_text,
            output_text=output_text,
            stdout=output_text,
            stderr="",
        ),
        input_model,
    )
    assert np.allclose(_dotted_value(result, result_path), expected, atol=1.0e-6)


@pytest.mark.cclib_orca
def test_live_hf_single_point():
    if not qcng.get_program("cclib-orca", check=False).found():
        pytest.skip("cclib-orca is unavailable")

    result = qcng.compute(
        _input("hf", "sto-3g"),
        "cclib-orca",
        raise_error=True,
        return_version=1,
        task_config={"ncores": 1, "memory": 1.0},
    )

    assert result.success is True
    assert isinstance(result.return_result, float)
