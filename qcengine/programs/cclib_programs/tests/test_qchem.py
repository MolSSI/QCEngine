import os
import re
from pathlib import Path

import numpy as np
import pytest
import qcengine as qcng
from qcelemental.models.v2 import AtomicInput

from qcengine.config import TaskConfig
from qcengine.programs.cclib_programs import cclib_qchem
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


def _input(driver, method, basis, keywords=None, geometry=None):
    return AtomicInput(
        molecule={
            "symbols": ["H", "H"],
            "geometry": geometry or [0.0, 0.0, -0.7, 0.0, 0.0, 0.7],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        },
        specification={
            "driver": driver,
            "model": {"method": method, "basis": basis},
            "keywords": keywords or {},
        },
    )


# One small request per supported fixture version; fragments cover the generator
# settings that are unique to these representatives, not complete input snapshots.
INPUT_ROWS = [
    pytest.param(
        "5.1",
        _input("energy", "hf", "sto-3g", {"scf_convergence": 8, "thresh": 14}),
        (r"JOBTYPE sp", r"METHOD hf", r"BASIS sto-3g", r"MEM_TOTAL 1536", r"SCF_CONVERGENCE 8", r"THRESH 14"),
        id="5.1-energy-scalar-keywords",
    ),
    pytest.param(
        "5.4",
        _input("gradient", "b3lyp", "6-31g(d)", geometry=[0.0, 0.0, 0.0, 0.0, 0.0, 1.4]),
        (r"JOBTYPE force", r"METHOD b3lyp", r"BASIS 6-31g\(d\)", r"INPUT_BOHR TRUE", r"H 0.0 0.0 1.4"),
        id="5.4-gradient-coordinates",
    ),
    pytest.param(
        "6.0",
        _input("hessian", "mp2", "aug-cc-pvdz"),
        (r"JOBTYPE freq", r"METHOD mp2", r"BASIS aug-cc-pvdz", r"MEM_TOTAL 1536"),
        id="6.0-hessian-resources",
    ),
]


# HF/DFT SCF, MP2, and CCSD are retained for versions with eligible fixtures.
# Q-Chem 6.0 has only eligible HF solvent fixtures in this cclib corpus.
OUTPUT_ROWS = [
    pytest.param("5.1", "QChem/basicQChem5.1/C_bigbasis.out", "properties.return_energy", -37.6045426355, id="5.1-hf"),
    pytest.param("5.1", "QChem/basicQChem5.1/dvb_sp.out", "properties.return_energy", -382.3003981057, id="5.1-dft"),
    pytest.param("5.1", "QChem/basicQChem5.1/water_mp2.out", "properties.return_energy", -75.00228214, id="5.1-mp2"),
    pytest.param("5.1", "QChem/basicQChem5.1/water_ccsd.out", "properties.return_energy", -75.01768352, id="5.1-ccsd"),
    pytest.param("5.4", "QChem/basicQChem5.4/C_bigbasis.out", "properties.return_energy", -37.6045426406, id="5.4-hf"),
    pytest.param("5.4", "QChem/basicQChem5.4/dvb_sp.out", "properties.return_energy", -382.3003975729, id="5.4-dft"),
    pytest.param("5.4", "QChem/basicQChem5.4/water_mp2.out", "properties.return_energy", -75.00228214, id="5.4-mp2"),
    pytest.param("5.4", "QChem/basicQChem5.4/water_ccsd.out", "properties.return_energy", -75.01768352, id="5.4-ccsd"),
    pytest.param(
        "6.0",
        "QChem/basicQChem6.0/water_hf_solvent_onsager.out",
        "properties.return_energy",
        -74.9645827956,
        id="6.0-hf-solvent",
    ),
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
    input_text = cclib_qchem.build_input(input_model, TASK_CONFIG, "/resolved/qchem").input_text
    for fragment in fragments:
        assert re.search(fragment, input_text), f"Q-Chem {version} input lacks {fragment!r}"


@pytest.mark.parametrize("version,relative_output,result_path,expected", OUTPUT_ROWS)
def test_parsed_outputs(version, relative_output, result_path, expected):
    fixture = _fixture_path(relative_output)
    output_text = fixture.read_text(encoding="utf-8", errors="replace")
    definition = cclib_qchem.QCHEM_DEFINITION
    input_model = _fixture_input(fixture)
    result = _parse_and_convert(
        definition,
        ExecutionResult(
            process_success=True,
            executable="/fixture/qchem",
            input_filename=definition.input_filename,
            output_filename=definition.output_filename,
            input_text=cclib_qchem.build_input(input_model, TASK_CONFIG, "/fixture/qchem").input_text,
            output_text=output_text,
            stdout="",
            stderr="",
        ),
        input_model,
    )
    assert np.allclose(_dotted_value(result, result_path), expected, atol=1.0e-6)


@pytest.mark.cclib_qchem
def test_live_hf_single_point():
    if not qcng.get_program("cclib-qchem", check=False).found():
        pytest.skip("cclib-qchem is unavailable")

    result = qcng.compute(
        _input("energy", "hf", "sto-3g"),
        "cclib-qchem",
        raise_error=True,
        return_version=1,
        task_config={"ncores": 1, "memory": 1.0},
    )

    assert result.success is True
    assert isinstance(result.return_result, float)
