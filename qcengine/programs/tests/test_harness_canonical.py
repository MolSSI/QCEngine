"""
Tests the DQM compute dispatch module
"""
import pprint
import re
import tempfile
from pathlib import Path

import pytest
from qcelemental.models import AtomicInput

import qcengine as qcng
from qcengine.testing import has_program, using

_canonical_methods = [
    ("adcc", {"method": "adc2", "basis": "6-31G"}, {"n_triplets": 3}),
    ("cfour", {"method": "hf", "basis": "6-31G"}, {}),
    ("dftd3", {"method": "b3lyp-d3"}, {}),
    ("gamess", {"method": "hf", "basis": "n31"}, {"basis__NGAUSS": 6}),
    ("gcp", {"method": "hf3c"}, {}),
    ("mctc-gcp", {"method": "dft/sv"}, {}),
    ("molpro", {"method": "hf", "basis": "6-31G"}, {}),
    ("mopac", {"method": "PM6"}, {}),
    ("mp2d", {"method": "MP2-DMP2"}, {}),
    ("mrchem", {"method": "blyp"}, {"world_prec": 1.0e-3}),
    ("nwchem", {"method": "hf", "basis": "6-31G"}, {}),
    ("openmm", {"method": "openff-1.0.0", "basis": "smirnoff"}, {}),
    ("psi4", {"method": "hf", "basis": "6-31G"}, {}),
    ("qchem", {"method": "hf", "basis": "6-31G"}, {}),
    ("qcore", {"method": "pbe", "basis": "6-31G"}, {}),
    ("rdkit", {"method": "UFF"}, {}),
    ("terachem", {"method": "bad"}, {}),
    ("terachem_pbs", {"method": "b3lyp", "basis": "6-31G"}, {}),
    ("torchani", {"method": "ANI1x"}, {}),
    ("turbomole", {"method": "pbe", "basis": "6-31G"}, {}),
    ("xtb", {"method": "GFN2-xTB"}, {}),
]


def _get_molecule(program, method):
    if program in ["openmm", "terachem_pbs"]:
        return qcng.get_molecule("water")
    elif program == "gamess" and method == "ccsd(t)":
        return qcng.get_molecule("water")
    else:
        return qcng.get_molecule("hydrogen")


@pytest.mark.parametrize("memory_trickery", [True, False])
@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_local_options_memory_gib(program, model, keywords, memory_trickery):
    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    harness = qcng.get_program(program)
    molecule = _get_molecule(program, model["method"])

    # native keywords that CONTRADICT config.memory below
    misplaced_memory_spec = {
        "cfour": {"memory_size": 5},
        "nwchem": {"memory": 5},
        "gamess": {"system__mwords": 5},
        "psi4": {},  # no contradictory memory keyword in psi
    }

    if memory_trickery and harness._defaults["managed_memory"] is True:
        use_keywords = {**keywords, **misplaced_memory_spec[program]}
    else:
        use_keywords = keywords

    #  <<  Config

    config = qcng.config.get_config(
        hostname="something",
        local_options={
            "ncores": 1,
            "nnodes": 1,
            "memory": 1.555,
        },
    )

    #  <<  Run

    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=use_keywords)
    ret = qcng.compute(inp, program, raise_error=True, local_options=config.dict())
    pprint.pprint(ret.dict(), width=200)
    assert ret.success is True

    #  <<  Reference

    stdout_ref = {  # 1.555 GiB = 208708567 quad-words
        "cfour": "Allocated    1592 MB of main memory",
        "gamess": "208000000 WORDS OF MEMORY AVAILABLE",
        "nwchem": r"total    =  2087085\d\d doubles =   1592.3 Mbytes",  # doubles is quad-words. Mbytes is MiB
        "psi4": "1592 MiB Core",
    }

    #  <<  Test

    assert config.ncores == 1
    assert pytest.approx(config.memory, 0.1) == 1.555

    if harness._defaults["managed_memory"] is True:
        assert re.search(stdout_ref[program], ret.stdout), f"Memory pattern not found: {stdout_ref[program]}"


@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_local_options_scratch(program, model, keywords):
    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    harness = qcng.get_program(program)
    molecule = _get_molecule(program, model["method"])

    #  <<  Config

    scratch_directory = tempfile.mkdtemp(suffix="_" + program)

    config = qcng.config.get_config(
        hostname="something",
        local_options={
            "scratch_directory": scratch_directory,
            "scratch_messy": True,
        },
    )

    #  <<  Run

    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)
    ret = qcng.compute(inp, program, raise_error=True, local_options=config.dict())
    pprint.pprint(ret.dict(), width=200)
    assert ret.success is True

    #  <<  Reference

    stdout_ref = {
        "cfour": "University of Florida",  # freebie
        "dftd3": "Grimme",  # freebie
        "gamess": "IOWA STATE UNIVERSITY",  # freebie
        "gcp": "Grimme",  # freebie
        "mp2d": "Beran",  # freebie
        "nwchem": "E. Apra",  # freebie
        "psi4": rf"Scratch directory: {scratch_directory}/tmp\w+_psi_scratch/",
    }

    # a scratch file (preferrably output) expected after job if scratch not cleaned up
    scratch_sample = {
        "cfour": "*/NEWFOCK",
        "dftd3": "*/dftd3_geometry.xyz",  # no outfiles
        "gamess": "*/gamess.dat",
        "gcp": "*/gcp_geometry.xyz",  # no outfiles
        "mp2d": "*/mp2d_geometry",  # no outfiles
        "nwchem": "*/nwchem.db",
        "psi4": "*/psi.*.35",
    }

    #  <<  Test

    assert config.scratch_directory.endswith(program)

    if harness._defaults["scratch"] is True:
        sample_file = list(Path(scratch_directory).glob(scratch_sample[program]))
        assert len(sample_file) == 1, f"Scratch sample not found: {scratch_sample[program]} in {scratch_directory}"

        assert re.search(stdout_ref[program], ret.stdout), f"Scratch pattern not found: {stdout_ref[program]}"


@pytest.mark.parametrize("ncores", [1, 3])
@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_local_options_ncores(program, model, keywords, ncores):
    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    harness = qcng.get_program(program)
    molecule = _get_molecule(program, model["method"])

    #  <<  Config

    config = qcng.config.get_config(
        hostname="something",
        local_options={
            "ncores": ncores,
            "nnodes": 1,
        },
    )

    #  <<  Run

    inp = AtomicInput(molecule=molecule, driver="energy", model=model, keywords=keywords)
    ret = qcng.compute(inp, program, raise_error=True, local_options=config.dict())
    pprint.pprint(ret.dict(), width=200)
    assert ret.success is True

    #  <<  Reference

    stdout_ref = {
        "cfour": rf"Running with {ncores} threads/proc",
        "gamess": rf"MEMDDI DISTRIBUTED OVER\s+{ncores} PROCESSORS",
        # "gamess": rf"PARALLEL VERSION RUNNING ON\s+{ncores} PROCESSORS IN\s+1 NODES",  # no line for serial
        # nwchem is node_parallel only
        "psi4": rf"Threads:\s+{ncores}",
    }

    #  <<  Test

    assert config.ncores == ncores
    assert config.nnodes == 1

    if harness._defaults["thread_parallel"] is True:
        assert re.search(stdout_ref[program], ret.stdout), f"Thread pattern not found: {stdout_ref[program]}"
