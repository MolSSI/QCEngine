import pprint
import re

import pytest
import qcelemental as qcel
from qcelemental.models import AtomicInput

import qcengine as qcng
from qcengine.testing import has_program, using

from .test_canonical_config import _canonical_methods, _get_molecule


@pytest.mark.parametrize("native", ["none", "input", "all"])
@pytest.mark.parametrize("program, model, keywords", _canonical_methods)
def test_protocol_native(program, model, keywords, native):
    """Ensure native_files protocol implemented in harness.

    For harnesses, run minimal gradient calc with different protocol settings; check expected
      native/DSL files show up in ``AtomicResult``.

    New Harness Instructions
    ------------------------
    * Make sure minimal calc is in _canonical_methods above. This uses gradient for more files.
    * Add regex to ``input_ref`` & ``all_ref`` below to check content of input and any other file.
    * If this test doesn't work, implement or adjust ``native_files`` in your harness.

    """
    if not has_program(program):
        pytest.skip(f"Program '{program}' not found.")

    harness = qcng.get_program(program)
    molecule = _get_molecule(program, model["method"])

    #  <<  Config

    protocols = {
        "native_files": native,
    }
    config = qcng.config.get_config(
        hostname="something",
        task_config={
            "ncores": 1,
            "nnodes": 1,
        },
    )

    #  <<  Run

    inp = AtomicInput(molecule=molecule, driver="gradient", model=model, keywords=keywords, protocols=protocols)
    ret = qcng.compute(inp, program, raise_error=True, local_options=config.dict())
    pprint.pprint(ret.dict(), width=200)
    assert ret.success is True

    #  <<  Reference

    input_ref = {
        "cfour": rf"\*CFOUR\(BASIS=6-31G",
        "dftd3": rf"1.000000     1.261000     1.703000",
        "gamess": rf"\$basis gbasis=n31 ngauss=6 \$end",
        "gcp": rf"level HF3C",
        "mctc-gcp": rf"level DFT/SV",
        "mp2d": rf"--TT_a1=0.944 --TT_a2=0.48",
        "nwchem": rf"H library 6-31G",
        "psi4": rf'"driver": "gradient",',
    }

    all_ref = {
        "cfour": ("GRD", rf"1.0+\s+0.0+\s+0.0+\s+0.03"),
        "dftd3": ("dftd3_geometry.xyz", rf"H\s+0.0+\s+0.0+\s+0.34"),
        "gamess": ("gamess.dat", rf"CLOSED SHELL ORBITALS --- GENERATED AT"),
        "gcp": ("gcp_geometry.xyz", rf"H\s+0.0+\s+0.0+\s+0.34"),
        "mctc-gcp": ("gcp_geometry.xyz", rf"H\s+0.0+\s+0.0+\s+0.34"),
        "mp2d": ("mp2d_geometry", rf"H\s+0.0+\s+0.0+\s+0.34"),
        "nwchem": ("nwchem.grad", rf"0.0, 0.0, 0.03"),
        "psi4": ("psi4.grad", rf"1.0+\s+(-?)0.0+\s+(-?)0.0+\s+0.03"),
    }

    #  <<  Test

    if native == "none":
        if qcel.util.parse_version(qcel.__version__) < qcel.util.parse_version("0.25.0"):
            assert ret.native_files is None
        else:
            assert ret.native_files == {}
    elif native == "input":
        assert list(ret.native_files.keys()) == ["input"]

    if ret.native_files:
        assert "stdout" not in ret.native_files, f"Stdout found in native_files -- clean up the harness"
        assert "stderr" not in ret.native_files, f"Stderr found in native_files -- clean up the harness"
    assert "outfiles" not in ret.extras, f"Outfiles found in extras -- clean up the harness"

    if native in ["input", "all"]:
        assert re.search(
            input_ref[program], ret.native_files["input"]
        ), f"Input file pattern not found: {input_ref[program]}"
    if native == "all" and program != "psi4":  # allow psi4 once native_files PR merged
        fl, snip = all_ref[program]
        assert re.search(snip, ret.native_files[fl]), f"Ancillary file pattern not found in {fl}: {snip}"
