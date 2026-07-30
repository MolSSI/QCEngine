import json
import os
import subprocess
import sys
from dataclasses import replace
from types import SimpleNamespace

import pytest

import qcengine as qcng
import qcengine.programs.cclib as cclib_harness
from qcelemental.models.v2 import AtomicInput, BasisSet

from qcengine.config import TaskConfig
from qcengine.exceptions import InputError, ResourceError, UnknownError
from qcengine.programs.cclib import CCLibHarness
from qcengine.testing import uusing


def test_cclib_testing_registration(monkeypatch):
    from qcengine import testing

    for selector in ("cclib-qchem", "cclib-orca"):
        assert selector in testing._programs
        monkeypatch.setitem(testing._programs, selector, True)
        testing._using_cache.pop(selector, None)

        def marked_test():
            pass

        marked_test = testing.uusing(selector)(marked_test)
        marks = {mark.name: mark for mark in marked_test.pytestmark}
        assert set(marks) == {"skipif", "addon", selector}
        assert marks["skipif"].args == (False,)
        assert testing.has_program(selector) is True
        testing._using_cache.pop(selector, None)


def test_registered_instances_are_frozen_and_independent():
    qchem = qcng.get_program("CCLIB-QCHEM", check=False)
    orca = qcng.get_program("cclib-orca", check=False)

    assert type(qchem) is type(orca) is CCLibHarness
    assert (qchem.program, qchem.node_parallel) == ("qchem", False)
    assert (orca.program, orca.node_parallel) == ("orca", True)
    assert qchem is not orca
    assert qcng.get_program("qchem", check=False) is not qchem
    assert {"qchem", "cclib-qchem", "cclib-orca"} <= qcng.list_all_programs()

    with pytest.raises(Exception, match="frozen"):
        qchem.program = "orca"


def test_import_qcengine_does_not_import_external_cclib():
    script = r'''
import importlib.abc
import json
import sys

class RejectCCLib(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "cclib" or fullname.startswith("cclib."):
            raise RuntimeError(f"eager external import: {fullname}")
        return None

sys.meta_path.insert(0, RejectCCLib())
import qcengine
print(json.dumps(sorted(qcengine.list_all_programs())))
'''
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=os.fspath(os.path.dirname(qcng.__file__)),
        text=True,
        capture_output=True,
        check=True,
    )
    programs = set(json.loads(result.stdout))
    assert {"qchem", "cclib-qchem", "cclib-orca"} <= programs


def _qchem_probe_output(version="5.1"):
    return (
        "A Quantum Leap Into The Future Of Chemistry\n"
        f"Q-Chem {version} for Linux\n"
        "Thank you very much for using Q-Chem\n"
    )


def _orca_probe_output(version="6.0"):
    return (
        "                         O   R   C   A\n"
        f"Program Version {version}.0\n"
        "ORCA TERMINATED NORMALLY\n"
    )


def test_missing_executable_is_a_resource_error(monkeypatch):
    harness = CCLibHarness(name="cclib-orca", program="orca")
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: object())
    monkeypatch.setattr(cclib_harness, "which", lambda command: None, raising=False)

    assert harness.found() is False
    with pytest.raises(ResourceError, match="orca.*PATH"):
        harness.found(raise_error=True)


def _valid_qchem_environment(tmp_path):
    qc = tmp_path / "qc"
    qcaux = tmp_path / "qcaux"
    qc.mkdir()
    qcaux.mkdir()
    qcprog = tmp_path / "qcprog"
    qchem = tmp_path / "qchem"
    for executable in (qcprog, qchem):
        executable.write_text("#!/bin/sh\n")
        executable.chmod(0o755)
    return {"QC": str(qc), "QCAUX": str(qcaux), "QCPROG": str(qcprog)}, str(qchem)


@pytest.mark.parametrize(
    "changes,invalid",
    [
        ({"QC": None}, ["QC"]),
        ({"QCAUX": None}, ["QCAUX"]),
        ({"QCPROG": None}, ["QCPROG"]),
        ({"QC": "file", "QCAUX": "file", "QCPROG": "directory"}, ["QC", "QCAUX", "QCPROG"]),
    ],
)
def test_qchem_environment_preflight_aggregates_invalid_variables(tmp_path, changes, invalid):
    environment, executable = _valid_qchem_environment(tmp_path)
    ordinary_file = tmp_path / "ordinary"
    ordinary_file.write_text("data")
    ordinary_dir = tmp_path / "ordinary-dir"
    ordinary_dir.mkdir()
    replacements = {"file": str(ordinary_file), "directory": str(ordinary_dir)}
    for variable, value in changes.items():
        if value is None:
            environment.pop(variable)
        else:
            environment[variable] = replacements[value]

    with pytest.raises(ResourceError) as exc_info:
        cclib_harness._preflight_qchem(executable, environment)
    message = str(exc_info.value)
    for variable in invalid:
        assert variable in message
    assert "cclib-qchem" in message


def test_qchem_environment_preflight_requires_readable_and_executable_resources(tmp_path):
    environment, executable = _valid_qchem_environment(tmp_path)
    os.chmod(environment["QC"], 0o300)
    os.chmod(environment["QCAUX"], 0o300)
    os.chmod(environment["QCPROG"], 0o644)
    os.chmod(executable, 0o644)

    try:
        with pytest.raises(ResourceError) as exc_info:
            cclib_harness._preflight_qchem(executable, environment)
        message = str(exc_info.value)
        assert all(name in message for name in ("QC", "QCAUX", "QCPROG", "resolved qchem executable"))
    finally:
        os.chmod(environment["QC"], 0o700)
        os.chmod(environment["QCAUX"], 0o700)


def test_qchem_environment_preflight_builds_child_environment_without_inherited_qcscratch(tmp_path):
    environment, executable = _valid_qchem_environment(tmp_path)
    environment["UNCHANGED"] = "preserved"

    child = cclib_harness._preflight_qchem(executable, environment)

    assert child is not environment
    assert child["UNCHANGED"] == "preserved"
    assert child["QCSCRATCH"]
    assert "QCSCRATCH" not in environment


def test_found_lazily_loads_cclib_before_resolving_executable(monkeypatch):
    calls = []
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: calls.append("cclib") or object())
    monkeypatch.setattr(cclib_harness, "which", lambda command: calls.append(command) or None)
    harness = CCLibHarness(name="cclib-orca", program="orca")

    assert harness.found() is False
    assert calls == ["cclib", "orca"]


def test_found_checks_resources_in_required_order(monkeypatch):
    events = []
    harness = CCLibHarness(name="cclib-qchem", program="qchem")

    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: events.append("cclib") or object())
    monkeypatch.setattr(cclib_harness, "which", lambda command: events.append("path") or "/opt/qchem")
    definition = replace(
        cclib_harness._PROGRAM_DEFINITIONS["qchem"],
        preflight=lambda executable, environment: events.append("preflight") or environment,
    )
    monkeypatch.setitem(cclib_harness._PROGRAM_DEFINITIONS, "qchem", definition)
    monkeypatch.setattr(
        cclib_harness,
        "_probe_executable",
        lambda instance, executable, environment: events.append("identity/version") or "5.1",
    )

    assert harness.found(raise_error=True) is True
    assert events == ["cclib", "path", "preflight", "identity/version"]


@pytest.mark.parametrize(
    "stage,expected",
    [
        ("cclib", "cclib is unavailable"),
        ("path", "PATH"),
        ("preflight", "QCAUX"),
        ("probe", "identity"),
    ],
)
def test_found_false_suppresses_resource_failures_and_found_true_preserves_detail(monkeypatch, stage, expected):
    harness = CCLibHarness(name="cclib-qchem", program="qchem")
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: object())
    monkeypatch.setattr(cclib_harness, "which", lambda command: "/opt/qchem")
    definition = replace(
        cclib_harness._PROGRAM_DEFINITIONS["qchem"],
        preflight=lambda executable, environment: environment,
    )
    monkeypatch.setitem(cclib_harness._PROGRAM_DEFINITIONS, "qchem", definition)
    monkeypatch.setattr(cclib_harness, "_probe_executable", lambda *args: "5.1")

    if stage == "cclib":
        monkeypatch.setattr(
            cclib_harness,
            "_load_cclib_api",
            lambda: (_ for _ in ()).throw(ModuleNotFoundError("cclib is unavailable")),
        )
    elif stage == "path":
        monkeypatch.setattr(cclib_harness, "which", lambda command: None)
    elif stage == "preflight":
        failing_definition = replace(
            definition,
            preflight=lambda *args: (_ for _ in ()).throw(ResourceError("QCAUX is invalid")),
        )
        monkeypatch.setitem(cclib_harness._PROGRAM_DEFINITIONS, "qchem", failing_definition)
    else:
        monkeypatch.setattr(
            cclib_harness,
            "_probe_executable",
            lambda *args: (_ for _ in ()).throw(ResourceError("executable identity mismatch")),
        )

    assert harness.found() is False
    with pytest.raises(ResourceError, match=expected):
        harness.found(raise_error=True)


def test_available_programs_does_not_raise_for_unavailable_cclib(monkeypatch):
    monkeypatch.setattr(
        cclib_harness,
        "_load_cclib_api",
        lambda: (_ for _ in ()).throw(ModuleNotFoundError("cclib is unavailable")),
    )

    available = qcng.list_available_programs()
    assert "cclib-qchem" not in available
    assert "cclib-orca" not in available


def _atomic_input(driver="energy", method="hf", basis="sto-3g", keywords=None, molecule=None, protocols=None):
    if molecule is None:
        molecule = {
            "symbols": ["He"],
            "geometry": [0.0, 0.0, 0.0],
            "fix_com": True,
            "fix_orientation": True,
        }
    return AtomicInput(
        molecule=molecule,
        specification={
            "driver": driver,
            "model": {"method": method, "basis": basis},
            "keywords": {} if keywords is None else keywords,
            "protocols": {} if protocols is None else protocols,
        },
    )


def _task_config(ncores=4, memory=2.734375, scratch_directory="/scratch"):
    return TaskConfig(
        ncores=ncores,
        nnodes=1,
        memory=memory,
        scratch_directory=scratch_directory,
        retries=0,
        mpiexec_command=None,
    )


@pytest.mark.parametrize("driver", ["energy", "gradient", "hessian"])
@pytest.mark.parametrize("method", ["hf", "b3lyp", "bp86", "mp2", "ccsd", "HF", "B3lYp", "Mp2"])
def test_input_subset_accepts_supported_drivers_and_methods_without_changing_spelling(driver, method):
    input_model = _atomic_input(driver=driver, method=method, basis="STO-3G")

    assert cclib_harness._validate_input_subset(input_model) == (driver, method, "STO-3G")


@pytest.mark.parametrize(
    "input_model,match",
    [
        (_atomic_input(driver="properties"), "driver"),
        (_atomic_input(method="mp3"), "method"),
        (_atomic_input(basis=None), "basis"),
        (_atomic_input(basis=""), "basis"),
        (_atomic_input(basis="   "), "basis"),
        (
            _atomic_input(
                basis=BasisSet.model_construct(name="custom", center_data={}, atom_map=[]),
            ),
            "basis",
        ),
        (
            _atomic_input(
                molecule={
                    "symbols": ["He"],
                    "geometry": [0.0, 0.0, 0.0],
                    "real": [False],
                    "fix_com": True,
                    "fix_orientation": True,
                }
            ),
            "real|ghost",
        ),
    ],
)
def test_input_subset_rejects_unsupported_requests_before_generation(input_model, match):
    with pytest.raises(InputError, match=match):
        cclib_harness._validate_input_subset(input_model)


_QCHEM_RESERVED = {
    "JOBTYPE",
    "METHOD",
    "BASIS",
    "MEM_TOTAL",
    "INPUT_BOHR",
    "SCF_FINAL_PRINT",
    "PRINT_GENERAL_BASIS",
    "PRINT_ORBITALS",
    "MOLDEN_FORMAT",
}


def _water_input(driver="energy", method="mp2", basis="sto-3g", keywords=None):
    return _atomic_input(
        driver=driver,
        method=method,
        basis=basis,
        keywords=keywords,
        molecule={
            "symbols": ["O", "H", "H"],
            "geometry": [
                -0.0,
                0.0,
                0.22517858316070177,
                -1.4941103633283772,
                -0.0,
                -0.9007143324538345,
                1.4941103633283772,
                -0.0,
                -0.9007143324538345,
            ],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
            "validated": True,
            "fix_com": True,
            "fix_orientation": True,
        },
    )


@pytest.mark.parametrize("driver,jobtype", [("energy", "sp"), ("gradient", "force"), ("hessian", "freq")])
@pytest.mark.parametrize("method", ["hf", "b3lyp", "bp86", "mp2", "ccsd"])
def test_qchem_input_maps_supported_driver_and_method(driver, jobtype, method):
    job = cclib_harness._build_qchem_input(
        _atomic_input(driver=driver, method=method), _task_config(), "/opt/qchem"
    )

    assert f"JOBTYPE {jobtype}\n" in job.input_text
    assert f"METHOD {method}\n" in job.input_text


def test_qchem_input_exact_mp2_water_geometry_resources_defaults_and_job_metadata():
    job = cclib_harness._build_qchem_input(_water_input(), _task_config(), "/opt/qchem")

    expected = """$comment
QCEngine CCLibHarness
$end

$molecule
0 1
O -0.0 0.0 0.22517858316070177
H -1.4941103633283772 -0.0 -0.9007143324538345
H 1.4941103633283772 -0.0 -0.9007143324538345
$end

$rem
JOBTYPE sp
METHOD mp2
BASIS sto-3g
MEM_TOTAL 2800
INPUT_BOHR TRUE
SCF_FINAL_PRINT 2
PRINT_GENERAL_BASIS TRUE
PRINT_ORBITALS TRUE
MOLDEN_FORMAT FALSE
$end
"""
    assert job.input_text == expected
    assert job.command == ["/opt/qchem", "-nt", "4", "dispatch.in", "dispatch.out"]
    assert job.infiles == {"dispatch.in": expected}
    assert job.outfiles == ["dispatch.out"]
    assert (job.input_filename, job.output_filename, job.executable) == (
        "dispatch.in",
        "dispatch.out",
        "/opt/qchem",
    )


@pytest.mark.parametrize("charge,multiplicity", [(0, 1), (1, 2), (-1, 2)])
def test_qchem_input_emits_charge_and_multiplicity(charge, multiplicity):
    input_model = _atomic_input(
        molecule={
            "symbols": ["He"],
            "geometry": [0.0, 0.0, 0.0],
            "molecular_charge": charge,
            "molecular_multiplicity": multiplicity,
            "fix_com": True,
            "fix_orientation": True,
        }
    )
    job = cclib_harness._build_qchem_input(input_model, _task_config(), "/opt/qchem")

    assert f"$molecule\n{charge} {multiplicity}\n" in job.input_text


def test_qchem_input_renders_and_sorts_supported_scalar_keywords():
    keywords = {
        "zeta": "value",
        "a_bool_true": True,
        "bool_false": False,
        "an_int": 7,
        "a_float": 1.25,
    }
    job = cclib_harness._build_qchem_input(_atomic_input(keywords=keywords), _task_config(), "/opt/qchem")

    ordinary = [
        "AN_INT 7",
        "A_BOOL_TRUE TRUE",
        "A_FLOAT 1.25",
        "BOOL_FALSE FALSE",
        "ZETA value",
    ]
    positions = [job.input_text.index(line) for line in ordinary]
    assert positions == sorted(positions)


@pytest.mark.parametrize(
    "keywords",
    [
        {"bad": "line one\nline two"},
        {"bad\nkey": "value"},
        {"bad": None},
        {"bad": [1]},
        {"bad": {"nested": 1}},
        {"bad": float("nan")},
        {"bad": float("inf")},
        {"bad": float("-inf")},
    ],
)
def test_qchem_input_rejects_malformed_keyword_values(keywords):
    with pytest.raises(InputError, match="keyword"):
        cclib_harness._build_qchem_input(_atomic_input(keywords=keywords), _task_config(), "/opt/qchem")


@pytest.mark.parametrize("reserved", sorted(_QCHEM_RESERVED))
def test_qchem_input_rejects_reserved_keywords_case_insensitively(reserved):
    with pytest.raises(InputError, match="reserved"):
        cclib_harness._build_qchem_input(
            _atomic_input(keywords={reserved.swapcase(): "user"}), _task_config(), "/opt/qchem"
        )


@pytest.mark.parametrize("reserved", sorted(_QCHEM_RESERVED))
@pytest.mark.parametrize("alias", [" {}", "{} ", "\t{}", "{}\t"])
def test_qchem_input_keyword_name_rejects_whitespace_aliases_for_every_reserved_key(reserved, alias):
    malformed_key = alias.format(reserved.swapcase())

    with pytest.raises(InputError, match="keyword name"):
        cclib_harness._build_qchem_input(
            _atomic_input(keywords={malformed_key: "user"}), _task_config(), "/opt/qchem"
        )


@pytest.mark.parametrize(
    "malformed_key",
    [
        "",
        " ",
        "\t",
        "ordinary key",
        "ordinary\tkey",
        " ordinary",
        "ordinary ",
        "ordinary\x00key",
        "ordinary\x1fkey",
        "ordinary\x7fkey",
    ],
)
def test_qchem_input_keyword_name_rejects_empty_whitespace_and_control_characters(malformed_key):
    with pytest.raises(InputError, match="keyword name"):
        cclib_harness._build_qchem_input(
            _atomic_input(keywords={malformed_key: "value"}), _task_config(), "/opt/qchem"
        )


def test_qchem_input_rejects_case_insensitive_user_keyword_collisions():
    with pytest.raises(InputError, match="collision"):
        cclib_harness._build_qchem_input(
            _atomic_input(keywords={"thresh": 8, "THRESH": 10}), _task_config(), "/opt/qchem"
        )


def test_qchem_input_build_input_delegates_to_job(monkeypatch):
    harness = CCLibHarness(name="cclib-qchem", program="qchem")
    monkeypatch.setattr(cclib_harness, "which", lambda executable: "/resolved/qchem")

    built = harness.build_input(_atomic_input(), _task_config(scratch_directory="/tmp/scratch"))

    assert built == {
        "commands": ["/resolved/qchem", "-nt", "4", "dispatch.in", "dispatch.out"],
        "infiles": {"dispatch.in": built["infiles"]["dispatch.in"]},
        "outfiles": ["dispatch.out"],
        "scratch_directory": "/tmp/scratch",
    }


@pytest.mark.parametrize("driver,keyword", [("energy", None), ("gradient", "engrad"), ("hessian", "freq")])
@pytest.mark.parametrize("method", ["hf", "b3lyp", "bp86", "mp2", "ccsd"])
def test_orca_input_maps_supported_driver_and_method(driver, keyword, method):
    job = cclib_harness._build_orca_input(
        _atomic_input(driver=driver, method=method), _task_config(), "/opt/orca"
    )

    expected_line = f"! {method} sto-3g" + (f" {keyword}" if keyword else "")
    assert job.input_text.splitlines()[0] == expected_line


def test_orca_input_resources_defaults_geometry_charge_order_and_job_metadata():
    input_model = _atomic_input(
        method="B3lYp",
        basis="Def2-SVP",
        molecule={
            "symbols": ["He", "H"],
            "geometry": [0.0, 0.0, 0.0, 1.0, -2.0, 3.0],
            "molecular_charge": 1,
            "molecular_multiplicity": 1,
            "validated": True,
            "fix_com": True,
            "fix_orientation": True,
        },
    )
    job = cclib_harness._build_orca_input(input_model, _task_config(), "/opt/orca")

    expected_prefix = """! B3lYp Def2-SVP
%output
PrintLevel Normal
Print[P_Basis] 2
Print[P_MOs] 1
Print[P_Overlap] 1
Print[P_Hirshfeld] 1
end
%pal
nprocs 4
end
%MaxCore 700
* xyz 1 1
He 0.0 0.0 0.0
"""
    assert job.input_text.startswith(expected_prefix)
    coordinate_line = job.input_text.splitlines()[-2].split()
    assert coordinate_line[0] == "H"
    assert [float(value) for value in coordinate_line[1:]] == pytest.approx(
        [0.52917721067, -1.05835442134, 1.58753163201], abs=1.0e-12
    )
    assert job.input_text.endswith("*\n")
    assert job.command == ["/opt/orca", "dispatch.inp"]
    assert job.infiles == {"dispatch.inp": job.input_text}
    assert job.outfiles == []
    assert (job.input_filename, job.output_filename, job.executable) == (
        "dispatch.inp",
        "dispatch.out",
        "/opt/orca",
    )


def test_orca_input_preserves_simple_order_sorts_blocks_and_appends_output_body():
    keywords = {
        "simple": ["rks", "usesym", "TightSCF"],
        "blocks": {
            "scf": "MaxIter 200",
            "output": "PrintLevel Mini\nPrint[P_AtCharges_M] 1",
            "basis": "NewGTO H \"def2-TZVP\" end",
        },
    }
    job = cclib_harness._build_orca_input(
        _atomic_input(driver="gradient", keywords=keywords), _task_config(), "/opt/orca"
    )

    assert job.input_text.splitlines()[0] == "! hf sto-3g engrad rks usesym TightSCF"
    defaults_end = job.input_text.index("Print[P_Hirshfeld] 1")
    user_output = job.input_text.index("PrintLevel Mini")
    basis = job.input_text.index("%basis")
    scf = job.input_text.index("%scf")
    resources = job.input_text.index("%pal")
    assert defaults_end < user_output < basis < scf < resources


@pytest.mark.parametrize(
    "body",
    [
        "end\n* xyz 9 1\nH 0 0 0\n*",
        "  EnD trailing text  ",
        "\tEND # comment",
        "%pal\nnprocs 99\nend",
        "  %MAXCORE 9999 # comment",
        "%coords\nctyp xyz",
        " %scf\nMaxIter 999",
        "$new_job",
        "  $NEW_JOB # comment",
        " * xyz 0 1",
        "  * # coordinate delimiter",
    ],
)
def test_orca_input_rejects_block_body_outer_syntax_before_execution(body):
    with pytest.raises(InputError, match="block body|outer syntax"):
        cclib_harness._build_orca_input(
            _atomic_input(keywords={"blocks": {"output": body}}), _task_config(), "/opt/orca"
        )


@pytest.mark.parametrize(
    "body",
    [
        "MaxIter 200 # weekend schedule",
        "Print[ P_Basis ] 2 # percentage %pal is only text here",
        "SomeValue prefix_end suffix",
        "Label job$new_job_backup",
        "NewGTO H \"def2-TZVP\" end",
    ],
)
def test_orca_input_allows_reserved_substrings_that_are_not_first_tokens(body):
    job = cclib_harness._build_orca_input(
        _atomic_input(keywords={"blocks": {"scf": body}}), _task_config(), "/opt/orca"
    )

    assert body in job.input_text
    assert job.input_text.index(body) < job.input_text.index("%pal")


@pytest.mark.parametrize(
    "keywords",
    [
        {"unknown": []},
        {"simple": "rks"},
        {"simple": [""]},
        {"simple": ["   "]},
        {"simple": ["rks\n* xyz 9 9"]},
        {"simple": [None]},
        {"blocks": "output"},
        {"blocks": []},
        {"blocks": {"bad-name": "value"}},
        {"blocks": {"1bad": "value"}},
        {"blocks": {"scf": None}},
        {"blocks": {"scf": ["MaxIter 200"]}},
    ],
)
def test_orca_input_rejects_malformed_keyword_structures(keywords):
    with pytest.raises(InputError, match="ORCA"):
        cclib_harness._build_orca_input(_atomic_input(keywords=keywords), _task_config(), "/opt/orca")


@pytest.mark.parametrize("block", ["pal", "PAL", "maxcore", "MaxCore", "coords", "COORDS"])
def test_orca_input_rejects_reserved_resource_and_coordinate_blocks(block):
    with pytest.raises(InputError, match="reserved|coordinate"):
        cclib_harness._build_orca_input(
            _atomic_input(keywords={"blocks": {block: "user body"}}), _task_config(), "/opt/orca"
        )


@pytest.mark.parametrize("keywords", [{"coordinates": []}, {"xyz": "0 1"}, {"coords": {}}])
def test_orca_input_rejects_top_level_coordinate_injection(keywords):
    with pytest.raises(InputError, match="unknown|coordinate"):
        cclib_harness._build_orca_input(_atomic_input(keywords=keywords), _task_config(), "/opt/orca")


def test_orca_input_build_input_delegates_to_job(monkeypatch):
    harness = CCLibHarness(name="cclib-orca", program="orca")
    monkeypatch.setattr(cclib_harness, "which", lambda executable: "/resolved/orca")

    built = harness.build_input(_atomic_input(), _task_config(scratch_directory="/tmp/scratch"))

    assert built == {
        "commands": ["/resolved/orca", "dispatch.inp"],
        "infiles": {"dispatch.inp": built["infiles"]["dispatch.inp"]},
        "outfiles": [],
        "scratch_directory": "/tmp/scratch",
    }


def test_qchem_execution_selection_requests_primary_output_and_preserves_job(monkeypatch, tmp_path):
    job = cclib_harness._build_qchem_input(_atomic_input(), _task_config(), "/resolved/qchem")
    definition = replace(cclib_harness._PROGRAM_DEFINITIONS["qchem"], preflight=lambda exe, env: dict(env))
    inherited_environment = os.environ.copy()
    calls = []

    def fake_execute(command, infiles=None, outfiles=None, **kwargs):
        calls.append((command, infiles, outfiles, kwargs))
        return True, {
            "stdout": "launcher stdout",
            "stderr": "launcher stderr",
            "outfiles": {"dispatch.out": _qchem_probe_output()},
        }

    monkeypatch.setattr(cclib_harness, "execute", fake_execute)
    result = cclib_harness._execute_job(
        definition, job, _task_config(scratch_directory=str(tmp_path))
    )

    command, infiles, outfiles, kwargs = calls[0]
    assert command == job.command
    assert infiles == job.infiles
    assert outfiles == ["dispatch.out"]
    assert {key: kwargs["environment"][key] for key in inherited_environment} == inherited_environment
    assert result.output_text == _qchem_probe_output()
    assert (result.stdout, result.stderr) == ("launcher stdout", "launcher stderr")


def test_orca_execution_selection_uses_captured_stdout_and_preserves_job(monkeypatch, tmp_path):
    job = cclib_harness._build_orca_input(_atomic_input(), _task_config(), "/resolved/orca")
    definition = cclib_harness._PROGRAM_DEFINITIONS["orca"]
    inherited_environment = os.environ.copy()
    calls = []
    output = _orca_probe_output()

    def fake_execute(command, infiles=None, outfiles=None, **kwargs):
        calls.append((command, infiles, outfiles, kwargs))
        return True, {"stdout": output, "stderr": "", "outfiles": {}}

    monkeypatch.setattr(cclib_harness, "execute", fake_execute)
    result = cclib_harness._execute_job(
        definition, job, _task_config(scratch_directory=str(tmp_path))
    )

    command, infiles, outfiles, kwargs = calls[0]
    assert (command, infiles, outfiles) == (job.command, job.infiles, [])
    assert kwargs["environment"] == inherited_environment
    assert result.output_text == output
    assert result.stdout == output


def test_qchem_managed_scratch_overrides_inherited_qcscratch(monkeypatch, tmp_path):
    monkeypatch.setenv("QCSCRATCH", "/inherited/unmanaged")
    job = cclib_harness._build_qchem_input(_atomic_input(), _task_config(), "/resolved/qchem")
    definition = replace(cclib_harness._PROGRAM_DEFINITIONS["qchem"], preflight=lambda exe, env: dict(env))
    observed = {}

    def fake_execute(command, infiles=None, outfiles=None, **kwargs):
        qcscratch = kwargs["environment"]["QCSCRATCH"]
        observed["qcscratch"] = qcscratch
        assert os.path.isdir(qcscratch)
        assert os.path.commonpath([qcscratch, str(tmp_path)]) == str(tmp_path)
        assert kwargs["scratch_directory"] == qcscratch
        return True, {"stdout": "", "stderr": "", "outfiles": {"dispatch.out": _qchem_probe_output()}}

    monkeypatch.setattr(cclib_harness, "execute", fake_execute)
    cclib_harness._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    assert not os.path.exists(observed["qcscratch"])


def _qchem_execution_case():
    job = cclib_harness._build_qchem_input(_atomic_input(), _task_config(), "/resolved/qchem")
    definition = replace(cclib_harness._PROGRAM_DEFINITIONS["qchem"], preflight=lambda exe, env: dict(env))
    return definition, job


def _assert_execution_message(error, definition, job, stage, diagnostic):
    message = str(error)
    assert definition.selector in message
    assert job.executable in message
    assert stage in message
    assert cclib_harness._diagnostic_tail(diagnostic) in message


def test_diagnostic_tail_enforces_exact_line_and_character_bounds():
    lines = [f"diagnostic line {index}: " + ("x" * 120) for index in range(50)]
    text = "\n".join(lines)

    expected = "\n".join(lines[-40:])[-4000:]
    tail = cclib_harness._diagnostic_tail(text)

    assert tail == expected
    assert len(tail) == 4000
    assert len(tail.splitlines()) <= 40
    assert cclib_harness._diagnostic_tail("short output") == "short output"
    assert cclib_harness._diagnostic_tail("") == ""


def test_nonzero_execution_error_is_unknown_with_bounded_diagnostic(monkeypatch, tmp_path):
    definition, job = _qchem_execution_case()
    diagnostic = "\n".join(f"failure line {index}" for index in range(60))
    output = diagnostic + "\n" + definition.normal_termination
    monkeypatch.setattr(
        cclib_harness,
        "execute",
        lambda *args, **kwargs: (False, {"stdout": "", "stderr": "", "outfiles": {"dispatch.out": output}}),
    )

    with pytest.raises(cclib_harness.UnknownError) as exc_info:
        cclib_harness._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    _assert_execution_message(exc_info.value, definition, job, "execution", output)


@pytest.mark.parametrize("program", ["qchem", "orca"])
def test_missing_primary_output_is_unknown_output_selection_error(monkeypatch, tmp_path, program):
    diagnostic = f"{program} launcher produced no primary output"
    if program == "qchem":
        definition, job = _qchem_execution_case()
        process = {"stdout": diagnostic, "stderr": "", "outfiles": {"dispatch.out": None}}
    else:
        definition = cclib_harness._PROGRAM_DEFINITIONS["orca"]
        job = cclib_harness._build_orca_input(_atomic_input(), _task_config(), "/resolved/orca")
        process = {"stderr": diagnostic, "outfiles": {}}
    monkeypatch.setattr(cclib_harness, "execute", lambda *args, **kwargs: (True, process))

    with pytest.raises(cclib_harness.UnknownError) as exc_info:
        cclib_harness._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    _assert_execution_message(exc_info.value, definition, job, "output selection", diagnostic)
    assert exc_info.value.__cause__ is not None


def test_zero_exit_without_normal_termination_marker_is_unknown(monkeypatch, tmp_path):
    definition, job = _qchem_execution_case()
    output = "Q-Chem stopped before its farewell"
    monkeypatch.setattr(
        cclib_harness,
        "execute",
        lambda *args, **kwargs: (True, {"stdout": "", "stderr": "", "outfiles": {"dispatch.out": output}}),
    )

    with pytest.raises(cclib_harness.UnknownError) as exc_info:
        cclib_harness._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    _assert_execution_message(exc_info.value, definition, job, "termination", output)


@pytest.mark.parametrize(
    "diagnostic,variable",
    [
        ("Undefined environment variable QCAUX", "QCAUX"),
        ("QCFILE: Undefined variable.", "QCFILE"),
        ("Environment variable 'QC' must be defined", "QC"),
    ],
)
def test_qchem_undefined_environment_execution_error_is_resource_error(
    monkeypatch, tmp_path, diagnostic, variable
):
    definition, job = _qchem_execution_case()
    monkeypatch.setattr(
        cclib_harness,
        "execute",
        lambda *args, **kwargs: (False, {"stdout": "", "stderr": diagnostic, "outfiles": {"dispatch.out": None}}),
    )

    with pytest.raises(ResourceError) as exc_info:
        cclib_harness._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    _assert_execution_message(exc_info.value, definition, job, "execution", diagnostic)
    assert variable in str(exc_info.value)


@pytest.mark.parametrize("license_message", ["FlexNet failure", "license checkout failed", "unable to validate license"])
def test_license_execution_error_is_resource_error(monkeypatch, tmp_path, license_message):
    definition = cclib_harness._PROGRAM_DEFINITIONS["orca"]
    job = cclib_harness._build_orca_input(_atomic_input(), _task_config(), "/resolved/orca")
    monkeypatch.setattr(
        cclib_harness,
        "execute",
        lambda *args, **kwargs: (False, {"stdout": license_message, "stderr": "", "outfiles": {}}),
    )

    with pytest.raises(ResourceError) as exc_info:
        cclib_harness._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    _assert_execution_message(exc_info.value, definition, job, "execution", license_message)


def test_execution_exception_is_unknown_and_chained(monkeypatch, tmp_path):
    definition = cclib_harness._PROGRAM_DEFINITIONS["orca"]
    job = cclib_harness._build_orca_input(_atomic_input(), _task_config(), "/resolved/orca")
    original = OSError("scheduler launch failed")

    def fail_execute(*args, **kwargs):
        raise original

    monkeypatch.setattr(cclib_harness, "execute", fail_execute)
    with pytest.raises(cclib_harness.UnknownError) as exc_info:
        cclib_harness._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    _assert_execution_message(exc_info.value, definition, job, "execution", str(original))
    assert exc_info.value.__cause__ is original


def _fake_writer_output(driver="energy", method="hf", basis="sto-3g"):
    return {
        "schema_name": "qcschema_output",
        "schema_version": 1,
        "molecule": {
            "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 1.4],
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
            "schema_name": "qcschema_molecule",
            "schema_version": 2,
            "symbols": ["H", "H"],
            "validated": True,
        },
        "provenance": {"creator": "Q-Chem", "version": "6.2", "routine": "cclib.QCSchemaWriter"},
        "success": True,
        "error": None,
        "stdout": None,
        "stderr": None,
        "extras": {
            "atomcoords": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]]],
            "atomnos": [1, 1],
            "charge": 0,
            "metadata": {"package": "Q-Chem", "success": True},
            "mult": 1,
            "natom": 2,
            "nbasis": 2,
            "nmo": 2,
            "scfenergies": [-1.0],
        },
        "driver": driver,
        "keywords": {},
        "model": {"method": method, "basis": basis},
        "properties": {
            "calcinfo_nalpha": 1,
            "calcinfo_natom": 2,
            "calcinfo_nbasis": 2,
            "calcinfo_nbeta": 1,
            "calcinfo_nmo": 2,
            "return_energy": -1.0,
            "scf_total_energy": -1.0,
        },
        "return_result": -1.0,
    }


def _fake_conversion_case(
    writer_output=None, *, parser_failure=None, parsed=None, detected=True, mismatch=False, program="qchem"
):
    if parsed is None:
        parsed = SimpleNamespace(metadata={"package": "Q-Chem", "package_version": "6.2", "success": True})
    opened = []

    class FakeQChem:
        def __init__(self, source):
            self.source = source

        def parse(self):
            if parser_failure is not None:
                raise parser_failure
            return parsed

    class FakeORCA(FakeQChem):
        pass

    parser_class = FakeQChem if program == "qchem" else FakeORCA

    class WrongParser(parser_class):
        pass

    def ccopen(source):
        opened.append(source.read())
        source.seek(0)
        if not detected:
            return None
        return WrongParser(source) if mismatch else parser_class(source)

    class FakeWriter:
        def __init__(self, data):
            assert data is parsed

        def as_dict(self, validate=True):
            assert validate is False
            if isinstance(writer_output, BaseException):
                raise writer_output
            return _fake_writer_output() if writer_output is None else writer_output

    api = cclib_harness._CCLibAPI(
        version="1.9.fake",
        QCSchemaWriter=FakeWriter,
        ccopen=ccopen,
    )
    execution = cclib_harness._ExecutionResult(
        process_success=True,
        executable=f"/resolved/{program}",
        input_filename="dispatch.in" if program == "qchem" else "dispatch.inp",
        output_filename="dispatch.out",
        input_text="complete native input",
        output_text="complete parsed output",
        stdout="launcher stdout",
        stderr="launcher stderr",
    )
    definition = replace(cclib_harness._PROGRAM_DEFINITIONS[program], parser_type=lambda: parser_class)
    return api, definition, execution, opened


def _assert_conversion_failure(monkeypatch, match, **case):
    api, definition, execution, _ = _fake_conversion_case(**case)
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)
    with pytest.raises(UnknownError, match=match) as exc_info:
        cclib_harness._parse_and_convert(definition, execution, _atomic_input())
    assert definition.selector in str(exc_info.value)
    assert exc_info.value.__cause__ is not None
    return exc_info.value


def _track_parser_temporary_file(monkeypatch):
    real_named_temporary_file = cclib_harness.tempfile.NamedTemporaryFile
    state = {}

    class TrackedTemporaryFile:
        def __init__(self, temporary):
            self._temporary = temporary
            self.name = temporary.name

        @property
        def closed(self):
            return self._temporary.closed

        def write(self, value):
            return self._temporary.write(value)

        def flush(self):
            return self._temporary.flush()

        def close(self):
            return self._temporary.close()

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.close()

    def tracked_named_temporary_file(*args, **kwargs):
        state["kwargs"] = dict(kwargs)
        kwargs["delete"] = False
        temporary = TrackedTemporaryFile(real_named_temporary_file(*args, **kwargs))
        state["temporary"] = temporary
        state["path"] = temporary.name
        return temporary

    monkeypatch.setattr(cclib_harness.tempfile, "NamedTemporaryFile", tracked_named_temporary_file)
    return state


def test_parser_temporary_file_is_closed_before_reopen_and_deleted_after_success(monkeypatch):
    api, definition, execution, _ = _fake_conversion_case()
    state = _track_parser_temporary_file(monkeypatch)
    original_ccopen = api.ccopen

    def assert_closed_before_reopen(source):
        assert state["temporary"].closed
        assert source.name == state["path"]
        return original_ccopen(source)

    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: replace(api, ccopen=assert_closed_before_reopen))

    result = cclib_harness._parse_and_convert(definition, execution, _atomic_input())

    assert result.success is True
    assert state["kwargs"]["suffix"] == ".out"
    assert state["kwargs"]["delete"] is False
    assert not os.path.exists(state["path"])


@pytest.mark.parametrize(
    "case,match",
    [
        ({"mismatch": True}, "parser identity"),
        ({"parser_failure": ValueError("parser exploded")}, "parser parse"),
        ({"writer_output": RuntimeError("writer exploded")}, "QCSchema writer"),
    ],
)
def test_parser_temporary_file_is_deleted_after_parse_or_conversion_failure(monkeypatch, case, match):
    api, definition, execution, _ = _fake_conversion_case(**case)
    state = _track_parser_temporary_file(monkeypatch)
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)

    with pytest.raises(UnknownError, match=match):
        cclib_harness._parse_and_convert(definition, execution, _atomic_input())

    assert state["temporary"].closed
    assert not os.path.exists(state["path"])


def test_parser_auto_detection_failure_is_bounded_stage_aware_and_chained(monkeypatch):
    api, definition, execution, _ = _fake_conversion_case(detected=False)
    execution = replace(execution, output_text="\n".join(f"line {index}" for index in range(1000)))
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)

    with pytest.raises(UnknownError, match="parser auto-detection") as exc_info:
        cclib_harness._parse_and_convert(definition, execution, _atomic_input())

    assert "line 0" not in str(exc_info.value)
    assert "line 999" in str(exc_info.value)
    assert exc_info.value.__cause__ is not None


def test_parser_class_mismatch_is_rejected(monkeypatch):
    error = _assert_conversion_failure(monkeypatch, "parser identity", mismatch=True)
    assert "QChem" in str(error)


def test_parser_exception_is_unknown_and_chained(monkeypatch):
    original = ValueError("parser exploded")
    error = _assert_conversion_failure(monkeypatch, "parser parse", parser_failure=original)
    assert error.__cause__ is original


@pytest.mark.parametrize("metadata", [{}, {"success": False}, {"success": None}])
def test_parser_incomplete_result_is_rejected(monkeypatch, metadata):
    _assert_conversion_failure(monkeypatch, "parser result validation", parsed=SimpleNamespace(metadata=metadata))


@pytest.mark.parametrize("failure", [RuntimeError("writer exploded"), NotImplementedError("unsupported method")])
def test_writer_failure_is_unknown_and_chained(monkeypatch, failure):
    error = _assert_conversion_failure(monkeypatch, "writer", writer_output=failure)
    assert error.__cause__ is failure


def test_incomplete_writer_fields_are_rejected_before_validation(monkeypatch):
    output = _fake_writer_output()
    output.pop("properties")
    _assert_conversion_failure(monkeypatch, "writer output", writer_output=output)


def test_writer_cclib_harness_extra_collision_is_rejected_without_overwrite(monkeypatch):
    output = _fake_writer_output()
    existing = {"writer_owned": "must survive"}
    output["extras"]["cclib_harness"] = existing
    api, definition, execution, _ = _fake_conversion_case(output)
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)

    with pytest.raises(UnknownError, match="writer augmentation") as exc_info:
        cclib_harness._parse_and_convert(definition, execution, _atomic_input())

    assert output["extras"]["cclib_harness"] is existing
    assert output["extras"]["cclib_harness"] == {"writer_owned": "must survive"}
    assert isinstance(exc_info.value.__cause__, ValueError)
    assert "cclib_harness" in str(exc_info.value.__cause__)
    assert definition.selector in str(exc_info.value)


def test_non_mapping_writer_extras_is_bounded_chained_augmentation_error(monkeypatch):
    output = _fake_writer_output()
    output["extras"] = 7
    api, definition, execution, _ = _fake_conversion_case(output)
    execution = replace(execution, output_text="\n".join(f"writer line {index}" for index in range(100)))
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)

    with pytest.raises(UnknownError, match="writer augmentation") as exc_info:
        cclib_harness._parse_and_convert(definition, execution, _atomic_input())

    assert isinstance(exc_info.value.__cause__, TypeError)
    assert "writer line 0" not in str(exc_info.value)
    assert "writer line 99" in str(exc_info.value)
    assert definition.selector in str(exc_info.value)


def test_qcschema_v1_validation_failure_is_unknown_and_chained(monkeypatch):
    output = _fake_writer_output()
    output["return_result"] = "not-an-energy"
    _assert_conversion_failure(monkeypatch, "QCSchema v1 validation", writer_output=output)


@pytest.mark.parametrize(
    "requested,parsed",
    [("hf", "RHF"), ("hf", "UHF"), ("mp2", "RMP2"), ("mp2", "UMP2"), ("ccsd", "RCCSD"), ("ccsd", "UCCSD")],
)
def test_successful_conversion_preserves_identity_data_and_aliases(monkeypatch, requested, parsed):
    output = _fake_writer_output(driver="energy", method=parsed, basis="STO-3G")
    api, definition, execution, opened = _fake_conversion_case(output)
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)
    input_model = _atomic_input(method=requested, basis="sto-3g")

    result = cclib_harness._parse_and_convert(definition, execution, input_model)

    assert result.schema_version == 2
    assert result.input_data == input_model
    assert list(result.molecule.symbols) == ["H", "H"]
    assert result.stdout == execution.output_text
    assert result.provenance.creator == "Q-Chem"
    assert result.provenance.version == "6.2"
    assert opened == [execution.output_text]
    assert set(result.extras) == set(output["extras"]) | {"cclib_harness"}
    for key, value in output["extras"].items():
        assert result.extras[key] == value
    assert result.extras["cclib_harness"] == {
        "selector": definition.selector,
        "cclib_version": "1.9.fake",
        "parser": definition.expected_parser,
        "executable": execution.executable,
    }


def test_orca_conversion_preserves_writer_provenance_and_single_dispersion_value(monkeypatch):
    output = _fake_writer_output()
    output["provenance"] = {"creator": "ORCA", "version": "6.0.1", "routine": "cclib.QCSchemaWriter"}
    output["extras"]["dispersionenergies"] = [-0.001]
    api, definition, execution, opened = _fake_conversion_case(output, program="orca")
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)

    result = cclib_harness._parse_and_convert(definition, execution, _atomic_input())

    assert opened == [execution.output_text]
    assert result.provenance.creator == "ORCA"
    assert result.provenance.version == "6.0.1"
    assert result.stdout == execution.output_text
    assert result.extras["dispersionenergies"] == [-0.001]
    assert list(result.extras).count("dispersionenergies") == 1
    assert result.extras["cclib_harness"]["parser"] == "ORCA"


def test_basis_identity_comparison_trims_and_casefolds_without_relabeling(monkeypatch):
    requested_basis = "  sTo-3G  "
    parsed_basis = "\tSTo-3g "
    output = _fake_writer_output(basis=parsed_basis)
    api, definition, execution, _ = _fake_conversion_case(output)
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)
    original_validate = cclib_harness._validate_v1_atomic_result
    validated_output = {}

    def capture_validated_output(value):
        validated_output.update(value)
        return original_validate(value)

    monkeypatch.setattr(cclib_harness, "_validate_v1_atomic_result", capture_validated_output)
    input_model = _atomic_input(basis=requested_basis)

    result = cclib_harness._parse_and_convert(definition, execution, input_model)

    assert validated_output["model"]["basis"] == parsed_basis
    assert result.input_data.specification.model.basis == requested_basis


def test_basis_identity_real_mismatch_after_trimming_is_rejected(monkeypatch):
    output = _fake_writer_output(basis="  6-31G  ")
    api, definition, execution, _ = _fake_conversion_case(output)
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)

    with pytest.raises(UnknownError, match="basis mismatch"):
        cclib_harness._parse_and_convert(definition, execution, _atomic_input(basis="  STO-3G  "))


@pytest.mark.parametrize(
    "field,requested,parsed",
    [
        ("driver", "energy", "gradient"),
        ("method", "hf", "b3lyp"),
        ("method", "hf", "ROHF"),
        ("basis", "sto-3g", "6-31g"),
    ],
)
def test_parsed_identity_mismatch_is_rejected_without_relabeling(monkeypatch, field, requested, parsed):
    values = {"driver": "energy", "method": "hf", "basis": "sto-3g"}
    values[field] = parsed
    output = _fake_writer_output(**values)
    api, definition, execution, _ = _fake_conversion_case(output)
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)
    inputs = {"driver": "energy", "method": "hf", "basis": "sto-3g"}
    inputs[field] = requested

    with pytest.raises(UnknownError, match=f"{field} mismatch"):
        cclib_harness._parse_and_convert(definition, execution, _atomic_input(**inputs))


@pytest.mark.parametrize(
    "protocol,expected",
    [("none", set()), ("input", {"input"}), ("all", {"input", "dispatch.out"})],
)
def test_native_files_protocols_and_complete_stdout(monkeypatch, protocol, expected):
    api, definition, execution, _ = _fake_conversion_case()
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)

    result = cclib_harness._parse_and_convert(
        definition, execution, _atomic_input(protocols={"native_files": protocol})
    )

    assert set(result.native_files) == expected
    assert "stdout" not in result.native_files
    assert "stderr" not in result.native_files
    assert "outfiles" not in result.extras
    assert result.stdout == execution.output_text
    if protocol != "none":
        assert result.native_files["input"] == execution.input_text
    if protocol == "all":
        assert result.native_files["dispatch.out"] == execution.output_text


def test_public_dispatch_returns_v1_while_direct_harness_preserves_v2_request(monkeypatch):
    api, definition, execution, _ = _fake_conversion_case()
    input_model = _atomic_input(protocols={"native_files": "input"})
    job = cclib_harness._Job(
        command=[execution.executable],
        infiles={execution.input_filename: execution.input_text},
        outfiles=[execution.output_filename],
        input_filename=execution.input_filename,
        output_filename=execution.output_filename,
        input_text=execution.input_text,
        executable=execution.executable,
    )
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)
    monkeypatch.setattr(cclib_harness, "which", lambda name: execution.executable)
    monkeypatch.setattr(cclib_harness, "_probe_executable", lambda *args, **kwargs: "6.2")
    monkeypatch.setattr(cclib_harness, "_execute_job", lambda actual_definition, actual_job, config: execution)
    monkeypatch.setitem(
        cclib_harness._PROGRAM_DEFINITIONS,
        "qchem",
        replace(
            definition,
            generator=lambda input_data, config, executable: job,
            preflight=lambda executable, environment: dict(environment),
        ),
    )

    direct = qcng.get_program("cclib-qchem", check=False).compute(input_model, _task_config())
    public = qcng.compute(
        input_model,
        "cclib-qchem",
        raise_error=True,
        task_config={"ncores": 1, "memory": 1.0},
        return_version=1,
        return_dict=False,
    )

    assert direct.schema_version == 2
    assert direct.input_data == input_model
    assert public.schema_version == 1
    assert public.driver.value == input_model.specification.driver.value
    assert public.model.method == input_model.specification.model.method
    assert public.model.basis == input_model.specification.model.basis


@pytest.mark.parametrize("protocol,retains_input", [("none", False), ("all", True)])
def test_native_files_v1_validation_incompatibility_uses_only_documented_input_fallback(
    monkeypatch, protocol, retains_input
):
    api, definition, execution, _ = _fake_conversion_case()
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)
    original_validate = cclib_harness._validate_v1_atomic_result

    def reject_native_files(output):
        if "native_files" in output:
            raise ValueError("native_files is not a permitted v1 field")
        return original_validate(output)

    monkeypatch.setattr(cclib_harness, "_validate_v1_atomic_result", reject_native_files)

    result = cclib_harness._parse_and_convert(
        definition, execution, _atomic_input(protocols={"native_files": protocol})
    )

    assert not result.native_files
    if retains_input:
        assert result.extras["cclib_harness"]["native_input"] == execution.input_text
    else:
        assert "native_input" not in result.extras["cclib_harness"]


@uusing("cclib-qchem")
def test_live_cclib_qchem_water_mp2_energy():
    result = qcng.compute(
        _water_input(),
        "cclib-qchem",
        raise_error=True,
        return_version=1,
    )

    assert result.success is True
    assert result.return_result == pytest.approx(-75.00228214, abs=1.0e-6)


_REAL_CCLIB_FIXTURES = [
    pytest.param("qchem", "QChem/basicQChem5.1/water_mp2.out", id="qchem-water-mp2"),
    pytest.param(
        "qchem",
        "QChem/basicQChem5.1/dvb_dispersion_bp86_d3zero.out",
        id="qchem-bp86-energy",
    ),
    pytest.param("qchem", "QChem/basicQChem5.1/water_ir.out", id="qchem-water-ir"),
    pytest.param("qchem", "QChem/basicQChem5.1/water_ccsd.out", id="qchem-water-ccsd"),
    pytest.param("orca", "ORCA/basicORCA6.0/water_mp2.out", id="orca-water-mp2"),
    pytest.param("orca", "ORCA/basicORCA6.0/dvb_sp_hf.out", id="orca-hf-energy"),
    pytest.param(
        "orca",
        "ORCA/basicORCA6.0/dvb_ir.out",
        id="orca-dvb-ir-deferred",
        marks=pytest.mark.skip(
            reason=(
                "deferred: cclib's successful ORCA DFT parse lacks metadata.functional, "
                "so QCSchemaWriter raises KeyError"
            )
        ),
    ),
    pytest.param("orca", "ORCA/basicORCA6.0/water_ccsd.out", id="orca-water-ccsd"),
]


@pytest.mark.parametrize("program,relative_path", _REAL_CCLIB_FIXTURES)
def test_real_cclib_fixture_parse_and_conversion(program, relative_path):
    source_root = os.environ.get("CCLIB_SOURCE_ROOT")
    if source_root is None:
        pytest.skip("CCLIB_SOURCE_ROOT is not set")
    data_root = os.path.join(source_root, "data")
    if not os.path.isdir(data_root):
        pytest.skip(f"cclib fixture data directory is absent: {data_root}")
    fixture = os.path.join(data_root, relative_path)
    if not os.path.isfile(fixture):
        pytest.skip(f"cclib fixture is absent: {fixture}")

    api = cclib_harness._load_cclib_api()
    initial_parser = api.ccopen(fixture)
    assert type(initial_parser) is cclib_harness._PROGRAM_DEFINITIONS[program].parser_type()
    try:
        parsed = initial_parser.parse()
    finally:
        initial_parser.inputfile.close()
    assert parsed.metadata["success"] is True
    writer_output = api.QCSchemaWriter(parsed).as_dict(validate=False)
    initial_v1 = cclib_harness._validate_v1_atomic_result(writer_output)

    input_model = AtomicInput(
        molecule=initial_v1.molecule.convert_v(2),
        specification={
            "driver": writer_output["driver"],
            "model": writer_output["model"],
            "keywords": {},
            "protocols": {"native_files": "none"},
        },
    )
    with open(fixture, encoding="utf-8", errors="replace") as handle:
        output_text = handle.read()
    definition = cclib_harness._PROGRAM_DEFINITIONS[program]
    execution = cclib_harness._ExecutionResult(
        process_success=True,
        executable=f"/fixture/{program}",
        input_filename=definition.input_filename,
        output_filename=definition.output_filename,
        input_text=f"fixture input for {relative_path}",
        output_text=output_text,
        stdout=output_text if program == "orca" else "",
        stderr="",
    )

    result = cclib_harness._parse_and_convert(definition, execution, input_model)

    assert result.success is True
    assert result.schema_version == 2
    assert result.input_data == input_model
    assert result.molecule == initial_v1.molecule.convert_v(2)
    assert result.stdout == output_text
    assert result.provenance.creator == initial_v1.provenance.creator
    assert result.provenance.version == initial_v1.provenance.version
    assert result.input_data.specification.driver.value == writer_output["driver"].lower()
    assert result.input_data.specification.model.method.lower() == writer_output["model"]["method"].lower()
    assert result.input_data.specification.model.basis.lower() == writer_output["model"]["basis"].lower()
    flat_extras = {key: value for key, value in result.extras.items() if key != "cclib_harness"}
    assert flat_extras
    assert result.extras["cclib_harness"] == {
        "selector": definition.selector,
        "cclib_version": api.version,
        "parser": definition.parser_name,
        "executable": execution.executable,
    }
