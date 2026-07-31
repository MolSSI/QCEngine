import inspect
import json
import os
import subprocess
import sys
from dataclasses import replace
from types import SimpleNamespace

import pytest

import qcengine as qcng
import qcengine.programs.cclib_programs.base as cclib_base
import qcengine.programs.cclib_programs.cclib_orca as cclib_orca
import qcengine.programs.cclib_programs.cclib_qchem as cclib_qchem
from qcelemental.models.v2 import AtomicInput, BasisSet
from qcelemental.util import parse_version

from qcengine.programs.cclib_programs import ORCACCLibHarness, QChemCCLibHarness

from qcengine.config import TaskConfig
from qcengine.exceptions import InputError, ResourceError, UnknownError
from qcengine.programs.cclib_programs.base import CCLibHarness
from qcengine.testing import using, uusing


def test_registration_and_independent_concrete_instances(monkeypatch):
    from qcengine import testing

    qchem = qcng.get_program("CCLIB-QCHEM", check=False)
    orca = qcng.get_program("cclib-orca", check=False)

    for selector in ("cclib-qchem", "cclib-orca"):
        assert selector in testing._programs
        monkeypatch.setitem(testing._programs, selector, True)
        testing._using_cache.pop(selector, None)

        def marked_test():
            pass

        marked = testing.uusing(selector)(marked_test)
        marks = {mark.name: mark for mark in marked.pytestmark}
        assert set(marks) == {"skipif", "addon", selector}
        assert testing.has_program(selector) is True
        testing._using_cache.pop(selector, None)

    assert type(qcng.get_program("cclib-qchem", check=False)) is QChemCCLibHarness
    assert type(qcng.get_program("cclib-orca", check=False)) is ORCACCLibHarness
    assert (qchem.name, qchem.node_parallel) == ("cclib-qchem", False)
    assert (orca.name, orca.node_parallel) == ("cclib-orca", True)
    assert qchem is not orca
    assert qcng.get_program("qchem", check=False) is not qchem
    assert {"qchem", "cclib-qchem", "cclib-orca"} <= qcng.list_all_programs()

    with pytest.raises(Exception, match="frozen"):
        qchem.name = "cclib-orca"


def test_import_qcengine_does_not_import_external_cclib():
    script = r"""
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
"""
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
    return "                         O   R   C   A\n" f"Program Version {version}.0\n" "ORCA TERMINATED NORMALLY\n"


def test_missing_executable_is_a_resource_error(monkeypatch):
    harness = ORCACCLibHarness()
    monkeypatch.setattr(cclib_base, "_load_cclib_api", lambda: object())
    monkeypatch.setattr(cclib_base, "which", lambda command: None, raising=False)

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
        ({}, []),
        ({"QC": None}, ["QC"]),
        ({"QCAUX": None}, ["QCAUX"]),
        ({"QCPROG": None}, ["QCPROG"]),
        ({"QC": "file", "QCAUX": "file", "QCPROG": "directory"}, ["QC", "QCAUX", "QCPROG"]),
        ({"permissions": True}, ["QC", "QCAUX", "QCPROG", "resolved qchem executable"]),
    ],
)
def test_qchem_environment_preflight(tmp_path, changes, invalid):
    environment, executable = _valid_qchem_environment(tmp_path)
    ordinary_file = tmp_path / "ordinary"
    ordinary_file.write_text("data")
    ordinary_dir = tmp_path / "ordinary-dir"
    ordinary_dir.mkdir()
    replacements = {"file": str(ordinary_file), "directory": str(ordinary_dir)}
    if changes.pop("permissions", False):
        os.chmod(environment["QC"], 0o300)
        os.chmod(environment["QCAUX"], 0o300)
        os.chmod(environment["QCPROG"], 0o644)
        os.chmod(executable, 0o644)
    else:
        for variable, value in changes.items():
            if value is None:
                environment.pop(variable)
            else:
                environment[variable] = replacements[value]

    try:
        if invalid:
            with pytest.raises(ResourceError) as exc_info:
                cclib_qchem.preflight(executable, environment)
            message = str(exc_info.value)
            assert all(variable in message for variable in invalid)
            assert "cclib-qchem" in message
        else:
            environment["UNCHANGED"] = "preserved"
            child = cclib_qchem.preflight(executable, environment)
            assert child is not environment
            assert child["UNCHANGED"] == "preserved"
            assert child["QCSCRATCH"]
            assert "QCSCRATCH" not in environment
    finally:
        if invalid and "permissions" not in changes:
            for name in ("QC", "QCAUX"):
                path = environment.get(name)
                if path and os.path.isdir(path):
                    os.chmod(path, 0o700)


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


@pytest.mark.parametrize(
    "input_model,match",
    [
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
def test_input_fields_reject_invalid_common_structural_requests(input_model, match):
    with pytest.raises(InputError, match=match):
        cclib_base._input_fields(input_model)


@pytest.mark.parametrize(
    "builder,field,value",
    [
        pytest.param(cclib_qchem.build_input, "method", "", id="qchem-blank-method"),
        pytest.param(cclib_qchem.build_input, "basis", "   ", id="qchem-blank-basis"),
        pytest.param(
            cclib_qchem.build_input,
            "method",
            "hf\n$end\n$molecule\n9 1\nH 0 0 0",
            id="qchem-coordinate-section",
        ),
        pytest.param(
            cclib_qchem.build_input,
            "basis",
            "sto-3g\r\nMEM_TOTAL 999999",
            id="qchem-resource-override",
        ),
        pytest.param(cclib_qchem.build_input, "method", "hf\n@@@", id="qchem-multi-job"),
        pytest.param(cclib_qchem.build_input, "method", "hf\x00evil", id="qchem-control"),
        pytest.param(
            cclib_orca.build_input,
            "method",
            "hf\n%pal\nnprocs 999\nend",
            id="orca-resource-override",
        ),
        pytest.param(
            cclib_orca.build_input,
            "basis",
            "sto-3g\r\n* xyz 9 1\nH 0 0 0\n*",
            id="orca-coordinate-section",
        ),
        pytest.param(cclib_orca.build_input, "method", "hf\n$new_job", id="orca-multi-job"),
        pytest.param(cclib_orca.build_input, "basis", "sto-3g\x1bevil", id="orca-control"),
        pytest.param(cclib_orca.build_input, "basis", "sto-3g\t%pal", id="orca-whitespace"),
    ],
)
def test_native_method_and_basis_tokens_reject_structure_escape(builder, field, value):
    values = {"method": "hf", "basis": "sto-3g"}
    values[field] = value

    with pytest.raises(InputError, match=field):
        builder(_atomic_input(**values), _task_config(), "/resolved/program")


def test_native_method_and_basis_tokens_normalize_edges_and_preserve_punctuation():
    qchem = cclib_qchem.build_input(
        _atomic_input(method="  wb97x-d3(0)  ", basis="  6-31+g(d,p)  "),
        _task_config(),
        "/opt/qchem",
    )
    assert "METHOD wb97x-d3(0)\n" in qchem.input_text
    assert "BASIS 6-31+g(d,p)\n" in qchem.input_text

    orca = cclib_orca.build_input(
        _atomic_input(method="  dlpno-ccsd(t)  ", basis="  def2-svp/c  "),
        _task_config(),
        "/opt/orca",
    )
    assert orca.input_text.startswith("! dlpno-ccsd(t) def2-svp/c\n")


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


def test_qchem_exact_input_and_driver_mapping():
    job = cclib_qchem.build_input(_water_input(), _task_config(), "/opt/qchem")

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

    for driver, jobtype in (("energy", "sp"), ("gradient", "force"), ("hessian", "freq")):
        mapped = cclib_qchem.build_input(_atomic_input(driver=driver, method="wb97x-d"), _task_config(), "/opt/qchem")
        assert f"JOBTYPE {jobtype}\n" in mapped.input_text
        assert "METHOD wb97x-d\n" in mapped.input_text
    with pytest.raises(InputError, match="cclib-qchem.*driver"):
        cclib_qchem.build_input(_atomic_input(driver="properties"), _task_config(), "/opt/qchem")


def test_qchem_scalar_keyword_rendering_and_malformed_reserved_rejection():
    keywords = {
        "zeta": "value",
        "a_bool_true": True,
        "bool_false": False,
        "an_int": 7,
        "a_float": 1.25,
    }
    job = cclib_qchem.build_input(_atomic_input(keywords=keywords), _task_config(), "/opt/qchem")

    ordinary = [
        "AN_INT 7",
        "A_BOOL_TRUE TRUE",
        "A_FLOAT 1.25",
        "BOOL_FALSE FALSE",
        "ZETA value",
    ]
    positions = [job.input_text.index(line) for line in ordinary]
    assert positions == sorted(positions)

    malformed = [
        {"bad": "line one\nline two"},
        {"bad\nkey": "value"},
        {"bad": None},
        {"bad": [1]},
        {"bad": {"nested": 1}},
        {"bad": float("nan")},
        {"": "value"},
        {"ordinary key": "value"},
        {"ordinary\x00key": "value"},
        {"$end": "value"},
        {"$molecule": "value"},
        {"-leading": "value"},
        {"nonascii_é": "value"},
    ]
    for candidate in malformed:
        with pytest.raises(InputError, match="keyword"):
            cclib_qchem.build_input(_atomic_input(keywords=candidate), _task_config(), "/opt/qchem")
    for reserved in _QCHEM_RESERVED:
        with pytest.raises(InputError, match="reserved"):
            cclib_qchem.build_input(_atomic_input(keywords={reserved.swapcase(): "user"}), _task_config(), "/opt/qchem")
        with pytest.raises(InputError, match="keyword name"):
            cclib_qchem.build_input(_atomic_input(keywords={f" {reserved}": "user"}), _task_config(), "/opt/qchem")
    with pytest.raises(InputError, match="collision"):
        cclib_qchem.build_input(_atomic_input(keywords={"thresh": 8, "THRESH": 10}), _task_config(), "/opt/qchem")


def test_orca_exact_input_and_energy_only_policy():
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
    job = cclib_orca.build_input(input_model, _task_config(), "/opt/orca")

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
    assert cclib_orca.build_input(
        _atomic_input(method="dlpno-ccsd(t)"), _task_config(), "/opt/orca"
    ).input_text.startswith("! dlpno-ccsd(t) sto-3g")
    for driver in ("gradient", "hessian", "properties"):
        with pytest.raises(InputError, match="cclib-orca.*only.*energy"):
            cclib_orca.build_input(_atomic_input(driver=driver), _task_config(), "/opt/orca")


def test_orca_deterministic_blocks_and_malformed_reserved_rejection():
    keywords = {
        "simple": ["rks", "usesym", "TightSCF"],
        "blocks": {
            "scf": "MaxIter 200",
            "output": "PrintLevel Mini\nPrint[P_AtCharges_M] 1",
            "basis": 'NewGTO H "def2-TZVP" end',
        },
    }
    job = cclib_orca.build_input(_atomic_input(driver="energy", keywords=keywords), _task_config(), "/opt/orca")

    assert job.input_text.splitlines()[0] == "! hf sto-3g rks usesym TightSCF"
    defaults_end = job.input_text.index("Print[P_Hirshfeld] 1")
    user_output = job.input_text.index("PrintLevel Mini")
    basis = job.input_text.index("%basis")
    scf = job.input_text.index("%scf")
    resources = job.input_text.index("%pal")
    assert defaults_end < user_output < basis < scf < resources

    harmless = [
        "MaxIter 200 # weekend schedule",
        "Print[ P_Basis ] 2 # percentage %pal is only text here",
        "SomeValue prefix_end suffix",
        "Label job$new_job_backup",
    ]
    for body in harmless:
        rendered = cclib_orca.build_input(
            _atomic_input(keywords={"blocks": {"scf": body}}), _task_config(), "/opt/orca"
        )
        assert body in rendered.input_text

    outer_syntax = [
        "end\n* xyz 9 1\nH 0 0 0\n*",
        "  EnD trailing text  ",
        "%pal\nnprocs 99\nend",
        "  %MAXCORE 9999 # comment",
        "%coords\nctyp xyz",
        "$new_job",
        " * xyz 0 1",
    ]
    for body in outer_syntax:
        with pytest.raises(InputError, match="block body|outer syntax"):
            cclib_orca.build_input(_atomic_input(keywords={"blocks": {"output": body}}), _task_config(), "/opt/orca")

    malformed = [
        {"unknown": []},
        {"simple": "rks"},
        {"simple": [""]},
        {"simple": ["rks\n* xyz 9 9"]},
        {"blocks": "output"},
        {"blocks": {"bad-name": "value"}},
        {"blocks": {"scf": None}},
        {"coordinates": []},
    ]
    for candidate in malformed:
        with pytest.raises(InputError, match="ORCA|unknown|coordinate"):
            cclib_orca.build_input(_atomic_input(keywords=candidate), _task_config(), "/opt/orca")
    for block in ("pal", "maxcore", "coords"):
        with pytest.raises(InputError, match="reserved|coordinate"):
            cclib_orca.build_input(
                _atomic_input(keywords={"blocks": {block: "user body"}}), _task_config(), "/opt/orca"
            )


def test_primary_output_selection(monkeypatch, tmp_path):
    job = cclib_qchem.build_input(_atomic_input(), _task_config(), "/resolved/qchem")
    definition = replace(cclib_qchem.QCHEM_DEFINITION, preflight=lambda exe, env: dict(env))
    inherited_environment = os.environ.copy()
    calls = []

    def fake_execute(command, infiles=None, outfiles=None, **kwargs):
        calls.append((command, infiles, outfiles, kwargs))
        return True, {
            "stdout": "launcher stdout",
            "stderr": "launcher stderr",
            "outfiles": {"dispatch.out": _qchem_probe_output()},
        }

    monkeypatch.setattr(cclib_base, "execute", fake_execute)
    result = cclib_base._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    command, infiles, outfiles, kwargs = calls[0]
    assert command == job.command
    assert infiles == job.infiles
    assert outfiles == ["dispatch.out"]
    assert {key: kwargs["environment"][key] for key in inherited_environment if key != "QCSCRATCH"} == {
        key: value for key, value in inherited_environment.items() if key != "QCSCRATCH"
    }
    assert result.output_text == _qchem_probe_output()
    assert (result.stdout, result.stderr) == ("launcher stdout", "launcher stderr")

    orca_job = cclib_orca.build_input(_atomic_input(), _task_config(), "/resolved/orca")
    orca_output = _orca_probe_output()
    monkeypatch.setattr(
        cclib_base,
        "execute",
        lambda *args, **kwargs: (True, {"stdout": orca_output, "stderr": "", "outfiles": {}}),
    )
    orca_result = cclib_base._execute_job(
        cclib_orca.ORCA_DEFINITION, orca_job, _task_config(scratch_directory=str(tmp_path))
    )
    assert orca_result.output_text == orca_output
    assert orca_result.stdout == orca_output


def test_managed_scratch_and_qcscratch_isolation(monkeypatch, tmp_path):
    monkeypatch.setenv("QCSCRATCH", "/inherited/unmanaged")
    job = cclib_qchem.build_input(_atomic_input(), _task_config(), "/resolved/qchem")
    definition = replace(cclib_qchem.QCHEM_DEFINITION, preflight=lambda exe, env: dict(env))
    observed = {}

    def fake_execute(command, infiles=None, outfiles=None, **kwargs):
        qcscratch = kwargs["environment"]["QCSCRATCH"]
        observed["qcscratch"] = qcscratch
        assert os.path.isdir(qcscratch)
        assert os.path.commonpath([qcscratch, str(tmp_path)]) == str(tmp_path)
        assert kwargs["scratch_directory"] == qcscratch
        return True, {"stdout": "", "stderr": "", "outfiles": {"dispatch.out": _qchem_probe_output()}}

    monkeypatch.setattr(cclib_base, "execute", fake_execute)
    cclib_base._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    assert not os.path.exists(observed["qcscratch"])
    source = inspect.getsource(cclib_base).casefold()
    assert "qcscratch" not in source
    assert "qchem" not in source
    assert "orca" not in source


def _qchem_execution_case():
    job = cclib_qchem.build_input(_atomic_input(), _task_config(), "/resolved/qchem")
    definition = replace(cclib_qchem.QCHEM_DEFINITION, preflight=lambda exe, env: dict(env))
    return definition, job


def _assert_execution_message(error, definition, job, stage, diagnostic):
    message = str(error)
    assert definition.selector in message
    assert job.executable in message
    assert stage in message
    if len(diagnostic) <= 4000 and len(diagnostic.splitlines()) <= 40:
        assert diagnostic in message


def test_diagnostic_tail_has_independent_line_and_character_boundaries():
    definition, job = _qchem_execution_case()
    line_diagnostic = "\n".join(f"boundary-line-{index:02d}" for index in range(42))

    with pytest.raises(UnknownError) as line_error:
        cclib_base._raise_execution_failure(definition, job, "execution", line_diagnostic)
    line_tail = str(line_error.value).split("Diagnostic tail:\n", 1)[1]
    assert line_tail.startswith("boundary-line-02\n")
    assert "boundary-line-01" not in line_tail
    assert line_tail.endswith("boundary-line-41")
    assert line_tail.count("\n") == 39

    character_diagnostic = "omitted-prefix:" + ("x" * 4000)
    with pytest.raises(UnknownError) as character_error:
        cclib_base._raise_execution_failure(definition, job, "execution", character_diagnostic)
    character_tail = str(character_error.value).split("Diagnostic tail:\n", 1)[1]
    assert character_tail == "x" * 4000
    assert "omitted-prefix" not in character_tail


@pytest.mark.parametrize(
    "case",
    ["nonzero", "missing-output", "termination", "qchem-environment", "license", "exception"],
)
def test_bounded_execution_and_resource_failure_classification(monkeypatch, tmp_path, case):
    definition, job = _qchem_execution_case()
    diagnostic = "\n".join(f"failure line {index}: {'x' * 120}" for index in range(60))
    expected_error = UnknownError
    stage = "execution"
    original = None

    if case == "nonzero":
        process = (False, {"stdout": "", "stderr": "", "outfiles": {"dispatch.out": diagnostic}})
    elif case == "missing-output":
        process = (True, {"stdout": diagnostic, "stderr": "", "outfiles": {"dispatch.out": None}})
        stage = "output selection"
    elif case == "termination":
        diagnostic = "Q-Chem stopped before its farewell"
        process = (True, {"stdout": "", "stderr": "", "outfiles": {"dispatch.out": diagnostic}})
        stage = "termination"
    elif case == "qchem-environment":
        diagnostic = "Undefined environment variable QCAUX"
        process = (False, {"stdout": "", "stderr": diagnostic, "outfiles": {"dispatch.out": None}})
        expected_error = ResourceError
    elif case == "license":
        definition = cclib_orca.ORCA_DEFINITION
        job = cclib_orca.build_input(_atomic_input(), _task_config(), "/resolved/orca")
        diagnostic = "license checkout failed"
        process = (False, {"stdout": diagnostic, "stderr": "", "outfiles": {}})
        expected_error = ResourceError
    else:
        definition = cclib_orca.ORCA_DEFINITION
        job = cclib_orca.build_input(_atomic_input(), _task_config(), "/resolved/orca")
        original = OSError("scheduler launch failed")
        diagnostic = str(original)
        process = None

    def execute(*args, **kwargs):
        if original is not None:
            raise original
        return process

    monkeypatch.setattr(cclib_base, "execute", execute)
    with pytest.raises(expected_error) as exc_info:
        cclib_base._execute_job(definition, job, _task_config(scratch_directory=str(tmp_path)))

    _assert_execution_message(exc_info.value, definition, job, stage, diagnostic)
    if diagnostic.startswith("failure line"):
        diagnostic_tail = str(exc_info.value).split("Diagnostic tail:\n", 1)[1]
        assert len(diagnostic_tail) == 4000
        assert "failure line 29" not in diagnostic_tail
        assert diagnostic_tail.endswith("failure line 59: " + ("x" * 120))
    if case in {"missing-output", "exception"}:
        assert exc_info.value.__cause__ is not None


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

    api = cclib_base.CCLibAPI(
        version="1.9.fake",
        QCSchemaWriter=FakeWriter,
        ccopen=ccopen,
    )
    execution = cclib_base.ExecutionResult(
        process_success=True,
        executable=f"/resolved/{program}",
        input_filename="dispatch.in" if program == "qchem" else "dispatch.inp",
        output_filename="dispatch.out",
        input_text="complete native input",
        output_text="complete parsed output",
        stdout="launcher stdout",
        stderr="launcher stderr",
    )
    definitions = {"qchem": cclib_qchem.QCHEM_DEFINITION, "orca": cclib_orca.ORCA_DEFINITION}
    definition = replace(definitions[program], parser_type=lambda: parser_class)
    return api, definition, execution, opened


def _assert_conversion_failure(monkeypatch, match, **case):
    api, definition, execution, _ = _fake_conversion_case(**case)
    monkeypatch.setattr(cclib_base, "_load_cclib_api", lambda: api)
    with pytest.raises(UnknownError, match=match) as exc_info:
        cclib_base._parse_and_convert(definition, execution, _atomic_input())
    assert definition.selector in str(exc_info.value)
    assert exc_info.value.__cause__ is not None
    return exc_info.value


def _track_parser_temporary_file(monkeypatch):
    real_named_temporary_file = cclib_base.tempfile.NamedTemporaryFile
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

    monkeypatch.setattr(cclib_base.tempfile, "NamedTemporaryFile", tracked_named_temporary_file)
    return state


@pytest.mark.parametrize("case", ["success", "identity", "parse", "writer"])
def test_parser_temporary_file_cleanup_on_success_and_failure(monkeypatch, case):
    options = {
        "identity": {"mismatch": True},
        "parse": {"parser_failure": ValueError("parser exploded")},
        "writer": {"writer_output": RuntimeError("writer exploded")},
    }.get(case, {})
    api, definition, execution, _ = _fake_conversion_case(**options)
    state = _track_parser_temporary_file(monkeypatch)
    original_ccopen = api.ccopen

    def assert_closed_before_reopen(source):
        assert state["temporary"].closed
        assert source.name == state["path"]
        return original_ccopen(source)

    monkeypatch.setattr(cclib_base, "_load_cclib_api", lambda: replace(api, ccopen=assert_closed_before_reopen))
    if case == "success":
        assert cclib_base._parse_and_convert(definition, execution, _atomic_input()).success is True
    else:
        with pytest.raises(UnknownError):
            cclib_base._parse_and_convert(definition, execution, _atomic_input())

    assert state["kwargs"] == {"mode": "w", "suffix": ".out", "encoding": "utf-8", "delete": False}
    assert state["temporary"].closed
    assert not os.path.exists(state["path"])


@pytest.mark.parametrize("stage", ["auto-detection", "type-loading", "identity", "parse", "result-validation"])
def test_parser_detection_identity_and_parse_failure_classification(monkeypatch, stage):
    options = {}
    if stage == "auto-detection":
        options["detected"] = False
    elif stage == "identity":
        options["mismatch"] = True
    elif stage == "parse":
        options["parser_failure"] = ValueError("parser exploded")
    elif stage == "result-validation":
        options["parsed"] = SimpleNamespace(metadata={"success": False})
    api, definition, execution, _ = _fake_conversion_case(**options)
    if stage == "type-loading":
        original = ImportError("cclib parser module is unavailable")
        definition = replace(definition, parser_type=lambda: (_ for _ in ()).throw(original))
    execution = replace(execution, output_text="\n".join(f"line {index}" for index in range(1000)))
    monkeypatch.setattr(cclib_base, "_load_cclib_api", lambda: api)
    expected = {
        "auto-detection": "auto-detection",
        "type-loading": "type loading",
        "identity": "identity",
        "parse": "parse",
        "result-validation": "result validation",
    }[stage]

    with pytest.raises(UnknownError, match=f"parser {expected}") as exc_info:
        cclib_base._parse_and_convert(definition, execution, _atomic_input())

    assert "line 0" not in str(exc_info.value)
    assert "line 999" in str(exc_info.value)
    assert exc_info.value.__cause__ is not None


@pytest.mark.parametrize("case", ["required-field", "collision", "extras-type", "validation", "writer"])
def test_writer_required_fields_extras_collision_and_validation_failure(monkeypatch, case):
    output = _fake_writer_output()
    match = "writer output"
    writer_output = output
    if case == "required-field":
        output.pop("properties")
    elif case == "collision":
        output["extras"]["cclib_harness"] = {"writer_owned": "must survive"}
        match = "writer augmentation"
    elif case == "extras-type":
        output["extras"] = 7
        match = "writer augmentation"
    elif case == "validation":
        output["return_result"] = "not-an-energy"
        match = "QCSchema v1 validation"
    else:
        writer_output = RuntimeError("writer exploded")
        match = "QCSchema writer"

    api, definition, execution, _ = _fake_conversion_case(writer_output)
    monkeypatch.setattr(cclib_base, "_load_cclib_api", lambda: api)
    with pytest.raises(UnknownError, match=match) as exc_info:
        cclib_base._parse_and_convert(definition, execution, _atomic_input())

    assert definition.selector in str(exc_info.value)
    assert exc_info.value.__cause__ is not None
    if case == "collision":
        assert output["extras"]["cclib_harness"] == {"writer_owned": "must survive"}
        assert isinstance(exc_info.value.__cause__, ValueError)


@pytest.mark.parametrize(
    "program,requested,parsed,requested_basis,parsed_basis,creator,version",
    [
        pytest.param(
            "qchem",
            "wb97x-d",
            "wB97X-D3",
            "  sTo-3G  ",
            "\tSTo-3g ",
            "Q-Chem",
            "6.2",
            id="normalized_basis-qchem",
        ),
        pytest.param(
            "orca",
            "dlpno-ccsd(t)",
            "DLPNO-CCSD(T0)",
            "sto-3g",
            "STO-3G",
            "ORCA",
            "6.0.1",
            id="orca",
        ),
    ],
)
def test_successful_v1_to_v2_conversion_and_provenance_preservation(
    monkeypatch, program, requested, parsed, requested_basis, parsed_basis, creator, version
):
    output = _fake_writer_output(driver="energy", method=parsed, basis=parsed_basis)
    output["provenance"] = {"creator": creator, "version": version, "routine": "cclib.QCSchemaWriter"}
    if program == "orca":
        output["extras"]["dispersionenergies"] = [-0.001]
    api, definition, execution, opened = _fake_conversion_case(output, program=program)
    monkeypatch.setattr(cclib_base, "_load_cclib_api", lambda: api)
    original_validate = cclib_base._validate_v1_atomic_result
    validated_output = {}

    def capture_validated_output(value):
        validated_output.update(value)
        return original_validate(value)

    monkeypatch.setattr(cclib_base, "_validate_v1_atomic_result", capture_validated_output)
    input_model = _atomic_input(method=requested, basis=requested_basis)

    result = cclib_base._parse_and_convert(definition, execution, input_model)

    assert result.schema_version == 2
    assert result.input_data == input_model
    assert list(result.molecule.symbols) == ["H", "H"]
    assert result.stdout == execution.output_text
    assert result.provenance.creator == creator
    assert result.provenance.version == version
    assert validated_output["model"]["basis"] == parsed_basis
    assert result.input_data.specification.model.basis == requested_basis
    assert opened == [execution.output_text]
    assert set(result.extras) == set(output["extras"]) | {"cclib_harness"}
    for key, value in output["extras"].items():
        assert result.extras[key] == value
    assert result.extras["cclib_harness"] == {
        "selector": definition.selector,
        "cclib_version": "1.9.fake",
        "parser": definition.parser_name,
        "executable": execution.executable,
    }


@pytest.mark.parametrize(
    "field,requested,parsed",
    [
        ("driver", "energy", "gradient"),
        ("basis", "sto-3g", "6-31g"),
    ],
)
def test_parsed_driver_or_basis_mismatch_is_rejected_without_relabeling(monkeypatch, field, requested, parsed):
    values = {"driver": "energy", "method": "hf", "basis": "sto-3g"}
    values[field] = parsed
    output = _fake_writer_output(**values)
    api, definition, execution, _ = _fake_conversion_case(output)
    monkeypatch.setattr(cclib_base, "_load_cclib_api", lambda: api)
    inputs = {"driver": "energy", "method": "hf", "basis": "sto-3g"}
    inputs[field] = requested

    with pytest.raises(UnknownError, match=f"{field} mismatch"):
        cclib_base._parse_and_convert(definition, execution, _atomic_input(**inputs))


@pytest.mark.parametrize(
    "protocol,expected,v1_incompatible",
    [
        pytest.param("none", set(), False, id="none"),
        pytest.param("input", {"input"}, False, id="input"),
        pytest.param("all", {"input", "dispatch.out"}, False, id="all"),
        pytest.param("none", set(), True, id="v1_incompatible-none"),
        pytest.param("all", set(), True, id="v1_incompatible-all"),
    ],
)
def test_native_file_protocols_and_public_schema_conversion(monkeypatch, protocol, expected, v1_incompatible):
    api, definition, execution, _ = _fake_conversion_case()
    monkeypatch.setattr(cclib_base, "_load_cclib_api", lambda: api)
    if v1_incompatible:
        original_validate = cclib_base._validate_v1_atomic_result

        def reject_native_files(output):
            if "native_files" in output:
                raise ValueError("native_files is not a permitted v1 field")
            return original_validate(output)

        monkeypatch.setattr(cclib_base, "_validate_v1_atomic_result", reject_native_files)

    result = cclib_base._parse_and_convert(definition, execution, _atomic_input(protocols={"native_files": protocol}))

    native_files = result.native_files or {}
    assert set(native_files) == expected
    assert "stdout" not in native_files
    assert "stderr" not in native_files
    assert "outfiles" not in result.extras
    assert result.stdout == execution.output_text
    if protocol != "none" and not v1_incompatible:
        assert native_files["input"] == execution.input_text
    if protocol == "all" and not v1_incompatible:
        assert native_files["dispatch.out"] == execution.output_text
    harness_extras = result.extras["cclib_harness"]
    retains_fallback_input = v1_incompatible and protocol == "all"
    assert ("native_input" in harness_extras) is retains_fallback_input
    if retains_fallback_input:
        assert harness_extras["native_input"] == execution.input_text

    input_model = _atomic_input(protocols={"native_files": protocol})
    job = cclib_base.Job(
        command=[execution.executable],
        infiles={execution.input_filename: execution.input_text},
        outfiles=[execution.output_filename],
        input_filename=execution.input_filename,
        output_filename=execution.output_filename,
        input_text=execution.input_text,
        executable=execution.executable,
    )
    monkeypatch.setattr(cclib_base, "which", lambda name: execution.executable)
    monkeypatch.setattr(cclib_base, "_probe_executable", lambda *args, **kwargs: "6.2")
    monkeypatch.setattr(cclib_base, "_execute_job", lambda actual_definition, actual_job, config: execution)
    monkeypatch.setattr(
        QChemCCLibHarness,
        "definition",
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
    assert set(direct.native_files or {}) == expected
    direct_extras = direct.extras["cclib_harness"]
    assert ("native_input" in direct_extras) is retains_fallback_input
    if retains_fallback_input:
        assert direct_extras["native_input"] == execution.input_text
    assert public.schema_version == 1
    assert public.driver.value == input_model.specification.driver.value
    assert public.model.method == input_model.specification.model.method
    assert public.model.basis == input_model.specification.model.basis


@pytest.mark.parametrize(
    "selector,minimum",
    [
        pytest.param("cclib-qchem", "5.1", marks=[*using("cclib-qchem"), pytest.mark.cclib_qchem]),
        pytest.param("cclib-orca", "6.0", marks=[*using("cclib-orca"), pytest.mark.cclib_orca]),
    ],
)
def test_live_version_uses_available_software(selector, minimum):
    assert parse_version(qcng.get_program(selector).get_version()) >= parse_version(minimum)


@pytest.mark.cclib_qchem
@uusing("cclib-qchem")
def test_live_qchem_hessian():
    result = qcng.compute(
        _atomic_input(
            driver="hessian",
            method="hf",
            basis="sto-3g",
            molecule={
                "symbols": ["H", "H"],
                "geometry": [0.0, 0.0, -0.7, 0.0, 0.0, 0.7],
                "fix_com": True,
                "fix_orientation": True,
            },
        ),
        "cclib-qchem",
        raise_error=True,
        task_config={"ncores": 1, "memory": 1.0},
    )
    assert result.success is True
    assert result.return_result.shape == (6, 6)


@pytest.mark.parametrize(
    "selector,input_model,expected,absolute_tolerance",
    [
        pytest.param(
            "cclib-qchem",
            _water_input(),
            -75.00228214,
            1.0e-6,
            marks=[*using("cclib-qchem"), pytest.mark.cclib_qchem],
            id="qchem-water-mp2",
        ),
        # ORCA 6.x HF/STO-3G for He; 1e-8 Eh permits only final-digit output/parser variation.
        pytest.param(
            "cclib-orca",
            _atomic_input(),
            -2.8077839575,
            1.0e-8,
            marks=[*using("cclib-orca"), pytest.mark.cclib_orca],
            id="orca-he-hf",
        ),
    ],
)
def test_live_qchem_and_orca_energy_calculations(selector, input_model, expected, absolute_tolerance):
    result = qcng.compute(
        input_model,
        selector,
        raise_error=True,
        return_version=1,
        task_config={"ncores": 1, "memory": 1.0},
    )

    assert result.success is True
    assert isinstance(result.return_result, float)
    assert result.return_result == pytest.approx(expected, abs=absolute_tolerance)


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

    api = cclib_base._load_cclib_api()
    initial_parser = api.ccopen(fixture)
    definitions = {"qchem": cclib_qchem.QCHEM_DEFINITION, "orca": cclib_orca.ORCA_DEFINITION}
    assert type(initial_parser) is definitions[program].parser_type()
    try:
        parsed = initial_parser.parse()
    finally:
        initial_parser.inputfile.close()
    assert parsed.metadata["success"] is True
    writer_output = api.QCSchemaWriter(parsed).as_dict(validate=False)
    initial_v1 = cclib_base._validate_v1_atomic_result(writer_output)

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
    definition = {"qchem": cclib_qchem.QCHEM_DEFINITION, "orca": cclib_orca.ORCA_DEFINITION}[program]
    execution = cclib_base.ExecutionResult(
        process_success=True,
        executable=f"/fixture/{program}",
        input_filename=definition.input_filename,
        output_filename=definition.output_filename,
        input_text=f"fixture input for {relative_path}",
        output_text=output_text,
        stdout=output_text if program == "orca" else "",
        stderr="",
    )

    result = cclib_base._parse_and_convert(definition, execution, input_model)

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
