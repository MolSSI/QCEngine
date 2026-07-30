import copy
import importlib
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


@pytest.fixture
def fake_cclib_api():
    class FakeCCData:
        def __init__(self, attributes=None):
            self.attributes = attributes or {}

    class FakeQCSchemaWriter:
        calls = 0
        mutation = None

        def __init__(self, data):
            self.data = data

        def as_dict(self, validate=True):
            type(self).calls += 1
            assert validate is False
            assert self.data.attributes["atomcoords"] == [[[1.0, 0.0, 0.0]]]
            assert self.data.attributes["atomcharges"] == {"mulliken": [0.0]}
            output = {
                "schema_name": "qcschema_output",
                "schema_version": 1,
                "molecule": {
                    "geometry": [1.8897261255, 0.0, 0.0],
                    "molecular_charge": 0,
                    "molecular_multiplicity": 1,
                    "schema_name": "qcschema_molecule",
                    "schema_version": 2,
                    "symbols": ["He"],
                    "validated": True,
                },
                "provenance": {
                    "creator": "Synthetic",
                    "version": "1.0",
                    "routine": "fake.QCSchemaWriter",
                },
                "success": True,
                "error": None,
                "stdout": None,
                "stderr": None,
                "extras": {
                    "atomcharges": {"mulliken": [0.0]},
                    "atomcoords": [[[1.8897261255, 0.0, 0.0]]],
                    "atomnos": [2],
                    "charge": 0,
                    "homos": [0],
                    "mult": 1,
                    "natom": 1,
                    "nbasis": 1,
                    "nmo": 1,
                    "scfenergies": [-0.9040333652549597],
                    "scftargets": [[[1.0e-6]]],
                    "scfvalues": [[[1.0e-4], [1.0e-7]]],
                },
                "driver": "energy",
                "keywords": {},
                "model": {"method": "hf", "basis": "sto-3g"},
                "properties": {
                    "calcinfo_nalpha": 1,
                    "calcinfo_natom": 1,
                    "calcinfo_nbasis": 1,
                    "calcinfo_nbeta": 1,
                    "calcinfo_nmo": 1,
                    "return_energy": -0.9040333652549597,
                    "scf_iterations": 2,
                    "scf_total_energy": -0.9040333652549597,
                },
                "return_result": -0.9040333652549597,
            }
            if type(self).mutation == "wrong_units":
                output["molecule"]["geometry"][0] = 1.0
                output["extras"]["atomcoords"][0][0][0] = 1.0
            elif type(self).mutation == "missing_atomcoords":
                output["extras"].pop("atomcoords")
            elif type(self).mutation == "missing_atomcharges":
                output["extras"].pop("atomcharges")
            return output

    class FakeQChem:
        pass

    class FakeORCA:
        pass

    return cclib_harness._CCLibAPI(
        version="1.9.test",
        ccData=FakeCCData,
        QCSchemaWriter=FakeQCSchemaWriter,
        ccread=lambda source: source,
        QChem=FakeQChem,
        ORCA=FakeORCA,
    )


@pytest.fixture(autouse=True)
def clear_cclib_caches():
    cclib_harness.CCLibHarness.version_cache.clear()
    if hasattr(cclib_harness, "_cclib_compatibility_cache"):
        cclib_harness._cclib_compatibility_cache = None
    yield
    cclib_harness.CCLibHarness.version_cache.clear()
    if hasattr(cclib_harness, "_cclib_compatibility_cache"):
        cclib_harness._cclib_compatibility_cache = None


def test_cclib_missing_is_reported_by_compatibility_probe(monkeypatch):
    def missing():
        raise ModuleNotFoundError("No module named 'cclib'")

    monkeypatch.setattr(cclib_harness, "_load_cclib_api", missing)
    success, message = cclib_harness._check_cclib_compatibility()

    assert success is False
    assert "cclib" in message
    assert "not importable" in message


def test_cclib_compatibility_validates_v1_geometry_and_flat_extras(monkeypatch, fake_cclib_api):
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: fake_cclib_api)

    success, version = cclib_harness._check_cclib_compatibility()
    assert success is True
    assert version == "1.9.test"

    output = fake_cclib_api.QCSchemaWriter(
        cclib_harness._synthetic_ccdata(fake_cclib_api.ccData)
    ).as_dict(validate=False)
    result = cclib_harness._validate_v1_atomic_result(output)
    assert result.schema_version == 1
    assert result.molecule.geometry[0][0] == pytest.approx(1.8897261255, abs=2.0e-9)
    assert result.extras["atomcoords"][0][0][0] == pytest.approx(1.8897261255, abs=2.0e-9)
    assert result.extras["atomcharges"] == {"mulliken": [0.0]}


@pytest.mark.parametrize("mutation", ["wrong_units", "missing_atomcoords", "missing_atomcharges"])
def test_cclib_compatibility_rejects_incomplete_writer(monkeypatch, fake_cclib_api, mutation):
    fake_cclib_api.QCSchemaWriter.mutation = mutation
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: fake_cclib_api)

    success, message = cclib_harness._check_cclib_compatibility()

    assert success is False
    assert mutation.replace("_", " ").split()[0] in message.lower()


def test_cclib_compatibility_probe_is_cached_once_per_process(monkeypatch, fake_cclib_api):
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: fake_cclib_api)

    first = cclib_harness._check_cclib_compatibility()
    second = cclib_harness._check_cclib_compatibility()
    assert first == second
    assert first[0] is True
    assert fake_cclib_api.QCSchemaWriter.calls == 1


@pytest.mark.addon
def test_real_cclib_compatibility_when_installed():
    pytest.importorskip("cclib")

    success, version_or_message = cclib_harness._check_cclib_compatibility()
    assert success is True, version_or_message


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
    monkeypatch.setattr(cclib_harness, "_check_cclib_compatibility", lambda: (True, "1.9"))
    monkeypatch.setattr(cclib_harness, "which", lambda command: None, raising=False)

    assert harness.found() is False
    with pytest.raises(ResourceError, match="orca.*PATH"):
        harness.found(raise_error=True)


@pytest.mark.parametrize(
    "program,path,output,expected",
    [
        ("qchem", "/opt/qchem-5.1/bin/qchem", _qchem_probe_output("5.1"), "5.1"),
        ("qchem", "/opt/qchem-6.2/bin/qchem", _qchem_probe_output("6.2.2"), "6.2.2"),
        ("orca", "/opt/orca-6/bin/orca", _orca_probe_output("6.0"), "6.0.0"),
        ("orca", "/opt/orca-6.1/bin/orca", _orca_probe_output("6.1"), "6.1.0"),
    ],
)
def test_executable_identity_and_supported_version(monkeypatch, program, path, output, expected):
    calls = []

    def fake_execute(command, infiles=None, **kwargs):
        calls.append((command, infiles, kwargs))
        return True, {"stdout": output, "stderr": ""}

    monkeypatch.setattr(cclib_harness, "execute", fake_execute, raising=False)
    harness = CCLibHarness(name=f"cclib-{program}", program=program)

    assert cclib_harness._probe_executable(harness, path) == expected
    assert harness.version_cache[path] == expected
    assert len(calls) == 1
    if program == "qchem":
        assert "$rem" in next(iter(calls[0][1].values()))
    else:
        assert "HF STO-3G" in next(iter(calls[0][1].values()))


@pytest.mark.parametrize(
    "program,output,match",
    [
        ("qchem", "Q-Chem 6.2 for Linux\n", "identity"),
        ("qchem", _qchem_probe_output("5.0"), "requires.*5.1"),
        ("qchem", _qchem_probe_output("not-a-version"), "version"),
        ("orca", "Program Version 6.1.0\nORCA TERMINATED NORMALLY\n", "identity"),
        ("orca", _orca_probe_output("5.0"), "requires.*6.0"),
        ("orca", "O   R   C   A\nProgram Version unknown\nORCA TERMINATED NORMALLY\n", "version"),
        ("orca", "O   R   C   A\nProgram Version 6.1.0\n", "normal termination"),
    ],
)
def test_executable_identity_or_version_rejection(monkeypatch, program, output, match):
    monkeypatch.setattr(
        cclib_harness,
        "execute",
        lambda *args, **kwargs: (True, {"stdout": output, "stderr": ""}),
        raising=False,
    )
    harness = CCLibHarness(name=f"cclib-{program}", program=program)
    path = f"/opt/{program}"

    with pytest.raises(ResourceError, match=match):
        cclib_harness._probe_executable(harness, path)
    assert path not in harness.version_cache


def test_unrelated_orca_executable_is_rejected(monkeypatch):
    output = "Orca is a screen reader and magnifier for the GNOME desktop.\n"
    monkeypatch.setattr(
        cclib_harness,
        "execute",
        lambda *args, **kwargs: (True, {"stdout": output, "stderr": ""}),
        raising=False,
    )
    harness = CCLibHarness(name="cclib-orca", program="orca")

    with pytest.raises(ResourceError, match="identity"):
        cclib_harness._probe_executable(harness, "/usr/bin/orca")


def test_executable_probe_cache_is_keyed_by_resolved_path(monkeypatch):
    calls = []

    def fake_execute(command, infiles=None, **kwargs):
        calls.append(command[0])
        version = "5.1" if "first" in command[0] else "6.2"
        return True, {"stdout": _qchem_probe_output(version), "stderr": ""}

    monkeypatch.setattr(cclib_harness, "execute", fake_execute, raising=False)
    harness = CCLibHarness(name="cclib-qchem", program="qchem")

    assert cclib_harness._probe_executable(harness, "/opt/first/qchem") == "5.1"
    assert cclib_harness._probe_executable(harness, "/opt/first/qchem") == "5.1"
    assert cclib_harness._probe_executable(harness, "/opt/second/qchem") == "6.2"
    assert calls == ["/opt/first/qchem", "/opt/second/qchem"]


def test_get_version_returns_external_program_version(monkeypatch):
    harness = CCLibHarness(name="cclib-orca", program="orca")
    monkeypatch.setattr(cclib_harness, "which", lambda command: "/opt/orca", raising=False)
    monkeypatch.setattr(
        cclib_harness,
        "execute",
        lambda *args, **kwargs: (True, {"stdout": _orca_probe_output("6.1"), "stderr": ""}),
        raising=False,
    )

    assert harness.get_version() == "6.1.0"


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


def test_found_checks_resources_in_required_order(monkeypatch):
    events = []
    harness = CCLibHarness(name="cclib-qchem", program="qchem")

    monkeypatch.setattr(
        cclib_harness,
        "_check_cclib_compatibility",
        lambda: events.append("compatibility") or (True, "1.9"),
    )
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
    assert events == ["compatibility", "path", "preflight", "identity/version"]


@pytest.mark.parametrize(
    "stage,expected",
    [
        ("compatibility", "writer feature unavailable"),
        ("path", "PATH"),
        ("preflight", "QCAUX"),
        ("probe", "identity"),
    ],
)
def test_found_false_suppresses_resource_failures_and_found_true_preserves_detail(monkeypatch, stage, expected):
    harness = CCLibHarness(name="cclib-qchem", program="qchem")
    monkeypatch.setattr(cclib_harness, "_check_cclib_compatibility", lambda: (True, "1.9"))
    monkeypatch.setattr(cclib_harness, "which", lambda command: "/opt/qchem")
    definition = replace(
        cclib_harness._PROGRAM_DEFINITIONS["qchem"],
        preflight=lambda executable, environment: environment,
    )
    monkeypatch.setitem(cclib_harness._PROGRAM_DEFINITIONS, "qchem", definition)
    monkeypatch.setattr(cclib_harness, "_probe_executable", lambda *args: "5.1")

    if stage == "compatibility":
        monkeypatch.setattr(
            cclib_harness,
            "_check_cclib_compatibility",
            lambda: (False, "writer feature unavailable"),
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
        "_check_cclib_compatibility",
        lambda: (False, "cclib is not importable"),
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
        ccData=object,
        QCSchemaWriter=FakeWriter,
        ccread=ccopen,
        QChem=FakeQChem,
        ORCA=FakeORCA,
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
    return api, cclib_harness._PROGRAM_DEFINITIONS[program], execution, opened


def _assert_conversion_failure(monkeypatch, match, **case):
    api, definition, execution, _ = _fake_conversion_case(**case)
    monkeypatch.setattr(cclib_harness, "_load_cclib_api", lambda: api)
    with pytest.raises(UnknownError, match=match) as exc_info:
        cclib_harness._parse_and_convert(definition, execution, _atomic_input())
    assert definition.selector in str(exc_info.value)
    assert exc_info.value.__cause__ is not None
    return exc_info.value


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
    monkeypatch.setattr(cclib_harness, "_check_cclib_compatibility", lambda: (True, api.version))
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


def _qchem_demonstration_result():
    geometry = [
        -0.0,
        0.0,
        0.22517858316070177,
        -1.4941103633283772,
        -0.0,
        -0.9007143324538345,
        1.4941103633283772,
        -0.0,
        -0.9007143324538345,
    ]
    return {
        "schema_name": "qcschema_output",
        "schema_version": 1,
        "success": True,
        "driver": "energy",
        "model": {"method": "mp2", "basis": "sto-3g"},
        "molecule": {
            "symbols": ["O", "H", "H"],
            "geometry": geometry,
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        },
        "provenance": {
            "creator": "QChem",
            "version": "5.1.2",
            "routine": "cclib.io.qcschemawriter.QCSchemaWriter",
        },
        "properties": {
            "calcinfo_nbasis": 7,
            "calcinfo_nmo": 7,
            "calcinfo_nalpha": 5,
            "calcinfo_nbeta": 5,
            "calcinfo_natom": 3,
            "return_energy": -75.00228214,
            "scf_dipole_moment": [0.0, 0.0, -0.6584056190],
            "scf_total_energy": -74.9643287618,
            "scf_iterations": 6,
            "mp2_correlation_energy": -0.0379533782,
            "mp2_total_energy": -75.00228214,
        },
        "return_result": -75.00228214,
        "extras": {
            "atomcharges": {"mulliken": [-0.339215, 0.169607, 0.169607]},
            "atomcoords": [[geometry[index : index + 3] for index in range(0, 9, 3)]],
            "atomnos": [8, 1, 1],
            "homos": [4],
            "moenergies": [[-20.244, -1.251, -0.603, -0.445, -0.388, 0.571, 0.709]],
            "mosyms": [["A1", "A1", "B1", "A1", "B2", "A1", "B1"]],
            "mpenergies": [[-75.00228214]],
            "scfenergies": [-74.9643287618],
            "scftargets": [[1.0e-5]],
            "scfvalues": [[[0.398], [0.0668], [0.00822], [0.0016], [2.83e-5], [8.23e-6]]],
            "cclib_harness": {
                "selector": "cclib-qchem",
                "cclib_version": "1.9.dev",
                "parser": "QChem",
                "executable": "/resolved/qchem",
            },
        },
    }


def _orca_demonstration_result(version="6.0.1"):
    geometry = [
        3.372998617495434,
        2.385631834752722,
        0.9675114303425262,
        5.004442645304063,
        2.027541962061342,
        0.24874653962013937,
        2.2358634804056874,
        2.3750380300934055,
        -0.45133273917372035,
    ]
    result = {
        "schema_name": "qcschema_output",
        "schema_version": 1,
        "success": True,
        "driver": "energy",
        "model": {"method": "ccsd", "basis": "sto-3g"},
        "molecule": {
            "symbols": ["O", "H", "H"],
            "geometry": geometry,
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
        },
        "provenance": {
            "creator": "ORCA",
            "version": version,
            "routine": "cclib.io.qcschemawriter.QCSchemaWriter",
        },
        "properties": {
            "calcinfo_nbasis": 7,
            "calcinfo_nmo": 7,
            "calcinfo_nalpha": 5,
            "calcinfo_nbeta": 5,
            "calcinfo_natom": 3,
            "return_energy": -75.013487814,
            "scf_total_energy": -74.96357424008319,
            "ccsd_correlation_energy": -0.04991357391681104,
            "ccsd_total_energy": -75.013487814,
            "mp2_correlation_energy": 74.9277738680832,
            "mp2_total_energy": -0.035800372,
        },
        "return_result": -75.013487814,
        "extras": {
            "atomcharges": {
                "mulliken": [-0.329397, 0.164693, 0.164703],
                "lowdin": [-0.222995, 0.111495, 0.1115],
            },
            "atomcoords": [
                [
                    [3.372998615785793, 2.385631833543539, 0.9675114298521326],
                    [5.004442642767507, 2.0275419610336605, 0.2487465394940595],
                    [2.235863479272416, 2.375038028889592, -0.4513327389449575],
                ]
            ],
            "atomnos": [8, 1, 1],
            "ccenergies": [-75.013487814],
            "homos": [4],
            "moenergies": [
                [
                    -20.242268999999997,
                    -1.2657839999999998,
                    -0.615347,
                    -0.452279,
                    -0.39087700000000003,
                    0.60058,
                    0.736585,
                ]
            ],
            "mosyms": [["A", "A", "A", "A", "A", "A", "A"]],
            "mpenergies": [[-0.035800372]],
            "scfenergies": [-74.96357424008319],
            "scftargets": [[1.0e-6, 1.0e-5, 1.0e-6]],
            "scfvalues": [
                [
                    [0.0, 0.0263, 0.0744],
                    [-0.0179, 0.0225, 0.0624],
                    [-0.0127, 0.0157, 0.0433],
                    [-0.0087, 0.0367, 0.101],
                    [-0.0199, 0.00115, 0.00458],
                    [-8.88e-6, 0.0006, 0.00221],
                    [-1.69e-6, 0.000344, 0.00117],
                    [-2.64e-7, 2.42e-5, 6.93e-5],
                    [2.6405e-7, 6.9252e-5, 2.4156e-5],
                ]
            ],
            "cclib_harness": {
                "selector": "cclib-orca",
                "cclib_version": "1.9.dev",
                "parser": "ORCA",
                "executable": "/resolved/orca",
            },
        },
    }
    if str(version).startswith("6.1"):
        result["extras"].update(
            {
                "atomcharges": {
                    "mulliken": [-0.329397, 0.164693, 0.164703],
                    "lowdin": [-0.222995, 0.111495, 0.1115],
                    "hirshfeld": [-0.288291, 0.144144, 0.144146],
                },
                "moenergies": [[-20.242272, -1.265785, -0.615354, -0.452275, -0.390879, 0.600583, 0.736578]],
                "scfenergies": [-74.96357424464694],
                "scfvalues": [
                    [
                        [0.0, 0.0568, 0.0744],
                        [-0.0179, 0.0486, 0.0624],
                        [-0.0127, 0.0339, 0.0433],
                        [-0.0087, 0.0793, 0.101],
                        [-0.0199, 0.00247, 0.00458],
                        [-8.88e-6, 0.0013, 0.00221],
                        [-1.69e-6, 0.000744, 0.00117],
                        [-2.64e-7, 5.22e-5, 6.93e-5],
                        [2.6405e-7, 6.9252e-5, 5.2183e-5],
                    ]
                ],
            }
        )
    return result


def test_qchem_demonstration_comparison_covers_all_acceptance_criteria():
    demonstration = importlib.import_module("qchem_water_mp2")
    comparisons = demonstration.compare_result(_qchem_demonstration_result())

    assert comparisons
    assert all(line.startswith("PASS ") for line in comparisons), comparisons
    labels = "\n".join(comparisons)
    for required in (
        "return_result",
        "return_energy",
        "scf_total_energy",
        "mp2_total_energy",
        "mp2_correlation_energy",
        "scf_dipole_moment",
        "calcinfo_nbasis",
        "calcinfo_nmo",
        "calcinfo_nalpha",
        "calcinfo_nbeta",
        "calcinfo_natom",
        "scf_iterations",
        "schema identity",
        "model",
        "molecule",
        "provenance",
        "cclib-qchem metadata",
        "flat extras",
        "Mulliken charges",
        "seven MO energies",
        "six SCF rows",
    ):
        assert required in labels


@pytest.mark.parametrize(
    "path,bad_value,failed_label",
    [
        (("return_result",), -1.0, "return_result"),
        (("properties", "scf_dipole_moment"), [0.0, 0.0, 0.0], "scf_dipole_moment"),
        (("properties", "calcinfo_nbasis"), 8, "calcinfo_nbasis"),
        (("molecule", "symbols"), ["H", "O", "H"], "molecule"),
        (("provenance", "routine"), "other.writer", "provenance"),
        (("extras", "atomcharges", "mulliken"), [0.0, 0.0, 0.0], "Mulliken charges"),
        (("extras", "moenergies"), [[0.0] * 7], "seven MO energies"),
        (
            ("extras", "scfvalues"),
            [[[0.398], [0.0668], [0.00822], [0.0016], [3.03e-5], [8.23e-6]]],
            "six SCF rows",
        ),
    ],
)
def test_qchem_demonstration_comparison_reports_failures(path, bad_value, failed_label):
    demonstration = importlib.import_module("qchem_water_mp2")
    result = _qchem_demonstration_result()
    target = result
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = bad_value

    failures = [line for line in demonstration.compare_result(result) if line.startswith("FAIL ")]
    assert any(failed_label in line for line in failures)


def test_orca_demonstration_comparison_checks_ccsd_scf_provenance_extras_and_anomalous_mp2():
    demonstration = importlib.import_module("orca_water_ccsd")
    comparisons = demonstration.compare_result(_orca_demonstration_result())

    assert comparisons
    assert all(line.startswith("PASS ") for line in comparisons), comparisons
    labels = "\n".join(comparisons)
    for required in (
        "return_result",
        "ccsd_total_energy",
        "scf_total_energy",
        "ccsd_correlation_energy",
        "calcinfo_nbasis",
        "calcinfo_nmo",
        "calcinfo_nalpha",
        "calcinfo_nbeta",
        "calcinfo_natom",
        "provenance",
        "cclib-orca metadata",
        "flat CCSD extras",
        "orbital extras",
        "SCF extras",
        "anomalous MP2 writer fields present but not numerically endorsed",
    ):
        assert required in labels


@pytest.mark.parametrize("version", ["6.0.1", "6.1.1"])
@pytest.mark.parametrize(
    "path,bad_value,failed_label",
    [
        (("extras", "moenergies", 0, 3), -0.452277, "orbital extras"),
        (("extras", "mosyms", 0, 3), "B", "orbital extras"),
        (("extras", "scftargets", 0, 1), 1.2e-5, "SCF extras"),
        (("extras", "scfvalues", 0, 4, 1), 0.001152, "SCF extras"),
        (("extras", "atomcharges", "lowdin", 1), 0.111497, "flat CCSD extras"),
        (("extras", "atomcoords", 0, 1, 2), 0.2487485394940595, "flat CCSD extras"),
        (("extras", "atomnos", 2), 2, "flat CCSD extras"),
        (("extras", "ccenergies", 0), -75.013485814, "flat CCSD extras"),
    ],
)
def test_orca_demonstration_comparison_reports_exact_extra_failures(version, path, bad_value, failed_label):
    demonstration = importlib.import_module("orca_water_ccsd")
    result = _orca_demonstration_result(version)
    target = result
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = bad_value

    failures = [line for line in demonstration.compare_result(result) if line.startswith("FAIL ")]
    assert any(failed_label in line for line in failures)


@pytest.mark.parametrize("version", ["6.0.1", "6.0.9", "6.1.1", "6.1.99"])
def test_orca_demonstration_selects_strict_versioned_reference(version):
    demonstration = importlib.import_module("orca_water_ccsd")
    comparisons = demonstration.compare_result(_orca_demonstration_result(version))

    assert all(line.startswith("PASS ") for line in comparisons), comparisons
    assert any("ORCA extras reference version" in line for line in comparisons)


@pytest.mark.parametrize(
    "version",
    [
        None,
        "",
        "6",
        "6.1",
        "6.1.not-a-version",
        "6.1.1.2",
        " 6.1.1",
        "6.1.1 ",
        "v6.1.1",
        "6.1.1dev",
        "6.1.+1",
        "6.2.0",
        "not-a-version",
    ],
)
def test_orca_demonstration_rejects_missing_or_unsupported_reference_version(version):
    demonstration = importlib.import_module("orca_water_ccsd")
    result = _orca_demonstration_result()
    result["provenance"]["version"] = version

    failures = [line for line in demonstration.compare_result(result) if line.startswith("FAIL ")]
    assert any("ORCA extras reference version" in line for line in failures)


def test_demonstration_main_calls_compute_writes_complete_json_and_reports(monkeypatch, tmp_path, capsys):
    cases = [
        ("qchem_water_mp2", _qchem_demonstration_result(), "cclib-qchem", None),
        (
            "orca_water_ccsd",
            _orca_demonstration_result(),
            "cclib-orca",
            {"ncores": 4, "memory": 2.734375},
        ),
    ]
    for module_name, result, selector, task_config in cases:
        demonstration = importlib.import_module(module_name)
        calls = []

        def fake_compute(atomic_input, program, **kwargs):
            calls.append((atomic_input, program, kwargs))
            return copy.deepcopy(result)

        monkeypatch.setattr(demonstration.qcengine, "compute", fake_compute)
        output_path = tmp_path / f"{module_name}.json"
        assert demonstration.main(output_path=output_path) == 0
        assert json.loads(output_path.read_text()) == result
        assert len(calls) == 1
        atomic_input, actual_selector, kwargs = calls[0]
        assert actual_selector == selector
        assert kwargs["raise_error"] is True
        assert kwargs["return_version"] == 1
        if task_config is None:
            assert "task_config" not in kwargs
        else:
            assert kwargs["task_config"] == task_config
        assert atomic_input.specification.driver.value == "energy"
        assert atomic_input.specification.model.basis == "sto-3g"
        assert "FAIL " not in capsys.readouterr().out


def test_demonstration_serializes_qcschema_v1_models_without_pydantic_v2_mode():
    demonstration = importlib.import_module("qchem_water_mp2")
    result = _qchem_demonstration_result()

    class V1Result:
        def model_dump(self, **kwargs):
            if "mode" in kwargs:
                raise TypeError("dict() got an unexpected keyword argument 'mode'")
            return result

        def json(self):
            return json.dumps(result)

    assert demonstration._jsonable(V1Result()) == result


def test_demonstration_main_returns_nonzero_when_a_comparison_fails(monkeypatch, tmp_path, capsys):
    demonstration = importlib.import_module("qchem_water_mp2")
    result = _qchem_demonstration_result()
    result["return_result"] = 0.0
    monkeypatch.setattr(demonstration.qcengine, "compute", lambda *args, **kwargs: result)

    assert demonstration.main(output_path=tmp_path / "failed.json") != 0
    assert "FAIL return_result" in capsys.readouterr().out


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


@uusing("cclib-orca")
def test_live_cclib_orca_water_mp2_energy():
    ccsd_input = importlib.import_module("orca_water_ccsd").build_atomic_input()
    input_model = AtomicInput(
        molecule=ccsd_input.molecule,
        specification={
            "driver": "energy",
            "model": {"method": "mp2", "basis": "sto-3g"},
            "keywords": {},
        },
    )
    result = qcng.compute(
        input_model,
        "cclib-orca",
        raise_error=True,
        task_config={"ncores": 4, "memory": 2.734375},
        return_version=1,
    )

    assert result.success is True
    assert float(result.return_result) == pytest.approx(-74.999371925, abs=5.0e-6)


@uusing("cclib-orca")
def test_live_cclib_orca_water_ccsd_demonstration():
    demonstration = importlib.import_module("orca_water_ccsd")
    result = qcng.compute(
        demonstration.build_atomic_input(),
        "cclib-orca",
        raise_error=True,
        task_config={"ncores": 4, "memory": 2.734375},
        return_version=1,
        return_dict=True,
    )

    assert not [line for line in demonstration.compare_result(result) if line.startswith("FAIL ")]


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
    assert type(initial_parser) is getattr(api, cclib_harness._PROGRAM_DEFINITIONS[program].parser_name)
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
