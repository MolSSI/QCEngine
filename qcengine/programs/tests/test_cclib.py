import json
import os
import subprocess
import sys
from dataclasses import replace

import pytest

import qcengine as qcng
import qcengine.programs.cclib as cclib_harness
from qcelemental.models.v2 import AtomicInput, BasisSet

from qcengine.config import TaskConfig
from qcengine.exceptions import InputError, ResourceError
from qcengine.programs.cclib import CCLibHarness


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


def _atomic_input(driver="energy", method="hf", basis="sto-3g", keywords=None, molecule=None):
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
