import json
import os
import subprocess
import sys
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

import qcengine as qcng
import qcengine.programs.cclib as cclib_harness
from qcengine.exceptions import ResourceError
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


def _synthetic_cclib_data(ccdata):
    return ccdata(
        {
            "atomcoords": np.array([[[1.0, 0.0, 0.0]]]),
            "atomnos": np.array([2]),
            "atomcharges": {"mulliken": np.array([0.0])},
            "charge": 0,
            "homos": np.array([0]),
            "metadata": {
                "package": "Synthetic",
                "package_version": "1.0",
                "methods": ["HF"],
                "basis_set": "sto-3g",
                "success": True,
            },
            "mult": 1,
            "natom": 1,
            "nbasis": 1,
            "nmo": 1,
            "scfenergies": np.array([-24.6]),
            "scftargets": np.array([[[1.0e-6]]]),
            "scfvalues": np.array([[[1.0e-4], [1.0e-7]]]),
        }
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


def test_cclib_compatibility_validates_v1_geometry_and_flat_extras():
    success, version = cclib_harness._check_cclib_compatibility()
    assert success is True
    assert version

    api = cclib_harness._load_cclib_api()
    output = api.QCSchemaWriter(_synthetic_cclib_data(api.ccData)).as_dict(validate=False)
    from qcelemental.models.v1 import AtomicResult

    result = AtomicResult(**output)
    assert result.schema_version == 1
    assert result.molecule.geometry[0][0] == pytest.approx(1.8897261255, abs=2.0e-9)
    assert result.extras["atomcoords"][0][0][0] == pytest.approx(1.8897261255, abs=2.0e-9)
    assert result.extras["atomcharges"] == {"mulliken": [0.0]}


@pytest.mark.parametrize("mutation", ["wrong_units", "missing_atomcoords", "missing_atomcharges"])
def test_cclib_compatibility_rejects_incomplete_writer(monkeypatch, mutation):
    real_api = cclib_harness._load_cclib_api()

    class BrokenWriter:
        calls = 0

        def __init__(self, data):
            self.data = data

        def as_dict(self, validate=False):
            BrokenWriter.calls += 1
            output = real_api.QCSchemaWriter(self.data).as_dict(validate=validate)
            if mutation == "wrong_units":
                output["molecule"]["geometry"][0] = 1.0
                output["extras"]["atomcoords"][0][0][0] = 1.0
            elif mutation == "missing_atomcoords":
                output["extras"].pop("atomcoords")
            else:
                output["extras"].pop("atomcharges")
            return output

    monkeypatch.setattr(
        cclib_harness,
        "_load_cclib_api",
        lambda: SimpleNamespace(
            version=real_api.version,
            ccData=real_api.ccData,
            QCSchemaWriter=BrokenWriter,
            ccread=real_api.ccread,
            QChem=real_api.QChem,
            ORCA=real_api.ORCA,
        ),
    )
    success, message = cclib_harness._check_cclib_compatibility()

    assert success is False
    assert mutation.replace("_", " ").split()[0] in message.lower()


def test_cclib_compatibility_probe_is_cached_once_per_process(monkeypatch):
    real_api = cclib_harness._load_cclib_api()

    class CountingWriter(real_api.QCSchemaWriter):
        calls = 0

        def as_dict(self, validate=False):
            CountingWriter.calls += 1
            return super().as_dict(validate=validate)

    monkeypatch.setattr(
        cclib_harness,
        "_load_cclib_api",
        lambda: SimpleNamespace(
            version=real_api.version,
            ccData=real_api.ccData,
            QCSchemaWriter=CountingWriter,
            ccread=real_api.ccread,
            QChem=real_api.QChem,
            ORCA=real_api.ORCA,
        ),
    )

    first = cclib_harness._check_cclib_compatibility()
    second = cclib_harness._check_cclib_compatibility()
    assert first == second
    assert first[0] is True
    assert CountingWriter.calls == 1


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
