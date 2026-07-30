"""Run and validate the live cclib-backed Q-Chem water MP2 demonstration."""

import json
import math
import sys
from pathlib import Path
from typing import Any, List

import qcengine
from qcelemental.models.v2 import AtomicInput


_GEOMETRY = [
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
_MO_ENERGIES = [-20.244, -1.251, -0.603, -0.445, -0.388, 0.571, 0.709]


def build_atomic_input() -> AtomicInput:
    """Construct the demonstration input directly from QCSchema data."""
    return AtomicInput(
        molecule={
            "symbols": ["O", "H", "H"],
            "geometry": _GEOMETRY,
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
            "fix_com": True,
            "fix_orientation": True,
        },
        specification={
            "driver": "energy",
            "model": {"method": "mp2", "basis": "sto-3g"},
            "keywords": {},
        },
    )


def _get(result: Any, *path: str) -> Any:
    value = result
    try:
        for key in path:
            value = value[key] if isinstance(value, dict) else getattr(value, key)
    except (AttributeError, KeyError, TypeError):
        return None
    return value


def _close(actual: Any, expected: float, tolerance: float) -> bool:
    try:
        return math.isclose(float(actual), expected, rel_tol=0.0, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def _close_sequence(actual: Any, expected: List[float], tolerance: float) -> bool:
    try:
        return len(actual) == len(expected) and all(
            _close(value, reference, tolerance) for value, reference in zip(actual, expected)
        )
    except TypeError:
        return False


def _line(label: str, passed: bool) -> str:
    return f"{'PASS' if passed else 'FAIL'} {label}"


def compare_result(result: Any) -> List[str]:
    """Compare a QCSchema result with every approved demonstration criterion."""
    properties = _get(result, "properties") or {}
    extras = _get(result, "extras") or {}
    molecule = _get(result, "molecule") or {}
    model = _get(result, "model") or {}
    provenance = _get(result, "provenance") or {}

    checks = [
        ("return_result", _close(_get(result, "return_result"), -75.00228214, 1.0e-6)),
        ("properties.return_energy", _close(_get(properties, "return_energy"), -75.00228214, 1.0e-6)),
        (
            "properties.scf_total_energy",
            _close(_get(properties, "scf_total_energy"), -74.9643287618, 1.0e-6),
        ),
        ("properties.mp2_total_energy", _close(_get(properties, "mp2_total_energy"), -75.00228214, 1.0e-6)),
        (
            "properties.mp2_correlation_energy",
            _close(_get(properties, "mp2_correlation_energy"), -0.0379533782, 1.0e-6),
        ),
        (
            "properties.scf_dipole_moment",
            _close_sequence(_get(properties, "scf_dipole_moment"), [0.0, 0.0, -0.6584056190], 1.0e-5),
        ),
    ]
    for field, expected in (
        ("calcinfo_nbasis", 7),
        ("calcinfo_nmo", 7),
        ("calcinfo_nalpha", 5),
        ("calcinfo_nbeta", 5),
        ("calcinfo_natom", 3),
        ("scf_iterations", 6),
    ):
        checks.append((f"properties.{field}", _get(properties, field) == expected))

    schema_ok = (
        _get(result, "success") is True
        and _get(result, "schema_name") == "qcschema_output"
        and _get(result, "schema_version") == 1
        and str(_get(result, "driver")).lower() == "energy"
    )
    checks.append(("schema identity", schema_ok))
    checks.append(
        (
            "model mp2/sto-3g",
            str(_get(model, "method")).lower() == "mp2" and str(_get(model, "basis")).lower() == "sto-3g",
        )
    )
    molecule_ok = (
        list(_get(molecule, "symbols") or []) == ["O", "H", "H"]
        and _close(_get(molecule, "molecular_charge"), 0.0, 0.0)
        and _get(molecule, "molecular_multiplicity") == 1
        and _close_sequence(_get(molecule, "geometry"), _GEOMETRY, 1.0e-8)
    )
    checks.append(("molecule neutral singlet O/H/H geometry", molecule_ok))
    checks.append(
        (
            "provenance QChem cclib QCSchemaWriter",
            _get(provenance, "creator") == "QChem"
            and "qcschemawriter" in str(_get(provenance, "routine")).lower(),
        )
    )
    checks.append(
        (
            "cclib-qchem metadata",
            _get(extras, "cclib_harness", "selector") == "cclib-qchem",
        )
    )

    required_extras = {
        "atomcharges",
        "atomcoords",
        "atomnos",
        "homos",
        "moenergies",
        "mosyms",
        "mpenergies",
        "scfenergies",
        "scftargets",
        "scfvalues",
    }
    checks.append(("representative flat extras", required_extras <= set(extras)))
    checks.append(
        (
            "Mulliken charges",
            _close_sequence(_get(extras, "atomcharges", "mulliken"), [-0.339215, 0.169607, 0.169607], 1.0e-5),
        )
    )
    moenergies = _get(extras, "moenergies")
    checks.append(
        (
            "seven MO energies",
            isinstance(moenergies, (list, tuple))
            and len(moenergies) == 1
            and _close_sequence(moenergies[0], _MO_ENERGIES, 1.0e-3),
        )
    )
    scfvalues = _get(extras, "scfvalues")
    checks.append(
        (
            "six SCF rows",
            isinstance(scfvalues, (list, tuple))
            and len(scfvalues) == 1
            and isinstance(scfvalues[0], (list, tuple))
            and len(scfvalues[0]) == 6,
        )
    )
    return [_line(label, passed) for label, passed in checks]


def _jsonable(result: Any) -> Any:
    if isinstance(result, dict):
        return result
    if hasattr(result, "model_dump"):
        try:
            return result.model_dump(mode="json")
        except TypeError:
            pass
    return json.loads(result.json())


def main(output_path: Any = "qchem_water_mp2.result.json") -> int:
    result = qcengine.compute(
        build_atomic_input(),
        "cclib-qchem",
        raise_error=True,
        return_version=1,
    )
    data = _jsonable(result)
    Path(output_path).write_text(json.dumps(data, indent=2) + "\n")
    comparisons = compare_result(data)
    for comparison in comparisons:
        print(comparison)
    return 1 if any(line.startswith("FAIL ") for line in comparisons) else 0


if __name__ == "__main__":
    sys.exit(main())
