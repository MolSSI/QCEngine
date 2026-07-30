"""Run and validate the live cclib-backed ORCA water CCSD demonstration."""

import json
import math
import re
import sys
from pathlib import Path
from typing import Any, List

import qcengine
from qcelemental.models.v2 import AtomicInput


# The water_ccsd.inp Cartesian coordinates converted from Angstrom to bohr.
_GEOMETRY_BOHR = [
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
_ATOMCHARGES = {
    "mulliken": [-0.329397, 0.164693, 0.164703],
    "lowdin": [-0.222995, 0.111495, 0.1115],
}
_ATOMCOORDS = [
    [
        [3.372998615785793, 2.385631833543539, 0.9675114298521326],
        [5.004442642767507, 2.0275419610336605, 0.2487465394940595],
        [2.235863479272416, 2.375038028889592, -0.4513327389449575],
    ]
]
_ATOMNOS = [8, 1, 1]
_CCENERGIES = [-75.013487814]
_HOMOS = [4]
_ORBITAL_ENERGIES = [
    [
        -20.242268999999997,
        -1.2657839999999998,
        -0.615347,
        -0.452279,
        -0.39087700000000003,
        0.60058,
        0.736585,
    ]
]
_ORBITAL_SYMMETRIES = [["A", "A", "A", "A", "A", "A", "A"]]
_SCFENERGIES = [-74.96357424008319]
_SCFTARGETS = [[1.0e-6, 1.0e-5, 1.0e-6]]
_SCFVALUES = [
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
]
_ATOMCHARGES_61 = {
    "mulliken": [-0.329397, 0.164693, 0.164703],
    "lowdin": [-0.222995, 0.111495, 0.1115],
    "hirshfeld": [-0.288291, 0.144144, 0.144146],
}
_ORBITAL_ENERGIES_61 = [[-20.242272, -1.265785, -0.615354, -0.452275, -0.390879, 0.600583, 0.736578]]
_SCFENERGIES_61 = [-74.96357424464694]
_SCFVALUES_61 = [
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
]
_ORCA_EXTRA_REFERENCES = {
    (6, 0): {
        "flat": {
            "atomcharges": _ATOMCHARGES,
            "atomcoords": _ATOMCOORDS,
            "atomnos": _ATOMNOS,
            "ccenergies": _CCENERGIES,
        },
        "orbital": {
            "homos": _HOMOS,
            "moenergies": _ORBITAL_ENERGIES,
            "mosyms": _ORBITAL_SYMMETRIES,
        },
        "scf": {
            "scfenergies": _SCFENERGIES,
            "scftargets": _SCFTARGETS,
            "scfvalues": _SCFVALUES,
        },
    },
    (6, 1): {
        "flat": {
            "atomcharges": _ATOMCHARGES_61,
            "atomcoords": _ATOMCOORDS,
            "atomnos": _ATOMNOS,
            "ccenergies": _CCENERGIES,
        },
        "orbital": {
            "homos": _HOMOS,
            "moenergies": _ORBITAL_ENERGIES_61,
            "mosyms": _ORBITAL_SYMMETRIES,
        },
        "scf": {
            "scfenergies": _SCFENERGIES_61,
            "scftargets": _SCFTARGETS,
            "scfvalues": _SCFVALUES_61,
        },
    },
}
_ORCA_VERSION = re.compile(r"^([0-9]+)\.([0-9]+)\.([0-9]+)$")
_TASK_CONFIG = {"ncores": 4, "memory": 2.734375}


def build_atomic_input() -> AtomicInput:
    """Construct the historical ORCA CCSD calculation directly as QCSchema."""
    return AtomicInput(
        molecule={
            "symbols": ["O", "H", "H"],
            "geometry": _GEOMETRY_BOHR,
            "molecular_charge": 0,
            "molecular_multiplicity": 1,
            "fix_com": True,
            "fix_orientation": True,
        },
        specification={
            "driver": "energy",
            "model": {"method": "ccsd", "basis": "sto-3g"},
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


def _matches(actual: Any, expected: Any, tolerance: float) -> bool:
    if isinstance(expected, dict):
        return (
            isinstance(actual, dict)
            and set(actual) == set(expected)
            and all(_matches(actual[key], value, tolerance) for key, value in expected.items())
        )
    if isinstance(expected, (list, tuple)):
        return (
            isinstance(actual, (list, tuple))
            and len(actual) == len(expected)
            and all(_matches(value, reference, tolerance) for value, reference in zip(actual, expected))
        )
    if isinstance(expected, float):
        return _close(actual, expected, tolerance)
    return actual == expected


def _orca_reference(version: Any) -> Any:
    """Validate major.minor.patch provenance, then select the strict major.minor reference."""
    if not isinstance(version, str):
        return None
    match = _ORCA_VERSION.fullmatch(version)
    if match is None:
        return None
    major, minor, _patch = (int(component) for component in match.groups())
    return _ORCA_EXTRA_REFERENCES.get((major, minor))


def _line(label: str, passed: bool) -> str:
    return f"{'PASS' if passed else 'FAIL'} {label}"


def compare_result(result: Any) -> List[str]:
    """Compare a QCSchema result with the approved ORCA CCSD criteria."""
    properties = _get(result, "properties") or {}
    extras = _get(result, "extras") or {}
    molecule = _get(result, "molecule") or {}
    model = _get(result, "model") or {}
    provenance = _get(result, "provenance") or {}

    checks = [
        ("return_result", _close(_get(result, "return_result"), -75.013487814, 1.0e-6)),
        (
            "properties.ccsd_total_energy",
            _close(_get(properties, "ccsd_total_energy"), -75.013487814, 1.0e-6),
        ),
        (
            "properties.scf_total_energy",
            _close(_get(properties, "scf_total_energy"), -74.96357424008319, 1.0e-6),
        ),
        (
            "properties.ccsd_correlation_energy",
            _close(_get(properties, "ccsd_correlation_energy"), -0.04991357391681104, 1.0e-6),
        ),
    ]
    for field, expected in (
        ("calcinfo_nbasis", 7),
        ("calcinfo_nmo", 7),
        ("calcinfo_nalpha", 5),
        ("calcinfo_nbeta", 5),
        ("calcinfo_natom", 3),
    ):
        checks.append((f"properties.{field}", _get(properties, field) == expected))

    checks.append(
        (
            "schema/model/molecule",
            _get(result, "success") is True
            and _get(result, "schema_name") == "qcschema_output"
            and _get(result, "schema_version") == 1
            and str(_get(result, "driver")).lower() == "energy"
            and str(_get(model, "method")).lower() == "ccsd"
            and str(_get(model, "basis")).lower() == "sto-3g"
            and list(_get(molecule, "symbols") or []) == ["O", "H", "H"]
            and _get(molecule, "molecular_multiplicity") == 1
            and _close(_get(molecule, "molecular_charge"), 0.0, 0.0)
            and _close_sequence(_get(molecule, "geometry"), _GEOMETRY_BOHR, 1.0e-7),
        )
    )
    checks.append(
        (
            "provenance ORCA cclib QCSchemaWriter",
            _get(provenance, "creator") == "ORCA"
            and "qcschemawriter" in str(_get(provenance, "routine")).lower(),
        )
    )
    checks.append(
        (
            "cclib-orca metadata",
            _get(extras, "cclib_harness", "selector") == "cclib-orca"
            and _get(extras, "cclib_harness", "parser") == "ORCA"
            and bool(_get(extras, "cclib_harness", "cclib_version"))
            and bool(_get(extras, "cclib_harness", "executable")),
        )
    )

    version = _get(provenance, "version")
    reference = _orca_reference(version)
    checks.append((f"ORCA extras reference version {version!r}", reference is not None))
    if reference is None:
        flat_ok = orbital_ok = scf_ok = False
    else:
        flat_ok = all(_matches(_get(extras, key), expected, 1.0e-6) for key, expected in reference["flat"].items())
        orbital_ok = all(
            _matches(_get(extras, key), expected, 1.0e-6) for key, expected in reference["orbital"].items()
        )
        scf_ok = all(_matches(_get(extras, key), expected, 1.0e-6) for key, expected in reference["scf"].items())
    checks.append(("representative flat CCSD extras", flat_ok))
    checks.append(("representative orbital extras", orbital_ok))
    checks.append(("representative SCF extras", scf_ok))
    checks.append(
        (
            "anomalous MP2 writer fields present but not numerically endorsed",
            _get(properties, "mp2_correlation_energy") is not None
            and _get(properties, "mp2_total_energy") is not None
            and "mpenergies" in extras,
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


def main(output_path: Any = "orca_water_ccsd.result.json") -> int:
    result = qcengine.compute(
        build_atomic_input(),
        "cclib-orca",
        raise_error=True,
        task_config=_TASK_CONFIG,
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
