"""Run and validate the live cclib-backed ORCA water CCSD demonstration."""

import json
import math
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
_ORBITAL_ENERGIES = [-20.242269, -1.265784, -0.615347, -0.452279, -0.390877, 0.60058, 0.736585]
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

    ccenergies = _get(extras, "ccenergies")
    checks.append(
        (
            "representative flat CCSD extras",
            isinstance(ccenergies, (list, tuple))
            and len(ccenergies) >= 1
            and _close(ccenergies[-1], -75.013487814, 1.0e-6)
            and {"atomcharges", "atomcoords", "atomnos"} <= set(extras),
        )
    )
    moenergies = _get(extras, "moenergies")
    mosyms = _get(extras, "mosyms")
    checks.append(
        (
            "representative orbital extras",
            _get(extras, "homos") == [4]
            and isinstance(moenergies, (list, tuple))
            and len(moenergies) == 1
            and _close_sequence(moenergies[0], _ORBITAL_ENERGIES, 1.0e-3)
            and isinstance(mosyms, (list, tuple))
            and len(mosyms) == 1
            and len(mosyms[0]) == 7,
        )
    )
    scfvalues = _get(extras, "scfvalues")
    scftargets = _get(extras, "scftargets")
    scfenergies = _get(extras, "scfenergies")
    checks.append(
        (
            "representative SCF extras",
            isinstance(scfenergies, (list, tuple))
            and len(scfenergies) >= 1
            and _close(scfenergies[-1], -74.96357424008319, 1.0e-6)
            and isinstance(scftargets, (list, tuple))
            and bool(scftargets)
            and isinstance(scfvalues, (list, tuple))
            and bool(scfvalues)
            and isinstance(scfvalues[0], (list, tuple))
            and bool(scfvalues[0]),
        )
    )
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
