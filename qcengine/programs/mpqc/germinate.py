"""Map QCSchema method names onto MPQC wavefunction types and property classes."""

import copy
from typing import Any, Dict, Optional, Tuple

from ...exceptions import InputError

ENERGY = "Energy"
EXCITATION_ENERGY = "ExcitationEnergy"

#: Generic method name -> (wfn type, needs separate ref block, DEFAULT property,
#: extra wfn keys). The property here is only the default; a `property__type`
#: keyword may select any other property the wfn type provides.
#: Extra keys are deep-copied per call so callers cannot corrupt the registry.
_METHOD_MAP: Dict[str, Tuple[str, bool, str, Dict[str, Any]]] = {
    "hf": ("SD", False, ENERGY, {}),
    "mp2": ("MP2", True, ENERGY, {"method": "standard"}),
    "ccsd": ("CCSD", True, ENERGY, {}),
    "ccsd(t)": ("CCSD(T)", True, ENERGY, {}),
    "ccsdt": ("CCSDT", True, ENERGY, {}),
    "cck": ("CCk", True, ENERGY, {"k": 2}),
    "sci": ("sCI", True, ENERGY, {}),
    "eom-ccsd": ("EOM-CCSD", True, EXCITATION_ENERGY, {}),
    "eom-ip-ccsd": ("EOM-IP-CCSD", True, EXCITATION_ENERGY, {}),
    "eom-ea-ccsd": ("EOM-EA-CCSD", True, EXCITATION_ENERGY, {}),
    "eom-cck": ("CCk", True, EXCITATION_ENERGY, {"k": 2, "eom": {"manifold": "2h2p"}}),
}

#: MPQC-native wfn type strings accepted verbatim. Real (non-complex,
#: non-periodic), spin-free LCAO wavefunctions only. This is an allow-list on
#: purpose.
_NATIVE_TYPES = {
    # SCF
    "RHF",
    "SD",
    "DF-RHF",
    "Direct-RHF",
    "DirectDF-RHF",
    "DFJ-RHF",
    "DFJ-CADFK-RHF",
    "CLR-DFJ-CADFK-RHF",
    # MP2 family
    "MP2",
    "RMP2F12",
    "DF-RMP2F12",
    "RLaplaceMP2",
    "PaoPnoRMP2",
    # CC family
    "CCSD",
    "CCSD(T)",
    "CCSDT",
    "CCSD(F12)",
    "CCSD(T)F12",
    "CCk",
    "CCk(F12)",
    # excited state / CI
    "CIS",
    "EOM-CCSD",
    "EOM-IP-CCSD",
    "EOM-EA-CCSD",
    "sCI",
    # Green's function
    "GF2F12",
}

#: Native types that are their own SCF and so take no `ref`.
_SELF_REFERENCING = {
    "RHF",
    "SD",
    "DF-RHF",
    "Direct-RHF",
    "DirectDF-RHF",
    "DFJ-RHF",
    "DFJ-CADFK-RHF",
    "CLR-DFJ-CADFK-RHF",
}

#: Native types that can produce excitation energies. A capability set, not a
#: classification. Membership is what a `property__type: ExcitationEnergy`
#: request is validated against.
_EXCITED_STATE_TYPES = {"EOM-CCSD", "EOM-IP-CCSD", "EOM-EA-CCSD", "sCI", "CIS", "CCk"}

#: Native types that provide only excitation energies.
#: These default to ExcitationEnergy and refuse an Energy request.
_EXCITATION_ONLY_TYPES = {"EOM-CCSD", "EOM-IP-CCSD", "EOM-EA-CCSD", "CIS"}

_SUPPORTED_PROPERTIES = {ENERGY.lower(): ENERGY, EXCITATION_ENERGY.lower(): EXCITATION_ENERGY}

_NATIVE_BY_LOWER = {name.lower(): name for name in _NATIVE_TYPES}

_SUPPORTED_DRIVERS = {"energy", "properties"}


def muster_modelchem(
    method: str, driver: str, property_type: Optional[str] = None
) -> Tuple[str, bool, str, Dict[str, Any]]:
    """Translate a QCSchema method and driver into MPQC input ingredients.

    Parameters
    ----------
    method
        ``AtomicInput.specification.model.method``. Case-insensitive. An
        optional ``mpqc-`` prefix is stripped.
    driver
        ``AtomicInput.specification.driver``, either the ``DriverEnum`` or a
        plain string.
    property_type
        The ``property__type`` keyword, case-insensitive, or None to take the
        method's default. The method name does not fix the property: it only
        supplies a default. Several MPQC wavefunctions provide more than one
        property class, so an explicit request is validated against what the
        resolved ``wfn_type`` actually provides.

    Returns
    -------
    wfn_type
        MPQC ``wfn.type``, in MPQC's exact casing.
    needs_ref
        True when a separate ``scf`` block must be emitted and referenced.
    property_type
        ``"Energy"`` or ``"ExcitationEnergy"``, in MPQC's exact casing.
    wfn_extras
        Extra keys the method requires inside the ``wfn`` block.

    Raises
    ------
    InputError
        For an unsupported driver, a method that is unknown or out of scope, or
        a property the resolved wavefunction type cannot provide.
    """
    # DriverEnum subclasses str but its __str__ returns "DriverEnum.energy",
    # so unwrap .value before lowercasing. Accepts a plain str too.
    driver = getattr(driver, "value", driver).lower()
    if driver not in _SUPPORTED_DRIVERS:
        raise InputError(
            f"{driver} not implemented for MPQC. MPQC has no analytic derivatives, and the "
            f"finite-difference FDGradient property is out of scope for this harness. "
            f"Supported drivers: {sorted(_SUPPORTED_DRIVERS)}."
        )

    normalized = method.lower().strip()
    if normalized.startswith("mpqc-"):
        normalized = normalized[len("mpqc-") :]

    if normalized in _METHOD_MAP:
        wfn_type, needs_ref, default_property, extras = _METHOD_MAP[normalized]
        wfn_extras = copy.deepcopy(extras)
    elif normalized in _NATIVE_BY_LOWER:
        wfn_type = _NATIVE_BY_LOWER[normalized]
        needs_ref = wfn_type not in _SELF_REFERENCING
        # Default to excitation energies only for types that provide nothing
        # else. sCI and CCk provide both, so they default to a total energy.
        default_property = EXCITATION_ENERGY if wfn_type in _EXCITATION_ONLY_TYPES else ENERGY
        wfn_extras = {}
    else:
        raise InputError(
            f"Method '{method}' not available in the MPQC harness. Accepted generic names: "
            f"{sorted(_METHOD_MAP)}. Accepted MPQC-native types: {sorted(_NATIVE_TYPES)}. "
            f"Complex/periodic (z*), gamma-point, MRA, GPU, and DFT methods are out of scope."
        )

    return wfn_type, needs_ref, _resolve_property(wfn_type, default_property, property_type), wfn_extras


def _resolve_property(wfn_type: str, default_property: str, requested: Optional[str]) -> str:
    """Validate a requested property against what ``wfn_type`` provides.

    Returns ``default_property`` when nothing was requested. Otherwise the
    request is checked against the wavefunction's capabilities rather than
    against the method name, so a dual-property type such as sCI or CCk can
    serve either an ``Energy`` or an ``ExcitationEnergy``.
    """
    if requested is None:
        return default_property

    resolved = _SUPPORTED_PROPERTIES.get(str(requested).lower())
    if resolved is None:
        raise InputError(
            f"MPQC property '{requested}' not supported by this harness. "
            f"Supported: {sorted(_SUPPORTED_PROPERTIES.values())}."
        )

    if resolved == EXCITATION_ENERGY and wfn_type not in _EXCITED_STATE_TYPES:
        raise InputError(
            f"MPQC wavefunction type '{wfn_type}' cannot provide ExcitationEnergy. "
            f"Types that can: {sorted(_EXCITED_STATE_TYPES)}."
        )
    if resolved == ENERGY and wfn_type in _EXCITATION_ONLY_TYPES:
        raise InputError(
            f"MPQC wavefunction type '{wfn_type}' cannot provide Energy; it computes "
            f"excitation energies only. Use the underlying ground-state method instead."
        )

    return resolved
