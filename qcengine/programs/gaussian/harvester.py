"""Parse Gaussian log output using cclib with regex fallback."""

import io
import re
from typing import Optional, Tuple

import numpy as np
import qcelemental as qcel
from qcelemental.models import Molecule

from ..util import PreservingDict

# rq-bb0a5d96 rq-c20746b1
_NORMAL_TERMINATION_RE = re.compile(r"Normal termination of Gaussian")
_VERSION_RE = re.compile(r"Gaussian\s+(\d+),\s+Revision\s+([\w.]+),")

# eV → Hartree conversion factor
# rq-bb0a5d96
_EV_TO_HARTREE = qcel.constants.conversion_factor("eV", "hartree")


# rq-d9434480
def harvest(
    in_mol: Molecule,
    method: str,
    log_content: str,
) -> Tuple[PreservingDict, Optional[np.ndarray], Optional[np.ndarray], Molecule]:
    """Parse a Gaussian log file and return QC data for building an AtomicResult.

    Parameters
    ----------
    in_mol:
        The input molecule used for geometry alignment and NRE cross-check.
    method:
        Lowercase QCSchema method string (e.g. "hf", "mp2", "b3lyp").
    log_content:
        Full text of gaussian.log.

    Returns
    -------
    (qcvars, grad, hess, out_mol)
        qcvars:  PreservingDict of QC variables (energies, charges, etc.).
        grad:    Gradient array shape (3*natoms,) in Hartree/Bohr, or None.
        hess:    Hessian matrix shape (3*natoms, 3*natoms) in Hartree/Bohr², or None.
        out_mol: Output molecule aligned to in_mol.

    Raises
    ------
    ValueError
        If no coordinate information can be extracted from the log.
    ValueError
        If the parsed geometry's NRE deviates from in_mol's by more than 1e-3 Hartree.
    """
    # rq-bb0a5d96 — parse with cclib; exceptions propagate to caller (rq-bb775969)
    import cclib

    ccdata = cclib.io.ccopen(io.StringIO(log_content)).parse()

    qcvars = PreservingDict()
    method_lower = method.lower()

    # --- Geometry and output molecule ---
    # rq-8d6fc024 rq-a76c3d7b rq-b7c7b8cb
    if not hasattr(ccdata, "atomcoords") or len(ccdata.atomcoords) == 0:
        raise ValueError("No coordinate information could be extracted from the Gaussian log.")

    # atomcoords in Ångström → convert to Bohr for QCElemental
    angstrom_to_bohr = qcel.constants.conversion_factor("angstrom", "bohr")
    coords_bohr = ccdata.atomcoords[-1] * angstrom_to_bohr  # shape (natoms, 3)

    symbols = [qcel.periodictable.to_symbol(z) for z in ccdata.atomnos]

    calc_mol = Molecule(
        symbols=symbols,
        geometry=coords_bohr.flatten(),
        fix_com=True,
        fix_orientation=True,
    )

    natoms = len(ccdata.atomnos)
    qcvars["N ATOMS"] = str(natoms)

    # rq-0753cd03
    qcvars["NUCLEAR REPULSION ENERGY"] = str(round(calc_mol.nuclear_repulsion_energy(), 10))

    # rq-90604835 — cross-check NRE against input molecule
    if abs(calc_mol.nuclear_repulsion_energy() - in_mol.nuclear_repulsion_energy()) > 1.0e-3:
        raise ValueError(
            f"Gaussian output geometry (NRE: {calc_mol.nuclear_repulsion_energy():.6f}) "
            f"is inconsistent with the input molecule (NRE: {in_mol.nuclear_repulsion_energy():.6f})."
        )

    # Align output molecule to input molecule (same pattern as GAMESS harvester)
    # rq-8d6fc024
    return_mol, _ = in_mol.align(calc_mol, atoms_map=False, mols_align=True, run_mirror=True, verbose=0)

    # --- Energies ---
    # rq-37f40e92 — HF/SCF energy is always extracted from scfenergies
    scf_hartree = float(ccdata.scfenergies[-1]) * _EV_TO_HARTREE
    qcvars["HF TOTAL ENERGY"] = str(scf_hartree)

    # rq-929a262a rq-c5165f1f — method-specific variables
    if method_lower in ("hf", "scf") or _is_dft(method_lower):
        qcvars["CURRENT ENERGY"] = str(scf_hartree)

    elif method_lower == "mp2":
        # rq-578bb785
        mp2_hartree = float(ccdata.mpenergies[-1][-1]) * _EV_TO_HARTREE
        qcvars["MP2 TOTAL ENERGY"] = str(mp2_hartree)
        qcvars["MP2 CORRELATION ENERGY"] = str(mp2_hartree - scf_hartree)
        qcvars["CURRENT ENERGY"] = str(mp2_hartree)

    elif method_lower == "ccsd":
        # rq-1b7ae62e
        ccsd_hartree = float(ccdata.ccenergies[-1]) * _EV_TO_HARTREE
        qcvars["CCSD TOTAL ENERGY"] = str(ccsd_hartree)
        qcvars["CCSD CORRELATION ENERGY"] = str(ccsd_hartree - scf_hartree)
        qcvars["CURRENT ENERGY"] = str(ccsd_hartree)

    elif method_lower == "ccsd(t)":
        # rq-eddef743
        # cclib stores CCSD(T) energy in ccenergies for Gaussian
        ccsdt_hartree = float(ccdata.ccenergies[-1]) * _EV_TO_HARTREE
        qcvars["CCSD(T) TOTAL ENERGY"] = str(ccsdt_hartree)
        qcvars["CCSD(T) CORRELATION ENERGY"] = str(ccsdt_hartree - scf_hartree)
        # Also populate CCSD energy if cclib exposes it separately
        # (cclib may not separate CCSD and CCSD(T) — use what we have)
        qcvars["CURRENT ENERGY"] = str(ccsdt_hartree)

    else:
        # Unknown method: fall back to SCF energy
        qcvars["CURRENT ENERGY"] = str(scf_hartree)

    # --- Gradient (forces → negate for gradient) ---
    # rq-36736115 rq-3b0adf96
    grad: Optional[np.ndarray] = None
    if hasattr(ccdata, "grads") and ccdata.grads is not None and len(ccdata.grads) > 0:
        # ccdata.grads shape: (nsteps, natoms, 3) in Hartree/Bohr (forces = -dE/dR)
        # Gradient = -forces
        grad = (-ccdata.grads[-1]).flatten()

    # --- Hessian ---
    # rq-f2e32162 rq-344f3432
    hess: Optional[np.ndarray] = None
    if hasattr(ccdata, "hessian") and ccdata.hessian is not None:
        hess = np.array(ccdata.hessian)

    # --- Properties: dipole and Mulliken charges ---
    # rq-c7fd08b5 rq-943856cc rq-53f8ed12 rq-4c4a1542 rq-ffae3fa6
    if hasattr(ccdata, "moments") and ccdata.moments is not None and len(ccdata.moments) > 1:
        dipole_vec = np.array(ccdata.moments[1])
        qcvars["DIPOLE MOMENT"] = str(float(np.linalg.norm(dipole_vec)))

    if hasattr(ccdata, "atomcharges") and ccdata.atomcharges:
        mulliken = ccdata.atomcharges.get("mulliken")
        if mulliken is not None:
            qcvars["MULLIKEN CHARGES"] = list(mulliken)

    return qcvars, grad, hess, return_mol


def is_normal_termination(log_content: str) -> bool:
    """Return True if the log contains a normal termination marker.

    rq-90141f40 rq-12c41d5e rq-c20746b1
    """
    return bool(_NORMAL_TERMINATION_RE.search(log_content))


def _is_dft(method_lower: str) -> bool:
    """Return True for method strings that are not recognised wavefunction methods."""
    return method_lower not in ("hf", "scf", "mp2", "ccsd", "ccsd(t)")
