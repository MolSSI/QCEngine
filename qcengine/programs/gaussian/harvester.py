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

# eV → Hartree conversion factor (only used as fallback for cclib quantities
# that have no direct Hartree regex; direct log parsing is preferred for
# energies to avoid precision loss from cclib's slightly different eV factor).
# rq-bb0a5d96
_EV_TO_HARTREE = qcel.constants.conversion_factor("eV", "hartree")

# Regex patterns for parsing energies directly from Gaussian log in Hartree.
# These bypass cclib's eV round-trip which causes ~1e-6 precision loss due to
# a conversion-factor mismatch between cclib and qcelemental.
_SCF_DONE_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*([-\d.]+)")
_MP2_ENERGY_RE = re.compile(r"EUMP2\s*=\s*([-\d.DE+]+)")
_CCSD_T_ENERGY_RE = re.compile(r"CCSD\(T\)=\s*([-\d.DE+]+)")
_CCSD_ENERGY_RE = re.compile(r"DE\(Corr\)=\s*[-\d.]+\s+E\(CORR\)=\s*([-\d.]+)")


def _fortran_float(s: str) -> float:
    """Convert a Fortran-style D-exponent string (e.g. '-0.76D+02') to float."""
    return float(s.replace("D", "E").replace("d", "e"))


def _parse_force_constants(log_content: str, ndim: int) -> Optional[np.ndarray]:
    """Parse the lower-triangle 'Force constants in Cartesian coordinates' block
    from a Gaussian Freq log into a full symmetric (ndim, ndim) matrix.

    Returns None if the block is not found.
    """
    marker = "Force constants in Cartesian coordinates:"
    idx = log_content.find(marker)
    if idx == -1:
        return None

    # Extract the block: starts after marker, ends at "Leave Link"
    block_start = idx + len(marker)
    block_end = log_content.find("Leave Link", block_start)
    if block_end == -1:
        block_end = len(log_content)
    block = log_content[block_start:block_end]

    hess = np.zeros((ndim, ndim))
    # Parse lower-triangle entries: row lines have format "  <row>  <val1> <val2> ..."
    # Column header lines have format "       <col1>  <col2>  ..."
    col_offset = 0
    for line in block.splitlines():
        tokens = line.split()
        if not tokens:
            continue
        # Column header line: all tokens are integers
        try:
            cols = [int(t) for t in tokens]
            col_offset = cols[0] - 1  # 1-based → 0-based
            continue
        except ValueError:
            pass
        # Data line: first token is row index, rest are Fortran D-exponents
        try:
            row = int(tokens[0]) - 1  # 1-based → 0-based
            vals = [_fortran_float(t) for t in tokens[1:]]
        except (ValueError, IndexError):
            continue
        for j, val in enumerate(vals):
            c = col_offset + j
            hess[row, c] = val
            hess[c, row] = val  # symmetric

    return hess


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

    # Normalise method variants: strip Gaussian-specific sub-options like ",full"
    # so that energy-extraction logic sees the canonical method name.
    # e.g. "ccsd(t,full)" → "ccsd(t)"
    _METHOD_NORMALISE = {
        "ccsd(t,full)": "ccsd(t)",
    }
    method_lower = _METHOD_NORMALISE.get(method_lower, method_lower)

    # --- Geometry and output molecule ---
    # rq-8d6fc024 rq-a76c3d7b rq-b7c7b8cb
    if not hasattr(ccdata, "atomcoords") or len(ccdata.atomcoords) == 0:
        raise ValueError("No coordinate information could be extracted from the Gaussian log.")

    # atomcoords in Ångström → convert to Bohr for QCElemental
    angstrom_to_bohr = qcel.constants.conversion_factor("angstrom", "bohr")
    coords_bohr = ccdata.atomcoords[-1] * angstrom_to_bohr  # shape (natoms, 3)

    # rq-f1d712bb rq-66067019 rq-7411efa0 — cclib cannot distinguish ghost atoms
    # from real ones in Gaussian output (Gaussian marks ghosts via the "Atomic Type"
    # column = 1000, which cclib does not expose). Trust in_mol.real as the sole
    # source of ghost designation; sanity-check only that the parsed atom count
    # matches in_mol when ghosts are present.
    in_real = list(in_mol.real)
    has_ghosts = not all(in_real)
    if has_ghosts:
        if len(ccdata.atomnos) != len(in_mol.symbols):
            raise ValueError(
                f"Gaussian log atom count ({len(ccdata.atomnos)}) does not match the "
                f"input molecule's atom count ({len(in_mol.symbols)})."
            )
        symbols = list(in_mol.symbols)
        real_flags = in_real
    else:
        symbols = [qcel.periodictable.to_symbol(z) for z in ccdata.atomnos]
        real_flags = None

    # rq-f1d712bb — preserve ghost designation on calc_mol so its NRE excludes ghosts
    calc_mol_kwargs = dict(
        symbols=symbols,
        geometry=coords_bohr.flatten(),
        fix_com=True,
        fix_orientation=True,
    )
    if real_flags is not None:
        calc_mol_kwargs["real"] = real_flags
    calc_mol = Molecule(**calc_mol_kwargs)

    natoms = len(ccdata.atomnos)
    qcvars["N ATOMS"] = str(natoms)

    if hasattr(ccdata, "nbasis") and ccdata.nbasis is not None:
        qcvars["N BASIS FUNCTIONS"] = str(ccdata.nbasis)
    if hasattr(ccdata, "nmo") and ccdata.nmo is not None:
        qcvars["N MOLECULAR ORBITALS"] = str(ccdata.nmo)
    # cclib doesn't have nalpha/nbeta directly; derive from homos array
    if hasattr(ccdata, "homos") and ccdata.homos is not None:
        qcvars["N ALPHA ELECTRONS"] = str(ccdata.homos[0] + 1)
        nbeta = ccdata.homos[-1] + 1 if len(ccdata.homos) > 1 else ccdata.homos[0] + 1
        qcvars["N BETA ELECTRONS"] = str(nbeta)

    # rq-0753cd03
    qcvars["NUCLEAR REPULSION ENERGY"] = str(round(calc_mol.nuclear_repulsion_energy(), 10))

    # rq-90604835 rq-a69d0898 rq-e5a98119 — cross-check NRE against input molecule.
    # Molecule.nuclear_repulsion_energy() excludes ghost atoms (zero Z), so when
    # calc_mol.real mirrors in_mol.real both NREs sum over real atoms only.
    if abs(calc_mol.nuclear_repulsion_energy() - in_mol.nuclear_repulsion_energy()) > 1.0e-3:
        raise ValueError(
            f"Gaussian output geometry (NRE: {calc_mol.nuclear_repulsion_energy():.6f}) "
            f"is inconsistent with the input molecule (NRE: {in_mol.nuclear_repulsion_energy():.6f})."
        )

    # Frame considerations (same pattern as GAMESS harvester)
    # rq-8d6fc024 rq-f8ce079b
    align_kwargs = dict(atoms_map=False, mols_align=True, run_mirror=True, verbose=0)
    if has_ghosts:
        align_kwargs["generic_ghosts"] = True
    if in_mol.fix_com and in_mol.fix_orientation:
        # Impose input frame if important as signalled by fix_*=True
        return_mol = in_mol
        _, data = calc_mol.align(in_mol, **align_kwargs)
        mill = data["mill"]
    else:
        return_mol, _ = in_mol.align(calc_mol, **align_kwargs)
        mill = qcel.molutil.compute_scramble(
            len(in_mol.symbols), do_resort=False, do_shift=False, do_rotate=False, do_mirror=False
        )  # identity AlignmentMill

    # --- Energies ---
    # Parse energies directly from the Gaussian log in Hartree to avoid
    # precision loss from cclib's eV round-trip (cclib and qcelemental use
    # slightly different Hartree↔eV conversion factors, causing ~1e-6 error).
    # Fall back to cclib values if regex fails.

    # rq-37f40e92 — HF/SCF energy
    scf_match = _SCF_DONE_RE.search(log_content)
    if scf_match:
        scf_hartree = float(scf_match.group(1))
    else:
        scf_hartree = float(ccdata.scfenergies[-1]) * _EV_TO_HARTREE
    qcvars["HF TOTAL ENERGY"] = str(scf_hartree)
    qcvars["SCF TOTAL ENERGY"] = str(scf_hartree)
    qcvars["CURRENT REFERENCE ENERGY"] = str(scf_hartree)

    # rq-929a262a rq-c5165f1f — method-specific variables
    if method_lower in ("hf", "scf") or _is_dft(method_lower):
        qcvars["CURRENT ENERGY"] = str(scf_hartree)

    elif method_lower == "mp2":
        # rq-578bb785
        mp2_match = _MP2_ENERGY_RE.search(log_content)
        if mp2_match:
            mp2_hartree = _fortran_float(mp2_match.group(1))
        else:
            mp2_hartree = float(ccdata.mpenergies[-1][-1]) * _EV_TO_HARTREE
        qcvars["MP2 TOTAL ENERGY"] = str(mp2_hartree)
        qcvars["MP2 CORRELATION ENERGY"] = str(mp2_hartree - scf_hartree)
        qcvars["CURRENT ENERGY"] = str(mp2_hartree)

    elif method_lower == "ccsd":
        # rq-1b7ae62e
        # Find the last converged CCSD energy from the iteration lines
        ccsd_matches = list(_CCSD_ENERGY_RE.finditer(log_content))
        if ccsd_matches:
            ccsd_hartree = float(ccsd_matches[-1].group(1))
        else:
            ccsd_hartree = float(ccdata.ccenergies[-1]) * _EV_TO_HARTREE
        qcvars["CCSD TOTAL ENERGY"] = str(ccsd_hartree)
        qcvars["CCSD CORRELATION ENERGY"] = str(ccsd_hartree - scf_hartree)
        qcvars["CURRENT ENERGY"] = str(ccsd_hartree)

    elif method_lower == "ccsd(t)":
        # rq-eddef743
        ccsdt_match = _CCSD_T_ENERGY_RE.search(log_content)
        if ccsdt_match:
            ccsdt_hartree = _fortran_float(ccsdt_match.group(1))
        else:
            ccsdt_hartree = float(ccdata.ccenergies[-1]) * _EV_TO_HARTREE
        qcvars["CCSD(T) TOTAL ENERGY"] = str(ccsdt_hartree)
        qcvars["CCSD(T) CORRELATION ENERGY"] = str(ccsdt_hartree - scf_hartree)
        qcvars["CURRENT ENERGY"] = str(ccsdt_hartree)

    else:
        # Unknown method: fall back to SCF energy
        qcvars["CURRENT ENERGY"] = str(scf_hartree)

    # --- Gradient (forces → negate for gradient) ---
    # rq-36736115 rq-3b0adf96
    grad: Optional[np.ndarray] = None
    if hasattr(ccdata, "grads") and ccdata.grads is not None and len(ccdata.grads) > 0:
        # ccdata.grads shape: (nsteps, natoms, 3) in Hartree/Bohr (forces = -dE/dR)
        # Gradient = -forces; align_gradient expects (natom, 3)
        calc_grad = -ccdata.grads[-1]  # shape (natoms, 3)
        grad = mill.align_gradient(calc_grad).flatten()

    # --- Hessian ---
    # rq-f2e32162 rq-344f3432
    hess: Optional[np.ndarray] = None
    if hasattr(ccdata, "hessian") and ccdata.hessian is not None:
        # cclib hessian would be in standard orientation → transform with mill
        hess = mill.align_hessian(np.array(ccdata.hessian))
    else:
        # cclib does not parse force constants from Gaussian Freq output;
        # fall back to direct log parsing. Gaussian's "Force constants in
        # Cartesian coordinates" are in the INPUT orientation.
        parsed_hess = _parse_force_constants(log_content, 3 * natoms)
        if parsed_hess is not None:
            if in_mol.fix_com and in_mol.fix_orientation:
                # return_mol = in_mol (input frame), hessian already in input frame
                hess = parsed_hess
            else:
                # return_mol is in calc/standard frame; need to transform
                # hessian from input→standard. The mill from calc_mol.align(in_mol)
                # transforms calc→input, so we need input→calc.
                # rq-f8ce079b
                _, data = in_mol.align(calc_mol, **align_kwargs)
                input_to_calc_mill = data["mill"]
                hess = input_to_calc_mill.align_hessian(parsed_hess)

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
