import re
import json
import logging

# from decimal import Decimal
from typing import Tuple

import numpy as np

# import qcelemental as qcel
from qcelemental.models import Molecule
from qcelemental.models.results import AtomicResultProperties
from qcelemental.molparse import regex

from ..util import PreservingDict

logger = logging.getLogger(__name__)


def harvest_moldft_output(outtext: str) -> Tuple[PreservingDict, Molecule, list, str, str]:
    """Function to read an entire MADNESS output file.

    Read all of the different "line search" segments of a file and returns
    values from the last segment for which a geometry was written.

    Args:
        outtext (str): Output written to stdout
    Returns:
        - (PreservingDict) Variables extracted from the output file in the last complete step
        - (Molecule): Molecule from the last complete step
        - (list): Gradient from the last complete step
        - (str): Version string
        - (str): Error message, if any
    """

    # Loop over all steps
    pass_psivar = []
    pass_coord = []
    pass_grad = []
    # Write now we split at Converge
    counter = 1

    splits = re.split(r"Converged!", outtext, re.MULTILINE)[-2]
    final_outpass = re.split(r"Iteration", splits, re.MULTILINE)[-1]
    psivar, madcoord, madgrad, version, error = harvest_outfile_moldft_pass(final_outpass)

    return psivar, madcoord, madgrad, version, error


def harvest_outfile_moldft_pass(outtext):
    """Function to read Madness output file *outtext* and parse important
    quantum chemical information from it in

    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None
    version = ""
    error = ""  # TODO (wardlt): The error string is never used.

    NUMBER = r"(?x:" + regex.NUMBER + ")"
    # fmt: off

    # Process version
    mobj = re.search(
        r'^\s+' + r'MADNESS' + r'\s+' + r'(\d+.\d\d+.\d)' + r'\s' + r'multiresolution suite' + r'\s*$',
        outtext, re.MULTILINE)
    # fmt: on
    if mobj:
        logger.debug("matched version")
        version = mobj.group(1)

    # Process SCF
    # 1)Fail to converge (TODO Robert ask for failed convergence)
    # fmt: off
    mobj = re.search(r'^\s+' + r'(?:Calculation failed to converge)' + r'\s*$', outtext, re.MULTILINE)
    # fmt: on
    if mobj:
        logger.debug("failed to converge")

    # 2)Calculation converged
    else:
        OPTIONS = [r"exchange-correlation", r"nuclear-repulsion", r"total"]
        PSIVAR = ["EXCHANGE-CORRELATION", "NUCLEAR REPULSION ENERGY", "TOTAL SCF ENERGY"]
        # OPTIONS=[r'kinetic',r'nonlocal psp',r'nuclear attraction',r'coulomb',r'PCM',r'exchange-correlation',r'nuclear-repulsion',r'total']
        # PSIVAR=['KINETIC ENERGY','NONLOCAL PSP','NUCLEAR ATTRACTION ENERGY','COULOMB','PCM','EXCHANGE-CORRELATION','NUCLEAR REPULSION ENERGY','TOTAL SCF ENERGY']
        optDict = dict(zip(OPTIONS, PSIVAR))

        for var, VAR in optDict.items():
            mobj = re.search(r"^\s+" + var + r"\s*" + NUMBER + r"s*$", outtext, re.MULTILINE)
            if mobj:
                logger.debug("matched SCF")  ## not sure what this means
                psivar[VAR] = mobj.group(1)
    # Other options

    # Process CURRENT energies (TODO: needs better way)
    if "TOTAL SCF ENERGY" in psivar:
        psivar["CURRENT REFERENCE ENERGY"] = psivar["TOTAL SCF ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["TOTAL SCF ENERGY"]

    return psivar, psivar_coord, psivar_grad, version, error


def harvest_hessian(hess: str) -> np.ndarray:
    pass


#    """Parses the contents of the NWChem hess file into a hessian array.

#     Args:
#         hess (str): Contents of the hess file
#     Returns:
#         (np.ndarray) Hessian matrix as a 2D array
#     """

# Change the "D[+-]" notation of Fortran output to "E[+-]" used by Python
#     hess_conv = hess.replace("D", "E")

#     # Parse all of the float values
#     hess_tri = [float(x) for x in hess_conv.strip().splitlines()]

#     # The value in the Hessian matrix is the lower triangle printed row-wise (e.g., 0,0 -> 1,0 -> 1,1 -> ...)
#     n = int(np.sqrt(8 * len(hess_tri) + 1) - 1) // 2  # Size of the 2D matrix

#     # Add the lower diagonal
#     hess_arr = np.zeros((n, n))
#     hess_arr[np.tril_indices(n)] = hess_tri

#     # Transpose and then set the lower diagonal again
#     hess_arr = np.transpose(hess_arr)  # Numpy implementations might only change the ordering to column-major
#     hess_arr[np.tril_indices(n)] = hess_tri

#     return hess_arr.T  # So that the array is listed in C-order, needed by some alignment routines


def extract_formatted_properties(psivars: PreservingDict) -> AtomicResultProperties:
    """Get named properties out of the general variables extracted out of the result file

    Args:
        psivars (PreservingDict): Dictionary of the output results
    Returns:
        (AtomicResultProperties) Properties in a standard format
    """
    # TODO (wardlt): Get more of the named variables out of the NWChem file

    # Initialize the output
    output = dict()

    # Extract the Calc Info
    output.update(
        {
            "calcinfo_nbasis": psivars.get("N BASIS", None),  ## Not a thing in madness
            "calcinfo_nmo": psivars.get("N MO", None),  ## Number of Mo orbitals
            "calcinfo_natom": psivars.get("N ATOMS", None),  ## Get madness to print this out
            "calcinfo_nalpha": psivars.get("N ALPHA ELECTRONS", None),  ## TODO (figure out how to read)
            "calcinfo_nbeta": psivars.get("N BETA ELECTRONS", None),
        }
    )

    # Get the "canonical" properties
    output["return_energy"] = psivars["CURRENT ENERGY"]
    output["nuclear_repulsion_energy"] = psivars["NUCLEAR REPULSION ENERGY"]

    # Get the SCF properties
    output["scf_total_energy"] = psivars.get("TOTAL SCF ENERGY", None)
    output["scf_xc_energy"] = psivars.get("EXCHANGE-CORRELATION", None)
    # TODO AdrianH right madness to output these variables
    # output["scf_one_electron_energy"] = psivars.get("ONE-ELECTRON ENERGY", None)
    # output["scf_two_electron_energy"] = psivars.get("TWO-ELECTRON ENERGY", None)
    # output["scf_dispersion_correction_energy"] = psivars.get("DFT DISPERSION ENERGY", None)

    return AtomicResultProperties(**output)


def harvest(in_mol: Molecule, **outfiles) -> Tuple[PreservingDict, None, None, Molecule, str, str]:
    """Parses all the pieces of output from Madness: the stdout in
    *nwout* Scratch files are not yet considered at this moment.

    Args:
        in_mol (Molecule): Input molecule
        madout (str): Madness output molecule
        outfiles (dict): Dictionary of outfile files and their contents
    Returns:
        - (PreservingDict) Variables extracted from the output file in the last complete step
        - (None): Hessian from the last complete step (Not yet implemented)
        - (None): Gradient from the last complete step (Not yet implemented)
        - (Molecule): Molecule from the last complete step
        - (str): Version string
        - (str): Error message, if any
    """

    # Parse the Madness output
    out_psivar, out_mol, out_grad, version, error = harvest_moldft_output(outfiles["moldft"]["stdout"])
    print(outfiles)
    if "molresponse" in outfiles.keys():
        response_psi_var= harvest_response_file(outfiles["molresponse"]["stdout"])
        out_psivar.update(response_psi_var)

    # If available, read higher-accuracy gradients
    #  These were output using a Python Task in Madness to read them out of the database
    if outfiles.get("mad.grad") is not None:
        logger.debug("Reading higher-accuracy gradients")
        out_grad = json.loads(outfiles.get("mad.grad"))

    # If available, read the hessian
    out_hess = None
    if outfiles.get("mad.hess") is not None:
        out_hess = harvest_hessian(outfiles.get("mad.hess"))

    # Make sure the input and output molecules are the same
    if out_mol:
        if in_mol:
            if abs(out_mol.nuclear_repulsion_energy() - in_mol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValueError(
                    """Madness outfile (NRE: %f) inconsistent with MADNESS input (NRE: %f)."""
                    % (out_mol.nuclear_repulsion_energy(), in_mol.nuclear_repulsion_energy())
                )
    # else:
    #    raise ValueError("""No coordinate information extracted from Madness output.""")

    # If present, align the gradients and hessian with the original molecular coordinates
    #  Madness rotates the coordinates of the input molecule. `out_mol` contains the coordinates for the
    #  rotated molecule, which we can use to determine how to rotate the gradients/hessian
    # amol, data = out_mol.align(in_mol, atoms_map=False, mols_align=True, verbose=0)

    # mill = data["mill"]  # Retrieve tool with alignment routines

    # if out_grad is not None:
    #    out_grad = mill.align_gradient(np.array(out_grad).reshape(-1, 3))
    # if out_hess is not None:
    #    out_hess = mill.align_hessian(np.array(out_hess))

    return out_psivar, out_hess, out_grad, out_mol, version, error


def harvest_response_file(outtext):
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None
    version = ""
    error = ""  # TODO (wardlt): The error string is never used.
    pass_psivar = []
    pass_coord = []
    pass_grad = []
    # Write now we split at Converge
    counter = 1

    splits = re.split(r"Converged!", outtext, re.MULTILINE)
    print(splits)
    splits = splits[-1]
    data = re.split(r"Iteration", splits, re.MULTILINE)[-1]
    print(data)

    NUMBER = r"(?x:" + regex.NUMBER + ")"  # NUMBER
    NUMSPACE = NUMBER + r"\s*"  # NUMBER + SPACE

    OPTIONS = [
        r"Number of Response States:",
        r"Number of Ground States:",
        r"k =",
    ]
    PSIVAR = ["NUM STATES", "NUM ORBITALS", "K"]
    optDict = dict(zip(OPTIONS, PSIVAR))

    for var, VAR in optDict.items():
        mobj = re.search(r"^\s*" + var + r"\s*" + NUMBER + r"\s*$", outtext, re.MULTILINE)
        # print(mobj)
        if mobj:
            psivar[VAR] = mobj.group(1)
    # Grab the Orbital Energies  There are NUM ORBITALS
    num_states = int(psivar["NUM STATES"])
    num_orbitals = int(psivar["NUM ORBITALS"])

    print(num_states)
    print(num_orbitals)
    # print(NUMSPACE)
    NUMSPACEORB = str()
    for i in range(num_orbitals):
        NUMSPACEORB += NUMSPACE
    # print(NUMSPACEORB)

    var = r"Orbital Energies: \[\*\]"
    VAR = "ORBITAL ENERGIES"
    mobj = re.search(
        r"^\s*" + var + r"\s*" + NUMSPACEORB + r"$",
        outtext,
        re.MULTILINE,
    )
    # print(mobj)

    if mobj:
        oe_list = []
        for i in range(num_orbitals):
            oe_list.append(mobj.group(i + 1))

        psivar[VAR] = np.array(oe_list, dtype=float)

    psivar = grab_tensor(r"Ground state overlap:", "OVERLAP", num_orbitals, num_orbitals, psivar, outtext)
    psivar = grab_tensor(r"Ground state hamiltonian:", "HAMILTONIAN", num_orbitals, num_orbitals, psivar, outtext)
    psivar = grab_tensor(r"Polarizability Final", "POLARIZABILITY", num_states, num_states, psivar, data)
    return psivar


def grab_tensor(var, VAR, row, col, psivar, data):
    first_line = r"^\s*" + var + r"\s+"
    NUMBER = r"(?x:" + regex.NUMBER + ")"  # NUMBER
    NUMSPACE = NUMBER + r"\s*"  # NUMBER + SPACE
    # print(first_line)

    CAPTURE_LINE = str()
    for j in range(col):
        CAPTURE_LINE += NUMSPACE
    total = first_line
    for i in range(row):
        front = r"^\[" + str(i) + r",\*\]\s*"
        line = front + CAPTURE_LINE
        total += line
    #    print(line)

    mobj = re.search(
        total,
        data,
        re.MULTILINE,
    )
    # print(mobj)
    if mobj:
        oe_list = []
        for i in range(row):
            for j in range(col):
                oe_list.append(mobj.group(i + 1))
        tensor = np.array(oe_list)
        psivar[VAR] = tensor.reshape((row, col))
    return psivar
