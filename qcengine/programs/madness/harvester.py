import re
import json
import logging
from decimal import Decimal
from typing import Tuple

import numpy as np
import qcelemental as qcel
from qcelemental.models import Molecule
from qcelemental.models.results import AtomicResultProperties
from qcelemental.molparse import regex

from ..util import PreservingDict

logger = logging.getLogger(__name__)


def harvest_output(outtext: str) -> Tuple[PreservingDict, Molecule, list, str, str]:
    """Function to read an entire NWChem output file.

    Reads all of the different "line search" segments of a file and returns
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
    # TODO (wardlt): Is it only necessary to read the last two steps?
    pass_psivar = []
    pass_coord = []
    pass_grad = []
    for outpass in re.split(r" Line search:", outtext, re.MULTILINE):
        psivar, madcoord, madgrad, version, error = harvest_outfile_pass(outpass)
        pass_psivar.append(psivar)## all the variables extracted
        pass_coord.append(madcoord)
        pass_grad.append(madgrad)

    # Determine which segment contained the last geometry
    retindx = -1 if pass_coord[-1] else -2

    return pass_psivar[retindx], pass_coord[retindx], pass_grad[retindx], version, error


def harvest_outfile_pass(outtext):
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
        r'^\s+' + r'MADNESS' + r'\s+' + r'(\d+.\d\d+.\d)' +r'\s'+ r'multiresolution suite'+r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        logger.debug('matched version')
        version = mobj.group(1)

    # Process SCF
    # 1)Fail to converge (TODO Robert ask for failed convergence)
    mobj = re.search(r'^\s+' + r'(?:Calculation failed to converge)' + r'\s*$', outtext, re.MULTILINE)
    if mobj:
        logger.debug('failed to converge')

    # 2)Calculation converged
    else:
        OPTIONS=[r'kinetic',r'nonlocal psp',r'nuclear attraction',r'coulomb',r'PCM',r'exchange-correlation',r'nuclear-repulsion',r'total']
        PSIVAR=['KINETIC ENERGY','NONLOCAL PSP','NUCLEAR ATTRACTION ENERGY','COULOMB','PCM','EXCHANGE-CORRELATION','NUCLEAR REPULSION ENERGY','TOTAL SCF ENERGY']
        optDict=dict(zip(OPTIONS,PSIVAR)) 

        for var,VAR in optDict.items():
            mobj = re.search(r'^\s+' + var + r'\s*' + NUMBER + r's*$', outtext, re.MULTILINE)
            if mobj:
                logger.debug('matched SCF')## not sure what this means
                psivar[VAR] = mobj.group(1)

         # 3) Initial geometry
#         mobj = re.search(
#             r'^\s+' + r'Geometry' + r'.*' + r'\s*' + r'^\s+' + r'(?:-+)\s*' + r'\s+' + r'\n' + r'^\s' +
#             r'Output coordinates in ' + r'(.*?)' + r'\s' + r'\(scale by' + r'.*' + r'\s' + r'to convert to a\.u\.\)' +
#             r'\s+' + r'\n' + r'^\s+' + r'No\.\       Tag          Charge          X              Y              Z' +
#             r'\s*' + r'^\s+' + r'---- ---------------- ---------- -------------- -------------- --------------' +
#             r'\s*' +
#             r'((?:\s+([1-9][0-9]*)+\s+([A-Z][a-z]*)+\s+\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)'
#             + r'\s*$', outtext, re.MULTILINE | re.IGNORECASE)

#         if mobj:
#             logger.debug('matched geom')

#             # dinky molecule w/ charge and multiplicity
#             if mobj.group(1) == 'angstroms':
#                 molxyz = '%d \n%d %d tag\n' % (len(mobj.group(2).splitlines()), out_charge, out_mult
#                                                )  # unit = angstrom
#                 for line in mobj.group(2).splitlines():
#                     lline = line.split()
#                     molxyz += '%s %16s %16s %16s\n' % (lline[-5], lline[-3], lline[-2], lline[-1])
#                     # Jiyoung was collecting charge (-4)? see if this is ok for ghosts
#                     # Tag    ,    X,        Y,        Z
#                 psivar_coord = Molecule(validate=False,
#                                         **qcel.molparse.to_schema(qcel.molparse.from_string(
#                                             molxyz, dtype='xyz+', fix_com=True, fix_orientation=True)["qm"],
#                                                                   dtype=2))

#             else:  # unit = a.u.
#                 molxyz = '%d au\n%d %d tag\n' % (len(mobj.group(2).splitlines()), out_charge, out_mult)
#                 for line in mobj.group(2).splitlines():
#                     lline = line.split()
#                     molxyz += '%s %16s %16s %16s\n' % (int(float(lline[-4])), lline[-3], lline[-2], lline[-1])
#                     # Tag    ,    X,        Y,        Z
#                 psivar_coord = Molecule(validate=False,
#                                         **qcel.molparse.to_schema(qcel.molparse.from_string(
#                                             molxyz, dtype='xyz+', fix_com=True, fix_orientation=True)["qm"],
#                                                                   dtype=2))

#         # Process gradient
#         mobj = re.search(
#             r'^\s+' + r'.*' + r'ENERGY GRADIENTS' + r'\s*' + r'\s+' + r'\n' + r'^\s+' +
#             r'atom               coordinates                        gradient' + r'\s*' + r'^\s+' +
#             r'x          y          z           x          y          z' + r'\s*' +
#             r'((?:\s+([1-9][0-9]*)+\s+([A-Z][a-x]*)+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)'
#             + r'\s*$', outtext, re.MULTILINE)

#         if mobj:
#             logger.debug('matched molgrad')
#             atoms = []
#             psivar_grad = []
#             for line in mobj.group(1).splitlines():
#                 lline = line.split()  # Num, Tag, coord x, coord y, coord z, grad x, grad y, grad z
#                 # print (lline)
#                 if lline == []:
#                     pass
#                 else:
#                     atoms.append(lline[1])  # Tag
#                     psivar_grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])

#         # Process dipole (Properties)
#         mobj = re.search(
#             r'^\s+' + r'Dipole moment' + r'\s+' + NUMBER + r'\s+' + r'A\.U\.' + r'\s*' + 
#             r'^\s+' + r'DMX' + r'\s+' + NUMBER + r'.*' +
#             r'^\s+' + r'DMY' + r'\s+' + NUMBER + r'.*' +
#             r'^\s+' + r'DMZ' + r'\s+' + NUMBER + r'.*' +
#             r'^\s+' + r'.*' +
#             r'^\s+' + r'Total dipole' + r'\s+' + NUMBER + r'\s+' + r'A\.U\.' + r'\s*' + 
#             r'^\s+' + r'Dipole moment' + r'\s+' + NUMBER + r'\s' + r'Debye\(s\)' + r'\s*' +
#             r'^\s+' + r'DMX' + r'\s+' + NUMBER + r'.*' +
#             r'^\s+' + r'DMY' + r'\s+' + NUMBER + r'.*' +
#             r'^\s+' + r'DMZ' + r'\s+' + NUMBER + r'.*' +
#             r'^\s+' + r'.*' + 
#             r'^\s+' + r'Total dipole' + r'\s+' + NUMBER + r'\s' + r'DEBYE\(S\)' + r'\s*$', 
#             outtext, re.MULTILINE)

#         if mobj:
#             logger.debug('matched total dipole')

#             # UNIT = DEBYE(S)
#             psivar['CURRENT DIPOLE X'] = mobj.group(7)
#             psivar['CURRENT DIPOLE Y'] = mobj.group(8)
#             psivar['CURRENT DIPOLE Z'] = mobj.group(9)
#             # total?

#             # Process error code
#             mobj = re.search(
#                 r'^\s+' + r'current input line \:' + r'\s*' + r'^\s+' + r'([1-9][0-9]*)' + r'\:' + r'\s+' + r'(.*)' +
#                 r'\s*' + r'^\s+'
#                 r'------------------------------------------------------------------------' + r'\s*' + r'^\s+'
#                 r'------------------------------------------------------------------------' + r'\s*' + r'^\s+' +
#                 r'There is an error in the input file' + r'\s*$', outtext, re.MULTILINE)
#             if mobj:
#                 logger.debug('matched error')
#             # print (mobj.group(1)) #error line number
#             # print (mobj.group(2)) #error reason
#             psivar['NWCHEM ERROR CODE'] = mobj.group(1)
#             # TODO process errors into error var

#     # fmt: on

#     # Get the size of the basis sets, etc
#     mobj = re.search(r"No. of atoms\s+:\s+(\d+)", outtext, re.MULTILINE)
#     if mobj:
#         psivar["N ATOMS"] = mobj.group(1)
#     mobj = re.search(
#         r"No. of electrons\s+:\s+(\d+)\s+Alpha electrons\s+:\s+(\d+)\s+Beta electrons\s+:\s+(\d+)",
#         outtext,
#         re.MULTILINE,
#     )
#     if mobj:
#         psivar["N ALPHA ELECTRONS"] = mobj.group(2)
#         psivar["N BETA ELECTRONS"] = mobj.group(3)
#         if psivar["N ALPHA ELECTRONS"] == psivar["N BETA ELECTRONS"]:
#             # get HOMO and LUMO energy
#             mobj = re.search(
#                 r"Vector"
#                 + r"\s+"
#                 + r"%d" % (psivar["N ALPHA ELECTRONS"])
#                 + r"\s+"
#                 + r"Occ="
#                 + r".*"
#                 + r"\s+"
#                 + r"E="
#                 + r"([+-]?\s?\d+[.]\d+)"
#                 + r"[D]"
#                 + r"([+-]0\d)",
#                 outtext,
#                 re.MULTILINE,
#             )
#             if mobj:
#                 homo = float(mobj.group(1)) * (10 ** (int(mobj.group(2))))
#                 psivar["HOMO"] = np.array([round(homo, 10)])
#             mobj = re.search(
#                 r"Vector"
#                 + r"\s+"
#                 + r"%d" % (psivar["N ALPHA ELECTRONS"] + 1)
#                 + r"\s+"
#                 + r"Occ="
#                 + r".*"
#                 + r"\s+"
#                 + r"E="
#                 + r"([+-]?\s?\d+[.]\d+)"
#                 + r"[D]"
#                 + r"([+-]0\d)",
#                 outtext,
#                 re.MULTILINE,
#             )
#             if mobj:
#                 lumo = float(mobj.group(1)) * (10 ** (int(mobj.group(2))))
#                 psivar["LUMO"] = np.array([round(lumo, 10)])

#     mobj = re.search(r"AO basis - number of functions:\s+(\d+)\s+number of shells:\s+(\d+)", outtext, re.MULTILINE)
#     if mobj:
#         psivar["N MO"] = mobj.group(2)
#         psivar["N BASIS"] = mobj.group(1)

#     # Search for Center of charge
#     mobj = re.search(
#         r"Center of charge \(in au\) is the expansion point"
#         + r"\n"
#         + r"\s+"
#         + r"X\s+=\s+([+-]?\d+[.]\d+)"
#         + r"\s+"
#         + r"Y\s+=\s+([+-]?\d+[.]\d+)"
#         + r"\s+"
#         + r"Z\s+=\s+([+-]?\d+[.]\d+)",
#         outtext,
#         re.MULTILINE,
#     )
#     if mobj:
#         psivar["CENTER OF CHARGE"] = np.array([mobj.group(1), mobj.group(2), mobj.group(3)])

#     mobj = re.search(
#         r"Dipole moment"
#         + r".*?"
#         + r"A\.U\."
#         + r"\s+"
#         + r"DMX\s+([+-]?\d+[.]\d+)\s+"
#         + r"DMXEFC\s+[+-]?\d+[.]\d+\s+"
#         + r"DMY\s+([+-]?\d+[.]\d+)\s+"
#         + r"DMYEFC\s+[+-]?\d+[.]\d+\s+"
#         + r"DMZ\s+([+-]?\d+[.]\d+)\s+"
#         + r"DMZEFC\s+[+-]?\d+[.]\d+\s+"
#         + r"\-EFC\-"
#         + r".*?"
#         + r"A\.U\.\s+"
#         + r"Total dipole\s+([+-]?\d+[.]\d+\s+)",
#         outtext,
#         re.MULTILINE,
#     )
#     # + r"DMY\s+" + r"([+-]?\d+[.]\d+)", outtext, re.MULTILINE)
#     if mobj:
#         psivar["DIPOLE MOMENT"] = np.array([mobj.group(1), mobj.group(2), mobj.group(3)])
#         psivar["TOTAL DIPOLE MOMENT"] = mobj.group(4)

#     mobj = re.search(
#         r"Quadrupole moments in atomic units\s+"
#         + r"Component\s+"
#         + r"Electronic\+nuclear\s+"
#         + r"Point charges\s+"
#         + r"Total\s+"
#         + r"-+\s+"
#         + r"XX\s+([+-]?\d+[.]\d+)\s+"
#         + r".*\s+.*\s+"
#         + r"YY\s+([+-]?\d+[.]\d+)\s+"
#         + r".*\s+.*\s+"
#         + r"ZZ\s+([+-]?\d+[.]\d+)\s+"
#         + r".*\s+.*\s+"
#         + r"XY\s+([+-]?\d+[.]\d+)\s+"
#         + r".*\s+.*\s+"
#         + r"XZ\s+([+-]?\d+[.]\d+)\s+"
#         + r".*\s+.*\s+"
#         + r"YZ\s+([+-]?\d+[.]\d+)\s+",
#         outtext,
#         re.MULTILINE,
#     )

#     if mobj:
#         psivar["QUADRUPOLE MOMENT"] = np.array(
#             [mobj.group(1), mobj.group(2), mobj.group(3), mobj.group(4), mobj.group(5), mobj.group(6)]
#         )

#     # Process CURRENT energies (TODO: needs better way)
#     if "HF TOTAL ENERGY" in psivar:
#         psivar["SCF TOTAL ENERGY"] = psivar["HF TOTAL ENERGY"]
#         psivar["CURRENT REFERENCE ENERGY"] = psivar["HF TOTAL ENERGY"]
#         psivar["CURRENT ENERGY"] = psivar["HF TOTAL ENERGY"]

         psivar["CURRENT EXCITATION ENERGY"] = psivar["%s EXCITATION ENERGY" % (cc_name)]

    return psivar, psivar_coord, psivar_grad, version, error


def harvest_hessian(hess: str) -> np.ndarray: pass 
    """Parses the contents of the NWChem hess file into a hessian array.

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


# def extract_formatted_properties(psivars: PreservingDict) -> AtomicResultProperties:
#     """Get named properties out of the general variables extracted out of the result file

#     Args:
#         psivars (PreservingDict): Dictionary of the output results
#     Returns:
#         (AtomicResultProperties) Properties in a standard format
#     """
#     # TODO (wardlt): Get more of the named variables out of the NWChem file

#     # Initialize the output
#     output = dict()

#     # Extract the Calc Info
#     output.update(
#         {
#             "calcinfo_nbasis": psivars.get("N BASIS", None),
#             "calcinfo_nmo": psivars.get("N MO", None),
#             "calcinfo_natom": psivars.get("N ATOMS", None),
#             "calcinfo_nalpha": psivars.get("N ALPHA ELECTRONS", None),
#             "calcinfo_nbeta": psivars.get("N BETA ELECTRONS", None),
#         }
#     )

#     # Get the "canonical" properties
#     output["return_energy"] = psivars["CURRENT ENERGY"]
#     output["nuclear_repulsion_energy"] = psivars["NUCLEAR REPULSION ENERGY"]

#     # Get the SCF properties
#     output["scf_total_energy"] = psivars.get("HF TOTAL ENERGY", None)
#     output["scf_one_electron_energy"] = psivars.get("ONE-ELECTRON ENERGY", None)
#     output["scf_two_electron_energy"] = psivars.get("TWO-ELECTRON ENERGY", None)
#     output["scf_dispersion_correction_energy"] = psivars.get("DFT DISPERSION ENERGY", None)

#     # Get the MP2 properties
#     output["mp2_total_correlation_energy"] = psivars.get("MP2 CORRELATION ENERGY", None)
#     output["mp2_total_energy"] = psivars.get("MP2 TOTAL ENERGY", None)
#     output["mp2_same_spin_correlation_energy"] = psivars.get("MP2 SAME-SPIN CORRELATION ENERGY", None)
#     output["mp2_opposite_spin_correlation_energy"] = psivars.get("MP2 OPPOSITE-SPIN CORRELATION ENERGY", None)
#     return AtomicResultProperties(**output)


def harvest(in_mol: Molecule, nwout: str, **outfiles) -> Tuple[PreservingDict, None, None, Molecule, str, str]:
    """Parses all the pieces of output from NWChem: the stdout in
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
    out_psivar, out_mol, out_grad, version, error = harvest_output(madout)

    # If available, read higher-accuracy gradients
    #  These were output using a Python Task in NWChem to read them out of the database
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
                    """Madness outfile (NRE: %f) inconsistent with Psi4 input (NRE: %f)."""
                    % (out_mol.nuclear_repulsion_energy(), in_mol.nuclear_repulsion_energy())
                )
    else:
        raise ValueError("""No coordinate information extracted from Madness output.""")

    # If present, align the gradients and hessian with the original molecular coordinates
    #  Madness rotates the coordinates of the input molecule. `out_mol` contains the coordinates for the
    #  rotated molecule, which we can use to determine how to rotate the gradients/hessian
    amol, data = out_mol.align(in_mol, atoms_map=False, mols_align=True, verbose=0)

    mill = data["mill"]  # Retrieve tool with alignment routines

    if out_grad is not None:
        out_grad = mill.align_gradient(np.array(out_grad).reshape(-1, 3))
    if out_hess is not None:
        out_hess = mill.align_hessian(np.array(out_hess))

    return out_psivar, out_hess, out_grad, out_mol, version, error
