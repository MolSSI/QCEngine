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
        r'^\s+' + r'MADNESS' + r'\s+' +
        r'(\d+.\d\d+.\d)' + r'\s' + r'multiresolution suite' + r'\s*$',
        outtext, re.MULTILINE)
    # fmt: on
    if mobj:
        logger.debug("matched version")
        version = mobj.group(1)

    # Process SCF
    # 1)Fail to converge (TODO Robert ask for failed convergence)
    # fmt: off
    mobj = re.search(
        r'^\s+' + r'(?:Calculation failed to converge)' + r'\s*$', outtext, re.MULTILINE)
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
                logger.debug("matched SCF")  # not sure what this means
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


# gets calc info from psivars preservingDict
# before
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
    print("extract_formatted_properties", output)

    # Extract the Calc Info
    output.update(
        {
            # Not a thing in madness
            "calcinfo_nbasis": psivars.get("N BASIS", None),
            "calcinfo_nmo": psivars.get("N MO", None),  # Number of Mo orbitals
            # Get madness to print this out
            "calcinfo_natom": psivars.get("N ATOMS", None),
            # TODO (figure out how to read)
            "calcinfo_nalpha": psivars.get("N ALPHA ELECTRONS", None),
            "calcinfo_nbeta": psivars.get("N BETA ELECTRONS", None),
        }
    )

    # Get the "canonical" properties
    # output["return_energy"] = psivars["CURRENT ENERGY"]
    # output["nuclear_repulsion_energy"] = psivars["NUCLEAR REPULSION ENERGY"]

    # Get the SCF properties
    # output["scf_total_energy"] = psivars.get("TOTAL SCF ENERGY", None)
    # output["scf_xc_energy"] = psivars.get("EXCHANGE-CORRELATION", None)
    # TODO AdrianH right madness to output these variables
    # output["scf_one_electron_energy"] = psivars.get("ONE-ELECTRON ENERGY", None)
    # output["scf_two_electron_energy"] = psivars.get("TWO-ELECTRON ENERGY", None)
    # output["scf_dispersion_correction_energy"] = psivars.get("DFT DISPERSION ENERGY", None)

    return AtomicResultProperties(**output)


def harvest(in_mol: Molecule, outfiles) -> Tuple[PreservingDict, None, None, Molecule, str, str]:
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
    out_psivar = PreservingDict()
    # Parse the Madness output
    # This is a weird unpacking but i'm sure i'll find a more elegant way to do this later
    moldft_info = outfiles.get("moldft")
    moldft_outfiles = moldft_info.get("outfiles")
    # At this point scf prints a list of json outputs where each list refers to the scf at givin protocol
    # Here I load the scf_info and calc_info as json
    # scf_info = json.loads(moldft_outfiles.get("scf_info.json"))
    calc_info = json.loads(moldft_outfiles.get("calc_info.json"))
    # Write harvest scf_info and harvest calc_info
    out_calc_vars = harvest_calc_info(calc_info)
    # out_scf_vars = harvest_scf_info(scf_info)
    out_psivar.update(out_calc_vars)
    # out_psivar.update(out_scf_vars)

    if "molresponse" in outfiles.keys():
        molresponse_info = outfiles.get("moldft")
        molresponse_outfiles = molresponse_info.get("outfiles")
        # At this point scf prints a list of json outputs where each list refers to the scf at given protocol
        # Here I load the scf_info and calc_info as json
        response_info = json.loads(molresponse_outfiles.get("response_base.json"))
        response_params, response_data_dict = read_molrespone_json(response_info)

    Idontneed_vars, out_mol, out_grad, version, error = harvest_moldft_output(outfiles["moldft"]["stdout"])
    if "molresponse" in outfiles.keys():
        response_psi_var = harvest_response_file(outfiles["molresponse"]["stdout"])
        out_psivar.update(response_psi_var)
    # If available, read higher-accuracy gradients
    #  These were output using a Python Task in Madness to read them out of the database
    if outfiles.get("mad.grad") is not None:
        logger.debug("Reading higher-accuracy gradients")
        out_grad = json.loads(outfiles.get("mad.grad"))
    # If available, read the hessian
    # TODO read in the geometry outputs from a geometry optimization
    # out_mol = None
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
    # TODO create a madness json that outputs basic info like the version and github hash?

    return out_psivar, out_hess, out_grad, out_mol, version, error


# collect the scf_info json
# first iterate through single numbers
# Then iteration throught tensor values
# this format should work for response values in the future
def harvest_scf_info(scf_info):
    psivar = PreservingDict()

    # We only need to last set in the list
    scf_info = scf_info[-1][0]
    print("harvest scf info", scf_info)

    scf_number_vars = [
        "scf_one_electron_energy",
        "scf_two_electron_energy",
        "nuclear_repulsion_energy",
        "scf_vv10_energy",
        "scf_xc_energy",
        "scf_dispersion_correction_energy",
        "scf_total_energy",
        "scf_iterations",
    ]

    for var in scf_number_vars:
        if scf_info.get(var) is not None:
            psivar[var.upper()] = scf_info.get(var)

    scf_tensor_vars = ["scf_dipole_moment"]

    for var in scf_tensor_vars:
        if scf_info.get(var) is not None:
            psivar[var.upper()] = tensor_to_numpy(scf_info.get(var))

    return psivar


def harvest_calc_info(calc_info):
    psivar = PreservingDict()
    qcvars = ["calcinfo_nbasis", "calcinfo_nmo", "calcinfo_nalpha", "calcinfo_nbeta", "calcinfo_natom", "return_energy"]

    for var in qcvars:
        if calc_info.get(var) is not None:
            psivar[var.upper()] = calc_info.get(var)

    return psivar


def tensor_to_numpy(j):
    array = np.empty(j["size"])
    array[:] = j["vals"]
    print(tuple(j["dims"]))
    return np.reshape(array, tuple(j["dims"]))


def read_frequency_proto_iter_data(my_iter_data, num_states, num_orbitals):
    num_iters = len(my_iter_data)
    dres = np.empty((num_iters, num_states))
    res_X = np.empty((num_iters, num_states))
    res_Y = np.empty((num_iters, num_states))
    polar = np.empty((num_iters, 3, 3))
    for i in range(num_iters):
        dres[i, :] = tensor_to_numpy(my_iter_data[i]["density_residuals"])
        res_X[i, :] = tensor_to_numpy(my_iter_data[i]["res_X"])
        res_Y[i, :] = tensor_to_numpy(my_iter_data[i]["res_Y"])
        polar[i, :, :] = tensor_to_numpy(my_iter_data[i]["polar"])
    data = {}
    names = ["density_residuals", "res_X", "res_Y", "polar"]
    vals = [dres, res_X, res_Y, polar]
    for name, val in zip(names, vals):
        data[name] = val
    return data


def read_excited_proto_iter_data(my_iter_data, num_states, num_orbitals):
    num_iters = len(my_iter_data)
    dres = np.empty((num_iters, num_states))
    res_X = np.empty((num_iters, num_states))
    res_Y = np.empty((num_iters, num_states))
    omega = np.empty((num_iters, num_states))
    for i in range(num_iters):
        dres[i, :] = tensor_to_numpy(my_iter_data[i]["density_residuals"])
        res_X[i, :] = tensor_to_numpy(my_iter_data[i]["res_X"])
        res_Y[i, :] = tensor_to_numpy(my_iter_data[i]["res_Y"])
        omega[i, :] = tensor_to_numpy(my_iter_data[i]["omega"])
    data = {}
    names = ["density_residuals", "res_X", "res_Y", "omega"]
    vals = [dres, res_X, res_Y, omega]
    for name, val in zip(names, vals):
        data[name] = val
    return data


# input response_info json and returns a dict of response paramters
# and a list of dicts of numpy arrays holding response data
def read_molrespone_json(response_info):
    protocol_data = response_info["protocol_data"]
    response_parameters = response_info["response_parameters"]
    n_states = response_parameters["states"]
    n_orbitals = response_parameters["num_orbitals"]
    num_protos = len(protocol_data)
    protos = []
    proto_data = []
    for p in range(num_protos):
        protos.append(protocol_data[p]["proto"])
        iter_data = protocol_data[p]["iter_data"]
        if response_parameters["excited_state"]:
            proto_data.append(read_excited_proto_iter_data(iter_data, n_states, n_orbitals))
        else:
            proto_data.append(read_frequency_proto_iter_data(iter_data, n_states, n_orbitals))
    return response_parameters, proto_data


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


# Translate a madness tensor defined within json output to a numpy array


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
