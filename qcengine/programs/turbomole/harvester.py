from math import sqrt
import re
from decimal import Decimal

import numpy as np

from ..util import PreservingDict


def parse_decimal(regex, text, method="search"):
    with_float = re.compile(regex + r"([\d\-\.]+)")
    matches = getattr(with_float, method)(text)

    if method == "search":
        matches = [matches.groups()]
    return [(method, Decimal(energy)) for method, energy in matches]


def parse_reference_energy(stdout: str):
    """Parse stdout and return the energy of the reference wavefunction."""
    energy_dict = PreservingDict()

    # Total energy from dscf or ridft
    total_energy_re = re.compile(r"total energy\s+=\s+([\d\-\.]+)")
    mobj = total_energy_re.search(stdout)
    total_energy = Decimal(mobj[1])

    # Check for DFT, default to HF
    energy_key = "HF TOTAL ENERGY"
    dft_mobj = re.search("density functional", stdout)
    if dft_mobj:
        energy_key = "DFT TOTAL ENERGY"
    energy_dict[energy_key] = total_energy

    # Take into account energies from ricc2 runs. They will be different
    # from the HF energy.
    current_energy = total_energy
    energy_dict["CURRENT ENERGY"] = current_energy

    return energy_dict


def parse_ricc2(stdout: str):
    ricc2_dict = PreservingDict()

    # As CC2 starts from a MP2 guess that is also reported there may be
    # multiple matches for the following regex. Thats why we capture all
    # matches with 'findall'.
    matches = parse_decimal(r"Final (.+?) energy\s+:\s+", stdout, "findall")
    if len(matches) == 0:
        matches = parse_decimal(r"E(MP2)\s+:\s+", stdout, "search")

    ricc2_dict["CURRENT ENERGY"] = matches[-1][1]

    return ricc2_dict


def parse_gradient(gradient):
    grad_re = re.compile(
        r"\$grad.+"
        r"cycle =\s+(?P<cycle>\d+)\s+"
        r"(?P<energy_type>.+?) energy =\s+(?P<energy>[\d\.\-]+)\s+"
        r"\|dE/dxyz\| =\s+(?P<grad_norm>[\d\.]+)"
        r"(?P<coords_gradients>.+)\$end",
        re.DOTALL,
    )
    mobj = grad_re.match(gradient)

    # Commented out the variables that aren't returned so LGTM doesn't
    # complain.
    # energy_type = mobj.group("energy_type")
    # grad_norm = Decimal(mobj.group("grad_norm"))
    # energy = Decimal(mobj.group("energy"))
    coords_grad = mobj.group("coords_gradients")
    # cycle = int(mobj.group("cycle"))

    *_, grad = re.split("[a-z]{1,3}", coords_grad.strip())
    grad = np.array(grad.strip().replace("D", "E").split(), dtype=float)

    return grad


def parse_hessian(hessian):
    first_hessian = hessian.strip().split("$")[1]
    split = first_hessian.split()
    hess_type = split.pop(0)
    assert hess_type in ("nprhessian", "hessian")

    def is_float(str_):
        return "." in str_

    hess_items = [item for item in split if is_float(item)]
    coord_num = int(sqrt(len(hess_items)))
    assert coord_num**2 == len(hess_items)
    hessian = np.array(hess_items, dtype=float).reshape(-1, coord_num)

    return hessian


def harvest(input_model, stdout, **outfiles):
    qcvars = PreservingDict()

    ref_energy_dict = parse_reference_energy(stdout)
    qcvars.update(ref_energy_dict)

    if "R I C C 2" in stdout:
        ricc2_dict = parse_ricc2(stdout)
        qcvars.update(ricc2_dict)

    gradient = None
    if "gradient" in outfiles:
        gradient = parse_gradient(outfiles["gradient"])
        qcvars["N ATOMS"] = gradient.size // 3

    # Prefer unprojected 'nprhessian' over projected 'hessian'.
    hessian_text = outfiles.get("nprhessian", outfiles.get("hessian", None))

    if hessian_text is not None:
        hessian = parse_hessian(hessian_text)
        qcvars["N ATOMS"] = hessian.shape[0] // 3
    else:
        hessian = None

    return qcvars, gradient, hessian
