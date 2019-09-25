from decimal import Decimal
import re

from ..util import PreservingDict


def parse_decimal(regex, text, method="search"):
    with_float = re.compile(regex + "([\d\-\.]+)")
    matches = getattr(with_float, method)(text)

    if method == "search":
        matches = [mobj.groups()]
        # remaining_groups = list(mobj.groups())
        # energy = remaining_groups.pop()
        # return Decimal(energy), remaining_groups
    return [(method, Decimal(energy)) for method, energy in matches]


def parse_reference_energy(stdout):
    energy_dict = PreservingDict()

    # Total energy from dscf or ridft
    total_energy_re = re.compile('total energy\s+=\s+([\d\-\.]+)')
    mobj = total_energy_re.search(stdout)
    total_energy = Decimal(mobj[1])

    # Check for DFT, default to HF
    energy_key = 'HF TOTAL ENERGY'
    dft_mobj = re.search('density functional', stdout)
    if dft_mobj:
        energy_key = 'DFT TOTAL ENERGY'
    energy_dict[energy_key] = total_energy

    # Take into account energies from ricc2 runs. They will be different
    # from the HF energy.
    current_energy = total_energy
    energy_dict['CURRENT ENERGY'] = current_energy

    return energy_dict


def parse_ricc2(stdout):
    ricc2_dict = PreservingDict()

    # As CC2 starts from a MP2 guess that is also reported there may be
    # multiple matches for the following regex. Thats why we caputre all
    # matches with 'findall'.
    matches = parse_decimal("Final (.+?) energy\s+:\s+", stdout, "findall")
    ricc2_dict['CURRENT ENERGY'] = matches[-1][1]

    return ricc2_dict


def harvest(input_model, stdout, **outfiles):
    qcvars = PreservingDict()

    ref_energy_dict = parse_reference_energy(stdout)
    qcvars.update(ref_energy_dict)

    if "R I C C 2" in stdout:
        ricc2_dict = parse_ricc2(stdout)
        qcvars.update(ricc2_dict)

    gradient = None
    hessian = None
    return qcvars, gradient, hessian
