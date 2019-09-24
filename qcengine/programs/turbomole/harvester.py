from decimal import Decimal
import re

from ..util import PreservingDict


def harvest(input_model, stdout, **outfiles):
    qcvars = PreservingDict()

    # Total energy from dscf or ridft
    total_energy_re = re.compile('total energy\s+=\s+([\d\-\.]+)')
    mobj = total_energy_re.search(stdout)
    total_energy = Decimal(mobj[1])

    # Check for DFT, default to HF
    energy_key = 'HF TOTAL ENERGY'
    dft_mobj = re.search('density functional', stdout)
    if dft_mobj:
        energy_key = 'DFT TOTAL ENERGY'
    qcvars[energy_key] = total_energy

    # Take into account energies from ricc2 runs. They will be different
    # from the HF energy.
    current_energy = total_energy
    qcvars['CURRENT ENERGY'] = current_energy

    gradient = None
    hessian = None
    import pdb; pdb.set_trace()
    return qcvars, gradient, hessian
