"""
Calls the Molpro executable.
"""

from qcelemental.models import ComputeError, FailedOperation, Provenance, Result

from qcengine.units import ureg


def _format_input(input_model, config):
    input_file = []

    # Write header info
    input_file.append("!Title")
    memory_mw_core = int(config.memory * (1024 ** 3) / 8e6 / config.ncores)

    input_file.append("memory,{},M".format(memory_mw_core))

    input_file.append('')

    # Write the geom
    input_file.append('geometry={')
    for sym, geom in zip(input_model.molecule.symbols, input_model.molecule.geometry):
        s = "{:<4s} {:>{width}.{prec}f} {:>{width}.{prec}f} {:>{width}.{prec}f}".format(sym, *geom, width=14, prec=10)
        input_file.append(s)
    input_file.append('}')

    # write charge and multiplicity
    input_file.append('set,charge={}'.format(input_model.molecule.molecular_charge))
    input_file.append('set,multiplicity={}'.format(input_model.molecule.molecular_multiplicity))

    input_file.append('')

    # write the basis set
    input_file.append('basis={')
    input_file.append('default,{}'.format(input_model.model.basis))
    input_file.append('}')

    input_file.append('')

    # Write Molpro commands
    posthf_methods = ['mp2','ccsd','ccsd(t)']
    write_hf = input_model.model.method.lower() in posthf_methods
    if write_hf:
        input_file.append('{hf}')
    input_file.append('{{{:s}}}'.format(input_model.model.method))

    input_file.append('')

    # Write gradient call if asked for
    if input_model.driver == 'gradient':
        input_file.append('{force}')

    input_file = "\n".join(input_file)
    print(input_file)


def molpro(input_model, config):
    return _format_input(input_model, config)
