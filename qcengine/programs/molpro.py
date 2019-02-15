"""
Calls the Molpro executable.
"""

from qcelemental.models import ComputeError, FailedOperation, Provenance, Result

from qcengine.units import ureg

import xml.etree.ElementTree as ET


def _format_input(input_model, config):
    input_file = []

    # Write header info
    input_file.append("!Title")
    memory_mw_core = int(config.memory * (1024**3) / 8e6 / config.ncores)

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
    posthf_methods = ['mp2', 'ccsd', 'ccsd(t)']
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


def _parse_output(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    #print(root.tag)

    # TODO Try to grab the last total energy in the general case?
    #      - Would be useful for arbitrarly complicated input file
    #      - However, would need every different string used specify that (e.g. HF --> Energy, MP2 --> total energy)
    # FIXME Need to properly construct output_data to pass to Result()
    output_data = {}

    # The jobstep tag in Molpro contains output from commands (e.g. {hf}, {force})
    for jobstep in root.findall(
            '{http://www.molpro.net/schema/molpro-output}job/{http://www.molpro.net/schema/molpro-output}jobstep'):
        #print("jobstep.tag: ")

        if 'SCF' in jobstep.attrib['command']:
            # Grab properties (e.g. Energy and Dipole moment)
            for child in jobstep.findall('{http://www.molpro.net/schema/molpro-output}property'):
                if child.attrib['name'] == 'Energy':
                    output_data['scf_method'] = child.attrib['method']
                    output_data['scf_total_energy'] = float(child.attrib['value'])

                elif child.attrib['name'] == 'Dipole moment':
                    output_data['scf_dipole_moment'] = [float(x) for x in child.attrib['value'].split()]

        elif 'MP2' in jobstep.attrib['command']:
            # Grab properties (e.g. Energy and Dipole moment)
            for child in jobstep.findall('{http://www.molpro.net/schema/molpro-output}property'):
                if child.attrib['name'] == 'total energy':
                    output_data['mp2_method'] = child.attrib['method']
                    output_data['mp2_total_energy'] = float(child.attrib['value'])

                elif child.attrib['name'] == 'correlation energy':
                    output_data['mp2_total_correlation_energy'] = float(child.attrib['value'])

                elif child.attrib['name'] == 'Dipole moment':
                    output_data['mp2_dipole_moment'] = [float(x) for x in child.attrib['value'].split()]

        # Grab gradient
        elif 'FORCE' in jobstep.attrib['command']:
            # Grab properties (e.g. Energy and Dipole moment)
            for child in jobstep.findall('{http://www.molpro.net/schema/molpro-output}gradient'):
                print("gradient.attrib: ")
                print(child.attrib)

    print(output_data)

    #return Result(**output_data)


def molpro(input_model, config):
    return _format_input(input_model, config)
