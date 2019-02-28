"""
Calls the Molpro executable.
"""

import xml.etree.ElementTree as ET
from typing import Any, Dict, Optional

#from qcelemental.models import ComputeError, FailedOperation, Provenance, Result
from qcelemental.models import Result

from ..util import which
from .executor import ProgramExecutor


class MolproExecutor(ProgramExecutor):

    _defaults = {
        "name": "Molpro",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }

    class Config(ProgramExecutor.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        pass

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:
        input_file = []
        posthf_methods = {'mp2', 'ccsd', 'ccsd(t)'}

        # Write header info
        input_file.append("!Title")
        memory_mw_core = int(config.memory * (1024**3) / 8e6 / config.ncores)
        input_file.append("memory,{},M".format(memory_mw_core))
        input_file.append('')

        # TODO Hook for global options
        # input_file.append('{symmetry,nosym}')

        # Write the geom
        input_file.append('geometry={')
        for sym, geom in zip(input_model.molecule.symbols, input_model.molecule.geometry):
            s = "{:<4s} {:>{width}.{prec}f} {:>{width}.{prec}f} {:>{width}.{prec}f}".format(
                sym, *geom, width=14, prec=10)
            input_file.append(s)
        input_file.append('}')

        # Write charge and multiplicity
        input_file.append('set,charge={}'.format(input_model.molecule.molecular_charge))
        input_file.append('set,multiplicity={}'.format(input_model.molecule.molecular_multiplicity))
        input_file.append('')

        # Write the basis set
        input_file.append('basis={')
        input_file.append('default,{}'.format(input_model.model.basis))
        input_file.append('}')
        input_file.append('')

        # Start of Molpro Commands
        # Write energy call
        write_hf = input_model.model.method.lower() in posthf_methods
        if write_hf:
            input_file.append('{hf}')
        input_file.append('{{{:s}}}'.format(input_model.model.method))
        input_file.append('')

        # Write gradient call if asked for
        if input_model.driver == 'gradient':
            input_file.append('{force}')

        input_file = "\n".join(input_file)

        return {
            "commands": ["molpro", "dispatch.mol"],
            "infiles": {
                "dispatch.mol": input_file
            },
            "input_result": input_model.copy(deep=True)
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':
        tree = ET.ElementTree(ET.fromstring(outfiles["dispatch.xml"]))
        root = tree.getroot()
        # print(root.tag)

        # TODO Try to grab the last total energy in the general case?
        #      - Would be useful for arbitrarily complicated input file
        #      - However would need every different string used to specify energy (e.g. HF --> Energy, MP2 --> total energy)
        # TODO Read information from molecule tag
        #      - cml:molecule, cml:atomArray (?)
        #      - basisSet
        #      - orbitals
        output_data = {}
        name_space = {'molpro_uri': 'http://www.molpro.net/schema/molpro-output'}
        # TODO Create enum class to take jobstep.attrib['command'] and determine what method it is
        # methods_set = {'HF', 'RHF', 'MP2', 'CCSD'}

        mp2_map = {
            "total energy": "mp2_total_energy",
            "correlation energy": "mp2_total_correlation_energy",
            "singlet pair energy": "mp2_same_spin_correlation_energy",
            "triplet pair energy": "mp2_opposite_spin_correlation_energy",
        }

        properties = {}

        # The jobstep tag in Molpro contains output from commands (e.g. {hf}, {force})
        for jobstep in root.findall('molpro_uri:job/molpro_uri:jobstep', name_space):
            # print("jobstep.tag: ")
            # print(jobstep.tag)

            if 'SCF' in jobstep.attrib['command']:
                # Grab properties (e.g. Energy and Dipole moment)
                for child in jobstep.findall('molpro_uri:property', name_space):
                    if child.attrib['name'] == 'Energy':
                        # properties['scf_method'] = child.attrib['method']
                        properties['scf_total_energy'] = float(child.attrib['value'])

                    elif child.attrib['name'] == 'Dipole moment':
                        properties['scf_dipole_moment'] = [float(x) for x in child.attrib['value'].split()]

            elif 'MP2' in jobstep.attrib['command']:
                # Grab properties (e.g. Energy and Dipole moment)
                for child in jobstep.findall('molpro_uri:property', name_space):
                    if child.attrib['name'] in mp2_map:
                        properties[mp2_map[child.attrib['name']]] = float(child.attrib['value'])

            # Grab gradient
            # TODO Handle situation where there are multiple FORCE calls
            elif 'FORCE' in jobstep.attrib['command']:
                # Grab properties (e.g. Energy and Dipole moment)
                for child in jobstep.findall('molpro_uri:gradient', name_space):
                    # print("gradient.attrib: ")
                    # print(child.attrib)
                    # Stores gradient as a single list where the ordering is [1x, 1y, 1z, 2x, 2y, 2z, ...]
                    output_data['return_result'] = [float(x) for x in child.text.split()]

        # A _bad_ way of figuring the correct energy
        # TODO Maybe a better way would be to use the method specified from the input?
        if "return_result" not in output_data:
            if "mp2_total_energy" in properties:
                output_data["return_result"] = properties["mp2_total_energy"]
            elif "scf_total_energy" in properties:
                output_data["return_result"] = properties["scf_total_energy"]
            else:
                raise KeyError("Could not find SCF total energy")

        output_data["properties"] = properties
        output_data['schema_name'] = 'qcschema_output'

        # TODO Should only return True if Molpro calculation terminated properly
        output_data['success'] = True

        return Result(**{**input_model.dict(), **output_data})

    def found(self) -> bool:
        return which('dftd3', return_bool=True)
