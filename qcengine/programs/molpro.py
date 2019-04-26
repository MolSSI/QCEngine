"""
Calls the Molpro executable.
"""

import xml.etree.ElementTree as ET
from typing import Any, Dict, Optional

# from qcelemental.models import ComputeError, FailedOperation, Provenance, Result
from qcelemental.models import Result
from qcelemental.util import which

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

    # TODO Switch over to Jinja
    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:
        input_file = []
        posthf_methods = {'mp2', 'ccsd', 'ccsd(t)'}

        # Write header info
        input_file.append("!Title")
        memory_mw_core = int(config.memory * (1024 ** 3) / 8e6 / config.ncores)
        input_file.append("memory,{},M".format(memory_mw_core))
        input_file.append('')

        # Have no symmetry by default
        input_file.append('{symmetry,nosym}')
        input_file.append('')

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
            input_file.append('{HF}')
        input_file.append('{{{:s}}}'.format(input_model.model.method))

        # Write gradient call if asked for
        if input_model.driver == 'gradient':
            input_file.append('')
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

        # TODO Think of how to handle multiple calls of same command.
        #      Currently it will grab the last one in the file.
        # TODO Read information from molecule tag
        #      - cml:molecule, cml:atomArray (?)
        #      - basisSet
        #      - orbitals
        output_data = {}
        properties = {}
        name_space = {'molpro_uri': 'http://www.molpro.net/schema/molpro-output'}

        # SCF maps
        scf_energy_map = {"Energy": "scf_total_energy"}
        scf_dipole_map = {"Dipole moment": "scf_dipole_moment"}
        # scf_extras = {"method": "molpro_scf_method"}

        # MP2 maps
        mp2_energy_map = {
            "total energy": "mp2_total_energy",
            "correlation energy": "mp2_correlation_energy",
            "singlet pair energy": "mp2_singlet_pair_energy",
            "triplet pair energy": "mp2_triplet_pair_energy"
        }
        mp2_dipole_map = {"Dipole moment": "mp2_dipole_moment"}
        # mp2_extras = {}

        # CCSD maps
        ccsd_energy_map = {
            "total energy": "ccsd_total_energy",
            "correlation energy": "ccsd_correlation_energy",
            "singlet pair energy": "ccsd_singlet_pair_energy",
            "triplet pair energy": "ccsd_triplet_pair_energy"
        }
        ccsd_dipole_map = {"Dipole moment": "ccsd_dipole_moment"}
        # ccsd_extras = {}

        # Compiling the method maps
        supported_methods = {"HF", "RHF", "MP2", "CCSD"}
        energy_map = {
            "HF": scf_energy_map,
            "RHF": scf_energy_map,
            "MP2": mp2_energy_map,
            "CCSD": ccsd_energy_map
        }
        dipole_map = {
            "HF": scf_dipole_map,
            "RHF": scf_dipole_map,
            "MP2": mp2_dipole_map,
            "CCSD": ccsd_dipole_map
        }
        # extras_map = {
        #     "HF": scf_extras,
        #     "RHF": scf_extras,
        #     "MP2": mp2_extras,
        #     "CCSD": ccsd_extras
        # }

        # The jobstep tag in Molpro contains output from commands (e.g. {hf}, {force})
        for jobstep in root.findall('molpro_uri:job/molpro_uri:jobstep', name_space):

            # Remove the -SCF part of the command string when Molpro calls HF or KS
            command = jobstep.attrib['command']
            if '-SCF' in command:
                command = command[:-4]

            if command in supported_methods:
                for child in jobstep.findall('molpro_uri:property', name_space):
                    if child.attrib['name'] in energy_map[command]:
                        properties[energy_map[command][child.attrib['name']]] = float(child.attrib['value'])
                    elif child.attrib['name'] in dipole_map[command]:
                        properties[dipole_map[command][child.attrib['name']]] = [float(x) for x in
                                                                                 child.attrib['value'].split()]
            # Grab gradient
            elif 'FORCE' in jobstep.attrib['command']:
                # Grab properties (e.g. Energy and Dipole moment)
                for child in jobstep.findall('molpro_uri:gradient', name_space):
                    # Stores gradient as a single list where the ordering is [1x, 1y, 1z, 2x, 2y, 2z, ...]
                    output_data['return_result'] = [float(x) for x in child.text.split()]

        # Convert triplet and singlet pair correlation energies to opposite-spin and same-spin correlation energies
        if 'mp2_singlet_pair_energy' in properties and 'mp2_triplet_pair_energy' in properties:
            properties["mp2_same_spin_correlation_energy"] = (2.0 / 3.0) * properties['mp2_triplet_pair_energy']
            properties["mp2_opposite_spin_correlation_energy"] = (1.0 / 3.0) \
                                                                 * properties['mp2_triplet_pair_energy'] \
                                                                 + properties['mp2_singlet_pair_energy']
            del properties['mp2_singlet_pair_energy']
            del properties['mp2_triplet_pair_energy']

        if 'ccsd_singlet_pair_energy' in properties and 'ccsd_triplet_pair_energy' in properties:
            properties["ccsd_same_spin_correlation_energy"] = (2.0 / 3.0) * properties['ccsd_triplet_pair_energy']
            properties["ccsd_opposite_spin_correlation_energy"] = (1.0 / 3.0) \
                                                                  * properties['ccsd_triplet_pair_energy'] \
                                                                  + properties['ccsd_singlet_pair_energy']
            del properties['ccsd_singlet_pair_energy']
            del properties['ccsd_triplet_pair_energy']

        # A slightly more robust way of determining the final energy.
        # Throws an error if the energy isn't found for the method specified from the input_model.
        method = input_model.model.method
        if 'total energy' in energy_map[method]:
            if energy_map[method]['total energy'] in properties:
                final_energy = properties[energy_map[method]['total energy']]
            else:
                raise KeyError("Could not find {:s} total energy".format(method))
        elif 'Energy' in energy_map[method]:
            if energy_map[method]['Energy'] in properties:
                final_energy = properties[energy_map[method]['Energy']]
            else:
                raise KeyError("Could not find {:s} total energy".format(method))
        else:
            raise KeyError("Could not find {:s} total energy".format(method))

        # Replace return_result with final_energy if gradient wasn't called
        if "return_result" not in output_data:
            output_data["return_result"] = final_energy

        output_data["properties"] = properties
        output_data['schema_name'] = 'qcschema_output'
        output_data['success'] = True

        return Result(**{**input_model.dict(), **output_data})

    def found(self) -> bool:
        return which('molpro', return_bool=True)
