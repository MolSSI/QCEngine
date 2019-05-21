"""
Calls the Molpro executable.
"""

import xml.etree.ElementTree as ET
from typing import Any, Dict, Optional

from qcelemental.models import Result, FailedOperation
import qcengine.util as uti
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
    version_cache: Dict[str, str] = {}

    class Config(ProgramExecutor.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    def found(self, raise_error: bool=False) -> bool:
        return which('molpro', return_bool=True, raise_error=raise_error, raise_msg='Please install via https://www.molpro.net/')

    def get_version(self) -> str:
        self.found(raise_error=True)

        #which_prog = which('molpro')
        #if which_prog not in self.version_cache:
        #    success, output = execute([which_prog, "v.inp"], {"v.inp": ""})

        #    if success:
        #        for line in output["stdout"].splitlines():
        #            if 'GAMESS VERSION' in line:
        #                branch = ' '.join(line.strip(' *\t').split()[3:])
        #        self.version_cache[which_prog] = safe_version(branch)

        #return self.version_cache[which_prog]
        
    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Run Molpro
        """
        # Check if Molpro executable is found
        self.found(raise_error=True)

        # Check Molpro version
        # if parse_version(self.get_version()) < parse_version("1.5"):
        #     raise TypeError("Molpro version '{}' not understood".format(self.get_version()))

        # Setup the job
        job_inputs = self.build_input(input_data, config)

        # Run Molpro
        exe_success, proc = uti.execute(job_inputs["commands"],
                                        infiles=job_inputs["infiles"],
                                        outfiles=["dispatch.out", "dispatch.xml"],
                                        scratch_location=job_inputs["scratch_location"],
                                        timeout=None
                                        )

        # Determine whether the calculation succeeded
        output_data = {}
        if not exe_success:
            output_data["success"] = False
            output_data["error"] = {"error_type": "internal_error",
                                    "error_message": proc["stderr"]
                                    }
            return FailedOperation(
                success=output_data.pop("success", False), error=output_data.pop("error"), input_data=output_data)

        # If execution succeeded, collect results
        result = self.parse_output(proc["outfiles"], input_data)
        return result

    # TODO Additional features
    #   - allow separate xyz file to be provided
    #   - pass scratch location info to Molpro command (-d option)
    #   - pass number of cores to Molpro command (-n option)
    #   - add support of wfu files? (-W option)
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

        # TODO Try to grab the last total energy in the general case?
        #      - Would be useful for arbitrarily complicated input file
        #      - However would need every different string used to specify energy (e.g. HF --> Energy, MP2 --> total energy)
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
        scf_extras = {}

        # MP2 maps
        mp2_energy_map = {
            "total energy": "mp2_total_energy",
            "correlation energy": "mp2_correlation_energy",
            "singlet pair energy": "mp2_singlet_pair_energy",
            "triplet pair energy": "mp2_triplet_pair_energy"
        }
        mp2_dipole_map = {"Dipole moment": "mp2_dipole_moment"}
        mp2_extras = {}

        # CCSD maps
        ccsd_energy_map = {
            "total energy": "ccsd_total_energy",
            "correlation energy": "ccsd_correlation_energy",
            "singlet pair energy": "ccsd_singlet_pair_energy",
            "triplet pair energy": "ccsd_triplet_pair_energy"
        }
        ccsd_dipole_map = {"Dipole moment": "ccsd_dipole_moment"}
        ccsd_extras = {}

        # Compiling the method maps
        scf_maps = {
            "energy": scf_energy_map,
            "dipole": scf_dipole_map,
            "extras": scf_extras
        }
        mp2_maps = {
            "energy": mp2_energy_map,
            "dipole": mp2_dipole_map,
            "extras": mp2_extras
        }
        ccsd_maps = {
            "energy": ccsd_energy_map,
            "dipole": ccsd_dipole_map,
            "extras": ccsd_extras
        }
        supported_methods = {"HF": scf_maps, "RHF": scf_maps, "MP2": mp2_maps, "CCSD": ccsd_maps}

        # The jobstep tag in Molpro contains output from commands (e.g. {hf}, {force})
        for jobstep in root.findall('molpro_uri:job/molpro_uri:jobstep', name_space):

            # Remove the -SCF part of the command string when Molpro calls HF or KS
            command = jobstep.attrib['command']
            if '-SCF' in command:
                command = command[:-4]

            # Grab energies and dipole moment
            if command in supported_methods:
                for child in jobstep.findall('molpro_uri:property', name_space):
                    if child.attrib['name'] in supported_methods[command]['energy']:
                        properties[supported_methods[command]['energy'][child.attrib['name']]] = float(child.attrib['value'])
                    elif child.attrib['name'] in supported_methods[command]['dipole']:
                        properties[supported_methods[command]['dipole'][child.attrib['name']]] = [float(x) for x in
                                                                                                  child.attrib['value'].split()]
            # Grab gradient
            elif 'FORCE' in jobstep.attrib['command']:
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
        method_energy_map = supported_methods[method]['energy']
        if 'total energy' in method_energy_map and method_energy_map['total energy'] in properties:
            final_energy = properties[method_energy_map['total energy']]
        elif 'Energy' in method_energy_map and method_energy_map['Energy'] in properties:
            final_energy = properties[method_energy_map['Energy']]
        else:
            raise KeyError("Could not find {:s} total energy".format(method))

        # Replace return_result with final_energy if gradient wasn't called
        if "return_result" not in output_data:
            output_data["return_result"] = final_energy

        output_data["properties"] = properties
        output_data['schema_name'] = 'qcschema_output'
        output_data['success'] = True

        return Result(**{**input_model.dict(), **output_data})
