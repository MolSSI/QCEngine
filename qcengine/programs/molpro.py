"""
Calls the Molpro executable.
"""

import string
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional

from qcelemental.models import Result
from qcelemental.util import parse_version, safe_version, which

from ..exceptions import InputError, UnknownError
from ..util import execute
from .model import ProgramHarness


class MolproHarness(ProgramHarness):
    _defaults = {
        "name": "Molpro",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    # Set of implemented dft functionals in Molpro according to
    # https://www.molpro.net/info/release/doc/manual/node208.html (version 2019.2)
    # TODO This isn't extensive because it definitely doesn't include component functionals (e.g. LYP, VWN5)
    _dft_functionals = {
        "B", "B-LYP", "BLYP", "B-P", "BP86", "B-VWN", "B3LYP", "B3LYP3", "B3LYP5", "B88X",
        "B97", "B97R", "BECKE", "BH-LYP", "CS", "D", "HFB", "HFS", "LDA", "LSDAC", "LSDC",
        "KYP88", "MM05", "MM05-2X", "MM06", "MM06-2X", "MM06-L", "MM06-HF", "PBE", "PBE0",
        "PBE0MOL", "PBEREV", "PW91", "S", "S-VWN", "SLATER", "VS99", "VWN", "VWN80", "M08-HX",
        "M08-SO", "M11-L", "TPSS", "TPSSH", "M12HFC", "HJSWPBE", "HJSWPBEH", "TCSWPBE", "PBESOL"
    }

    # Currently supported methods in QCEngine for Molpro
    _scf_methods = {"HF", "RHF", "KS", "RKS"}
    _post_hf_methods = {'MP2', 'CCSD', 'CCSD(T)'}
    _supported_methods = {*_scf_methods, *_post_hf_methods}

    class Config(ProgramHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which('molpro',
                     return_bool=True,
                     raise_error=raise_error,
                     raise_msg='Please install via https://www.molpro.net/')

    def get_version(self) -> str:
        self.found(raise_error=True)

        name_space = {'molpro_uri': 'http://www.molpro.net/schema/molpro-output'}
        which_prog = which('molpro')
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "version.inp", "-d", ".", "-W", "."],
                                      infiles={"version.inp": ""},
                                      outfiles=["version.out", "version.xml"])

            if success:
                tree = ET.ElementTree(ET.fromstring(output["outfiles"]["version.xml"]))
                root = tree.getroot()
                version_tree = root.find('molpro_uri:job/molpro_uri:platform/molpro_uri:version', name_space)
                year = version_tree.attrib['major']
                minor = version_tree.attrib['minor']
                molpro_version = year + "." + minor
                self.version_cache[which_prog] = safe_version(molpro_version)

        return self.version_cache[which_prog]

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Run Molpro
        """
        # Check if Molpro executable is found
        self.found(raise_error=True)

        # Check Molpro version
        if parse_version(self.get_version()) < parse_version("2018.1"):
            raise TypeError("Molpro version '{}' not supported".format(self.get_version()))

        # Setup the job
        job_inputs = self.build_input(input_data, config)

        # Run Molpro
        proc = self.execute(job_inputs)

        # Return proc if it is type UnknownError for error propagation otherwise process the output
        if isinstance(proc, UnknownError):
            return proc
        else:
            # If execution succeeded, collect results
            result = self.parse_output(proc["outfiles"], input_data)
            return result

    def execute(self,
                inputs: Dict[str, Any],
                extra_infiles: Optional[List[str]] = None,
                extra_outfiles: Optional[List[str]] = None,
                as_binary: Optional[List[str]] = None,
                extra_commands=None,
                scratch_name: Optional[str] = None,
                scratch_messy: bool = False,
                timeout: Optional[int] = None):
        """
        For option documentation go look at qcengine/util.execute
        """

        # Collect all input files and update with extra_infiles
        infiles = inputs["infiles"]
        if extra_infiles is not None:
            infiles.update(extra_infiles)

        # Collect all output files and update with extra_outfiles
        outfiles = ["dispatch.out", "dispatch.xml", "dispatch.wfu"]
        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        # Replace commands with extra_commands if present
        commands = inputs["commands"]
        if extra_commands is not None:
            commands = extra_commands

        # Run the Molpro program
        exe_success, proc = execute(commands,
                                    infiles=infiles,
                                    outfiles=outfiles,
                                    as_binary=as_binary,
                                    scratch_name=scratch_name,
                                    scratch_directory=inputs["scratch_directory"],
                                    scratch_messy=scratch_messy,
                                    timeout=timeout)

        # Determine whether the calculation succeeded
        if not exe_success:
            return UnknownError(proc["stderr"])

        return proc

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:
        if template is None:
            input_file = []

            # Memory is in megawords per core for Molpro
            memory_mw_core = int(config.memory * (1024**3) / 8e6 / config.ncores)
            input_file.append("memory,{},M".format(memory_mw_core))
            # TODO Decide how I want this keyword to look/work
            # if input_model.keywords.get('wfu') is not None:
            #     input_file.append('file,2,{}'.format(input_model.keywords.get('wfu')))
            input_file.append('')

            # Write the geom
            xyz_block = input_model.molecule.to_string(dtype='molpro', units='Bohr')
            input_file.append(xyz_block)

            # Write the basis set
            input_file.append('basis={')
            input_file.append('default,{}'.format(input_model.model.basis))
            input_file.append('}')
            input_file.append('')

            # Write energy call
            energy_call = []
            # If post-hf method is called then make sure to write a HF call first
            if input_model.model.method.upper() in self._post_hf_methods:
                energy_call.append('{HF}')
            # If DFT call make sure to write {rks,method}
            if input_model.model.method.upper() in self._dft_functionals:
                energy_call.append('{{rks,{:s}}}'.format(input_model.model.method))
            else:
                energy_call.append('{{{:s}}}'.format(input_model.model.method))

            # Write appropriate driver call
            if input_model.driver == 'energy':
                input_file.extend(energy_call)
            elif input_model.driver == 'gradient':
                input_file.extend(energy_call)
                input_file.append('')
                input_file.append('{force}')
            else:
                raise InputError('Driver {} not implemented for Molpro.'.format(input_model.driver))

            input_file = "\n".join(input_file)
        else:
            # Some of the potential different template options
            # (A) ordinary build_input (need to define a base template)
            # (B) user wants to add stuff after normal template (A)
            # (C) user knows their domain language (doesn't use any QCSchema quantities)

            # # Build dictionary for substitute
            # sub_dict = {
            #     "method": input_model.model.method,
            #     "basis": input_model.model.basis,
            #     "charge": input_model.molecule.molecular_charge
            # }

            # Perform substitution to create input file
            str_template = string.Template(template)
            input_file = str_template.substitute()

        return {
            "commands": ["molpro", "dispatch.mol", "-d", ".", "-W", ".", "-n",
                         str(config.ncores)],
            "infiles": {
                "dispatch.mol": input_file
            },
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True)
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':
        tree = ET.ElementTree(ET.fromstring(outfiles["dispatch.xml"]))
        root = tree.getroot()
        # print(root.tag)

        # TODO Read information from molecule tag
        #      - cml:molecule, cml:atomArray (?)
        #      - basisSet
        #      - orbitals
        output_data = {}
        properties = {}
        name_space = {'molpro_uri': 'http://www.molpro.net/schema/molpro-output'}

        # Molpro commands map
        molpro_map = {
            "Energy": {
                "HF": "scf_total_energy",
                "RHF": "scf_total_energy",
                "KS": "scf_total_energy",
                "RKS": "scf_total_energy"
            },
            "total energy": {
                "MP2": "mp2_total_energy",
                "CCSD": "ccsd_total_energy",
                "CCSD(T)": "ccsd_prt_pr_total_energy"
            },
            "correlation energy": {
                "MP2": "mp2_correlation_energy",
                "CCSD": "ccsd_correlation_energy",
                "CCSD(T)": "ccsd_prt_pr_correlation_energy",  # Need both CCSD(T) and Total
                "Total": "ccsd_prt_pr_correlation_energy"  # Total corresponds to CCSD(T) correlation energy
            },
            "singlet pair energy": {
                "MP2": "mp2_singlet_pair_energy",
                "CCSD": "ccsd_singlet_pair_energy"
            },
            "triplet pair energy": {
                "MP2": "mp2_triplet_pair_energy",
                "CCSD": "ccsd_triplet_pair_energy"
            },
            "Dipole moment": {
                "HF": "scf_dipole_moment",
                "RHF": "scf_dipole_moment",
                "KS": "scf_dipole_moment",
                "RKS": "scf_dipole_moment",
                "MP2": "mp2_dipole_moment",
                "CCSD": "ccsd_dipole_moment",
                "CCSD(T)": "ccsd_prt_pr_dipole_moment"
            }
        }

        # Molpro variables map used for quantities not found in the command map
        molpro_var_map = {
            "_ENUC": "nuclear_repulsion_energy",
            "_DFTFUN": "scf_xc_energy"
            # "_EMP2_SCS": "scs_mp2_total_energy"
        }

        # Loop through each jobstep
        # The jobstep tag in Molpro contains output from commands (e.g. {hf}, {force})
        for jobstep in root.findall('molpro_uri:job/molpro_uri:jobstep', name_space):
            # Remove the -SCF part of the command string when Molpro calls HF or KS
            command = jobstep.attrib['command']
            if '-SCF' in command:
                command = command[:-4]
            # Grab energies and dipole moment
            if command in self._supported_methods:
                for child in jobstep.findall('molpro_uri:property', name_space):
                    prop_name = child.attrib['name']
                    prop_method = child.attrib['method']
                    value = child.attrib['value']
                    if prop_name in molpro_map:
                        if prop_method in molpro_map[prop_name]:
                            if prop_name == "Dipole moment":
                                properties[molpro_map[prop_name][prop_method]] = [float(x) for x in value.split()]
                            else:
                                properties[molpro_map[prop_name][prop_method]] = float(value)
            # Grab gradient
            elif 'FORCE' in command:
                for child in jobstep.findall('molpro_uri:gradient', name_space):
                    # Stores gradient as a single list where the ordering is [1x, 1y, 1z, 2x, 2y, 2z, ...]
                    output_data['return_result'] = [float(x) for x in child.text.split()]

        # Look for final energy in the molecule tag in case it's needed
        # Note: For the DFT case mol_method is the name of the functional plus R or U in front
        molecule = root.find('molpro_uri:job/molpro_uri:molecule', name_space)
        mol_method = molecule.attrib['method']
        mol_final_energy = float(molecule.attrib['energy'])
        # Loop over each variable under the variables tag to grab additional info from molpro_var_map
        for variable in molecule.findall('molpro_uri:variables/molpro_uri:variable', name_space):
            var_name = variable.attrib['name']
            if var_name in molpro_var_map:
                print(var_name, variable[0].text)
                molpro_var_map[var_name] = float(variable[0].text)

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

        # Grab the method from input
        method = input_model.model.method
        if method.upper() in self._dft_functionals:  # Determine if method is a DFT functional
            method = "RKS"

        # A slightly more robust way of determining the final energy.
        # Throws an error if the energy isn't found for the method specified from the input_model.
        if method in self._post_hf_methods and molpro_map['total energy'][method] in properties:
            final_energy = properties[molpro_map['total energy'][method]]
        elif method in self._scf_methods and molpro_map['Energy'][method] in properties:
            final_energy = properties[molpro_map['Energy'][method]]
        else:
            # Back up method for determining final energy if not already present in properties
            # Use the total energy from the molecule tag if it matches the input method
            if input_model.model.method in mol_method:
                final_energy = mol_final_energy
                if method in self._post_hf_methods:
                    properties[molpro_map['total energy'][method]] = mol_final_energy
                    properties[
                        molpro_map['correlation energy'][method]] = mol_final_energy - properties['scf_total_energy']
                elif method in self._scf_methods:
                    properties[molpro_map['Energy'][method]] = mol_final_energy
            else:
                raise KeyError("Could not find {:s} total energy".format(method))

        # Replace return_result with final_energy if gradient wasn't called
        if "return_result" not in output_data:
            output_data["return_result"] = final_energy

        # Final output_data assignments needed for the Result object
        output_data["properties"] = properties
        output_data['schema_name'] = 'qcschema_output'
        output_data['stdout'] = outfiles["dispatch.out"]
        output_data['success'] = True

        return Result(**{**input_model.dict(), **output_data})
