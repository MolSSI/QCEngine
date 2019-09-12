"""
Calls the entos executable.
"""

import string
import json
from typing import Any, Dict, List, Optional, Tuple

from qcelemental.models import Result
from qcelemental.util import parse_version, safe_version, which

from ..exceptions import InputError, UnknownError
from ..util import execute, popen
from .model import ProgramHarness


class EntosHarness(ProgramHarness):
    _defaults = {
        "name": "entos",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    # Energy commands that are currently supported
    _energy_commands = {"dft", "xtb"}

    # List of DFT functionals
    _dft_functionals = {
        "SLATER", "DIRAC", "SLATERD3", "DIRACD3", "VWN5", "VWN", "VWN1", "SVWN", "LDA", "BLYP", "BPW91",
        "BLYPD3", "B88", "PBEX", "PBERX", "PBEC", "LYP", "PW91", "P86", "PBE", "PBER", "PBED3", "PBERD3",
        "B3LYP3", "B3LYP", "B3LYP5", "PBE0", "PBE1PBE", "B3LYP3D3", "B3LYPD3", "B3LYP5D3", "PBE0D3",
        "PBE1PBED3", "CAMB3LYP", "WB97X", "CAMB3LYPD3", "WB97XD3"
    }

    # Available keywords for each of the energy commands
    _dft_keywords = {"df"}
    _xtb_keywords = {}
    _keyword_map = {
        "dft": _dft_keywords,
        "xtb": _xtb_keywords
    }

    class Config(ProgramHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which('entos', return_bool=True, raise_error=raise_error, raise_msg='Please install via XXX')

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which('entos')
        if which_prog not in self.version_cache:
            with popen([which_prog, '--version']) as exc:
                exc["proc"].wait(timeout=15)
            self.version_cache[which_prog] = safe_version(exc["stdout"].split()[2])

        return self.version_cache[which_prog]

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Run entos
        """
        # Check if entos executable is found
        self.found(raise_error=True)

        # Check entos version
        if parse_version(self.get_version()) < parse_version("0.5"):
            raise TypeError("entos version '{}' not supported".format(self.get_version()))

        # Setup the job
        job_inputs = self.build_input(input_data, config)

        # Run entos
        exe_success, proc = self.execute(job_inputs)

        # Determine whether the calculation succeeded
        if exe_success:
            # If execution succeeded, collect results
            result = self.parse_output(proc["outfiles"], input_data)
            return result
        else:
            # Return UnknownError for error propagation
            return UnknownError(proc["stderr"])

    def execute(self,
                inputs: Dict[str, Any],
                extra_infiles: Optional[Dict[str, str]] = None,
                extra_outfiles: Optional[List[str]] = None,
                extra_commands: Optional[List[str]] = None,
                scratch_name: Optional[str] = None,
                scratch_messy: bool = False,
                timeout: Optional[int] = None) -> Tuple[bool, Dict[str, Any]]:
        """
        For option documentation go look at qcengine/util.execute
        """

        # Collect all input files and update with extra_infiles
        infiles = inputs["infiles"]
        if extra_infiles is not None:
            infiles.update(extra_infiles)

        # Collect all output files and extend with with extra_outfiles
        outfiles = ["dispatch.out", "results.json"]
        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        # Replace commands with extra_commands if present
        commands = inputs["commands"]
        if extra_commands is not None:
            commands = extra_commands

        # Run the entos program
        exe_success, proc = execute(commands,
                                    infiles=infiles,
                                    outfiles=outfiles,
                                    scratch_name=scratch_name,
                                    scratch_directory=inputs["scratch_directory"],
                                    scratch_messy=scratch_messy,
                                    timeout=timeout)

        # Entos does not create an output file and only prints to stdout
        proc["outfiles"]["results.json"] = proc["stdout"]
        return exe_success, proc

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:

        # Write the geom xyz file with unit au
        xyz_file = input_model.molecule.to_string(dtype='xyz', units='Angstrom')

        # Create input file
        if template is None:

            # Determine the energy_command
            if input_model.model.method.upper() in self._dft_functionals:
                energy_command = "dft"
            # For now entos supports HF calculations through the dft energy_command with xc = HF
            elif input_model.model.method.upper() == "HF":
                energy_command = "dft"
            else:
                energy_command = input_model.model.method.lower()

            # Check method is supported
            if energy_command not in self._energy_commands:
                raise InputError(f'Energy method, {input_model.model.method}, not implemented for entos.')

            # Define base options for the energy command (options that can be taken directly from input_model)
            energy_options = {
                "dft": {
                    'xc': input_model.model.method.upper(),
                    'ao': input_model.model.basis,
                    'charge': input_model.molecule.molecular_charge,
                    'spin': float(input_model.molecule.molecular_multiplicity - 1),
                },
                "xtb": {}
            }

            # Resolve keywords (extra options) for the energy command
            caseless_keywords = {k.lower(): v for k, v in input_model.keywords.items()}
            energy_extra_options = {}
            for key in caseless_keywords.keys():
                if key in self._keyword_map[energy_command]:
                    energy_extra_options[key] = caseless_keywords[key]

            # Additional sub trees
            structure = {'structure': {'file': 'geometry.xyz'}}  # Structure sub tree
            print_results = {'print': {'results': True}}  # Print sub tree
            name_results = {'name': 'json_results'}

            # Create the input dictionary for a energy call
            if input_model.driver == 'energy':
                input_dict = {
                    energy_command: {
                        **structure, **energy_options[energy_command], **energy_extra_options, **name_results
                    },
                    **print_results
                }
            # Create the input dictionary for a gradient call
            elif input_model.driver == 'gradient':
                input_dict = {
                    'gradient': {
                        **structure,
                        energy_command: {
                            **energy_options[energy_command], **energy_extra_options
                        },
                        **name_results
                    },
                    **print_results
                }
            # TODO Add support for hessians
            # elif input_model.driver == 'hessian':
            else:
                raise NotImplementedError('Driver {} not implemented for entos.'.format(input_model.driver))

            # Write input file
            input_file = self.write_input_recursive(input_dict)
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
            #     "df_basis": input_model.keywords["df_basis"].upper(),
            #     "charge": input_model.molecule.molecular_charge
            # }

            # Perform substitution to create input file
            str_template = string.Template(template)
            input_file = str_template.substitute()

        return {
            "commands": ["entos", "-n", str(config.ncores), "-o", "dispatch.out", "--json-results", "dispatch.in"],
            "infiles": {
                "dispatch.in": input_file,
                "geometry.xyz": xyz_file
            },
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True)
        }

    def write_input_recursive(self, d: Dict[str, Any]) -> List:
        input_file = []
        for key, value in d.items():
            if isinstance(value, dict):
                input_file.append(key + '(')
                rec_input = self.write_input_recursive(value)
                indented_line = map(lambda x: "  " + x, rec_input)
                input_file.extend(indented_line)
                input_file.append(')')
            else:
                if isinstance(value, str):
                    input_file.append("{0} = '{1}'".format(key, value))
                elif isinstance(value, bool):
                    input_file.append("{0} = {1}".format(key, str(value).lower()))
                else:
                    input_file.append("{0} = {1}".format(key, value))
        return input_file

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':

        output_data = {}
        properties = {}

        xtb_map = {}
        dft_map = {"energy", "n_iter", "ao_basis", "density", "orbitals", "fock"}
        gradient_map = {"energy", "gradient"}

        # TODO Parse the results.json outfile instead of text parsing
        #      One issue is that depending on driver different levels of information is returned.
        # tmp = json.load(outfiles["results.json"])
        # entos_results = tmp["json_results"]

        # Parse the output file, collect properties and gradient
        output_lines = outfiles["dispatch.out"].split('\n')
        gradients = []
        natom = len(input_model.molecule.symbols)
        for idx, line in enumerate(output_lines):
            fields = line.split()
            if fields[:1] == ["energy:"]:
                properties["scf_total_energy"] = float(fields[-1])
            elif fields[:2] == ["Molecular", "Dipole:"]:
                properties["scf_dipole_moment"] = [float(x) for x in fields[2:5]]
            elif fields[:3] == ["SCF", "converged", "in"]:
                properties["scf_iterations"] = int(fields[3])
            elif fields == ["Gradient", "(hartree/bohr):"]:
                # Gradient is stored as (dE/dx1,dE/dy1,dE/dz1,dE/dx2,dE/dy2,...)
                for i in range(idx + 2, idx + 2 + natom):
                    grad = output_lines[i].strip('\n').split()[1:]
                    gradients.extend([float(x) for x in grad])

        if input_model.driver == 'gradient':
            if len(gradients) == 0:
                raise ValueError('Gradient not found.')
            else:
                output_data["return_result"] = gradients

        # Replace return_result with final_energy if gradient wasn't called
        if "return_result" not in output_data:
            if "scf_total_energy" in properties:
                output_data["return_result"] = properties["scf_total_energy"]
            else:
                raise KeyError("Could not find SCF total energy")

        output_data["properties"] = properties
        output_data['schema_name'] = 'qcschema_output'
        output_data['success'] = True

        return Result(**{**input_model.dict(), **output_data})
