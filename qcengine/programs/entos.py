"""
Calls the entos executable.
"""

from typing import Any, Dict, Optional, List

from qcelemental.models import Result, FailedOperation
from ..util import execute, popen
from qcelemental.util import which, safe_version, parse_version
from ..exceptions import UnknownError

from .model import ProgramHarness
import string


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
        proc = self.execute(job_inputs)

        if isinstance(proc, FailedOperation):
            return proc
        else:
            # If execution succeeded, collect results
            result = self.parse_output(proc["outfiles"], input_data)
            return result

    def execute(self, inputs, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None):

        if extra_outfiles is not None:
            outfiles = ["dispatch.out"].extend(extra_outfiles)
        else:
            outfiles = ["dispatch.out"]

        if extra_commands is not None:
            commands = extra_commands
        else:
            commands = inputs["commands"]

        exe_success, proc = execute(commands,
                                    infiles=inputs["infiles"],
                                    outfiles=outfiles,
                                    scratch_directory=inputs["scratch_directory"],
                                    scratch_name=scratch_name,
                                    timeout=timeout
                                    )
        proc["outfiles"]["dispatch.out"] = proc["stdout"]

        # Determine whether the calculation succeeded
        if not exe_success:
            return UnknownError(proc["stderr"])

        return proc

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:

        # Write the geom xyz file with unit au
        xyz_file = input_model.molecule.to_string(dtype='xyz', units='Angstrom')

        # Create input dictionary
        if template is None:
            structure = {'structure': {'file': 'geometry.xyz'}}
            dft_info = {'xc': input_model.model.method,
                        'ao': input_model.model.basis.upper(),
                        'df_basis': input_model.keywords["df_basis"].upper(),
                        'charge': input_model.molecule.molecular_charge
                        }
            print_results = {'print': {'results': True}}

            if input_model.driver == 'energy':
                input_dict = {'dft': {**structure, **dft_info},
                              **print_results
                              }
            # Write gradient call if asked for
            elif input_model.driver == 'gradient':
                input_dict = {'gradient': {**structure, 'dft': {**dft_info}},
                              **print_results
                              }
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
            "commands": ["entos", "-n", str(config.ncores), "dispatch.in"],
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
