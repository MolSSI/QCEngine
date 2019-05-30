"""
Calls the entos executable.
"""

from typing import Any, Dict, Optional, List

from qcelemental.models import Result, FailedOperation
from ..util import execute, popen
from qcelemental.util import which, safe_version, parse_version

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

    class Config(ProgramHarness.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    def found(self, raise_error: bool = False) -> bool:
        return which('entos', return_bool=True, raise_error=raise_error, raise_msg='Please install via XXX')

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which('entos')
        if which_prog not in self.version_cache:
            with popen([which_prog, '--version']) as exc:
                exc["proc"].wait(timeout=5)
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
        exe_success, proc = execute(job_inputs["commands"],
                                    infiles=job_inputs["infiles"],
                                    outfiles=["dispatch.out"],
                                    scratch_location=job_inputs["scratch_location"],
                                    timeout=None
                                    )
        proc["outfiles"]["dispatch.out"] = proc["stdout"]

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

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:

        # Write the geom xyz file with unit au
        xyz_file = input_model.molecule.to_string(dtype='xyz', units='Angstrom')

        # Write input file
        input_dict = {}

        # TODO Is there a way to require keywords?
        if input_model.driver == 'energy':
            input_dict = {'dft': {
                            'structure': {
                              'file': 'geometry.xyz'
                              },
                            'xc': input_model.model.method,
                            'ao': input_model.model.basis.upper(),
                            'df_basis': input_model.keywords["df_basis"].upper(),
                            'charge': input_model.molecule.molecular_charge
                            }
                          }

        # Write gradient call if asked for
        if input_model.driver == 'gradient':
            input_dict = {'gradient': {
                            'structure': {
                              'file': 'geometry.xyz'
                              },
                            'dft': {
                              'xc': input_model.model.method,
                              'ao': input_model.model.basis,
                              'df_basis': input_model.keywords["df_basis"].upper(),
                              'charge': input_model.molecule.molecular_charge
                               }
                             }
                          }

        # Write input file
        input_file = self.write_input_recursive(input_dict)
        input_file = "\n".join(input_file)

        return {
            "commands": ["entos", "-n", str(config.ncores), "dispatch.in"],
            "infiles": {
                "dispatch.in": input_file,
                "geometry.xyz": xyz_file
            },
            "scratch_location": config.scratch_directory,
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
            if fields[:2] == ["TOTAL", "ENERGY:"]:
                properties["scf_total_energy"] = float(fields[-1])
            elif fields[:2] == ["Molecular", "Dipole:"]:
                properties["scf_dipole_moment"] = [float(x) for x in fields[2:5]]
            elif fields[:3] == ["SCF", "converged", "in"]:
                properties["scf_iterations"] = int(fields[3])
            # elif fields == ["Gradient units are Hartree/Bohr"]:
            #     # Gradient is stored as (dE/dx1,dE/dy1,dE/dz1,dE/dx2,dE/dy2,...)
            #     for i in range(idx + 3, idx + 3 + natom):
            #         grad = output_lines[i].strip('\n').split()
            #         for x in grad:
            #             gradients.append(float(x))

        if len(gradients) > 0:
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
