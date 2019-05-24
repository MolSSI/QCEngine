"""
Calls the entos executable.
"""

from typing import Any, Dict, Optional

from qcelemental.models import Result, FailedOperation
from ..util import execute
from qcelemental.util import which

from .executor import ProgramExecutor


class EntosExecutor(ProgramExecutor):
    _defaults = {
        "name": "entos",
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
        return which('entos', return_bool=True, raise_error=raise_error, raise_msg='Please install via XXX')

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
        Run entos
        """
        # Check if entos executable is found
        self.found(raise_error=True)

        # Check entos version
        # if parse_version(self.get_version()) < parse_version("2018.1"):
        #     raise TypeError("entos version '{}' not understood".format(self.get_version()))

        # Setup the job
        job_inputs = self.build_input(input_data, config)

        # Run entos
        exe_success, proc = execute(job_inputs["commands"],
                                    infiles=job_inputs["infiles"],
                                    outfiles=["dispatch.out"],
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

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:
        input_file = []

        # Write header info
        input_file.append("!Title")
        input_file.append("memory,{},M".format(config.memory))
        input_file.append('')

        # Write the geom
        input_file.append('geometry={')
        for sym, geom in zip(input_model.molecule.symbols, input_model.molecule.geometry):
            s = "{:<4s} {:>{width}.{prec}f} {:>{width}.{prec}f} {:>{width}.{prec}f}".format(
                sym, *geom, width=14, prec=10)
            input_file.append(s)
        input_file.append('}')

        # Write gradient call if asked for
        if input_model.driver == 'gradient':
            input_file.append('')
            input_file.append('{force}')

        input_file = "\n".join(input_file)

        return {
            "commands": ["entos", "dispatch.in", "-n", str(config.ncores)],
            "infiles": {
                "dispatch.in": input_file
            },
            "scratch_location": config.scratch_directory,
            "input_result": input_model.copy(deep=True)
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':

        output_data = {}
        properties = {}

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
        supported_methods = {"HF": scf_maps, "RHF": scf_maps, "MP2": mp2_maps}

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
