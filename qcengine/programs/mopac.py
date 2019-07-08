"""
Calls the Psi4 executable.
"""
import json
import os
from typing import Any, Dict, Optional

from qcelemental.models import Result
from qcelemental.util import parse_version, safe_version, which

from .model import ProgramHarness
from ..exceptions import InputError, RandomError, ResourceError, UnknownError
from ..util import execute, popen, temporary_directory


class MopacHarness(ProgramHarness):

    _defaults = {
        "name": "MOPAC",
        "scratch": True, # Input/output file
        "thread_safe": True,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    def __init__(self, **kwargs):
        extras = {  # All units taken from within MOPAC
            "bohr_to_angstroms": 0.5291772083,
            "hartree_to_ev": 27.2113834,
            "ev_to_kcalmol": 23.060529,
        }
        extras["au_to_debye"] = 2.99792458e10 * 1.602176462e0 * 1e-10 * extras["bohr_to_angstroms"]

        # convert from kcal/mol to Hartree
        extras["energy_conversion"] = 1.0 / (extras["hartree_to_ev"] * extras["ev_to_kcalmol"])

        # convert from kcal/mol/angstrom to Hartree/bohr
        extras["force_conversion"] = extras["bohr_to_angstroms"] / (extras["hartree_to_ev"] * extras["ev_to_kcalmol"])

        # convert from Debye to atomic units
        extras["dipole_conversion"] = 1.0 / extras["au_to_debye"]

        kwargs["extras"] = extras
        super().__init__(**{**self._defaults, **kwargs})

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which('mopac',
                     return_bool=True,
                     raise_error=raise_error,
                     raise_msg='Please install via http://openmopac.net.')

    def get_version(self) -> str:
        self.found(raise_error=True)

        # Not really possible to pull at the moment, MolSSI will add a version ability
        return "2016"

    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs Psi4 in API mode
        """
        self.found(raise_error=True)

        exec_commands = self.build_input(input_model, config)

        output = self.execute(exec_commands)
        # print(output)
        print(output["outfiles"]["dispatch.out"])
        # print(output["outfiles"]["dispatch.aux"])
        # print(output["stderr"])

        ret = self.parse_output(output["outfiles"], input_model)


        raise Exception()

    def execute(self, inputs, extra_infiles=None, extra_outfiles=None, extra_commands=None, scratch_name=None,
        scratch_messy=False, timeout=None):
        """
        For option documentation go look at qcengine/util.execute
        """

        infiles = inputs["infiles"]
        if extra_infiles is not None:
            infiles.update(extra_infiles)

        outfiles = inputs["outfiles"]
        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        commands = inputs["commands"]
        if extra_commands is not None:
            commands = extra_commands

        exe_success, proc = execute(commands,
                                    infiles=infiles,
                                    outfiles=outfiles,
                                    scratch_directory=inputs["scratch_directory"],
                                    scratch_name=scratch_name,
                                    scratch_messy=scratch_messy,
                                    timeout=timeout,
                                    environment=inputs.get("environment", None)
                                    )

        # Determine whether the calculation succeeded
        if not exe_success:
            return UnknownError(proc["stderr"])

        return proc

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:

        if template is not None:
            raise KeyError("MOPAC does not currently support input templates.")

        method = input_model.model.method.lower()
        if method not in {"mndo", "am1", "pm3", "rm1", "mndod", "pm6", "pm6-d3", "pm6-dh+", "pm6-dh2", "pm6-dh2x", "pm6-d3h4", "pm6-3dh4x", "pm7", "pm7-ts"}:
            raise InputError(f"MOPAC does not have method: {method.upper()}")

        if input_model.driver not in ["energy", "gradient"]:
            raise InputError(f"MOPAC can only compute energies and gradients, found {input_model.driver}")

        input_file = []

        # 1SCF says not to compute a geometry optimization, always compute a gradient (free), and dump the aux file
        input_file.append(f"{method.upper()} 1SCF GRADIENTS AUX "
                          f"CHARGE={input_model.molecule.molecular_charge} "
                          f"MS={(input_model.molecule.molecular_multiplicity-1)/2}")
        input_file.append("")
        input_file.append("")

        mol = input_model.molecule
        geometry = mol.geometry * self.extras["bohr_to_angstroms"]

        for x in range(len(mol.symbols)):
            if mol.real[x] is False:
                continue

            input_file.append(f"{mol.symbols[x]}  {geometry[x, 0]:17.12f}  {geometry[x, 1]:17.12f}  {geometry[x, 2]:17.12f}")


        env = os.environ.copy()
        env["MKL_NUM_THREADS"] = str(config.ncores)
        env["OMP_NUM_THREADS"] = str(config.ncores)

        return {
            "commands": ["mopac", "dispatch.mop"],
            "infiles": {
                "dispatch.mop": "\n".join(input_file),
            },
            "outfiles": {
                "dispatch.out",
                "dispatch.aux",
            },
            "scratch_directory": config.scratch_directory,
            "environment": env,
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':

        print("======")
        # Setup the job
        data_dict = {}
        in_section = False

        last_key = None
        for line in outfiles["dispatch.aux"].splitlines():
            if ("START" in line) or ("END" in line) or ("#" in line):
                continue

            print(line)
            if "=" in line:
                key, value = line.split("=", 1)
                last_key = key
                print(key, value)
            else:
                pass


            # print(line)


        output_data = input_model.dict()



        return Result(**output_data)
