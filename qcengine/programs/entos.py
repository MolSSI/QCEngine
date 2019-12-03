"""
The entos QCEngine Harness
"""

import json
import string
from typing import Any, Dict, List, Optional, Set, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from ..config import TaskConfig
    from qcelemental.models import AtomicInput

from qcelemental.models import AtomicResult
from qcelemental.util import parse_version, safe_version, which

from ..exceptions import InputError, UnknownError
from ..util import execute, popen
from .model import ProgramHarness


class EntosHarness(ProgramHarness):
    _defaults: Dict[str, Any] = {
        "name": "entos",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    # List of DFT functionals
    _dft_functionals: Set[str] = {
        "SLATER",
        "DIRAC",
        "SLATERD3",
        "DIRACD3",
        "VWN5",
        "VWN",
        "VWN1",
        "SVWN",
        "LDA",
        "BLYP",
        "BPW91",
        "BLYPD3",
        "B88",
        "PBEX",
        "PBERX",
        "PBEC",
        "LYP",
        "PW91",
        "P86",
        "PBE",
        "PBER",
        "PBED3",
        "PBERD3",
        "B3LYP3",
        "B3LYP",
        "B3LYP5",
        "PBE0",
        "PBE1PBE",
        "B3LYP3D3",
        "B3LYPD3",
        "B3LYP5D3",
        "PBE0D3",
        "PBE1PBED3",
        "CAMB3LYP",
        "WB97X",
        "CAMB3LYPD3",
        "WB97XD3",
    }

    # Available keywords for each of the energy commands
    _scf_keywords_extra: Set[str] = {
        "ansatz",
        "df",
        "orbital_grad_threshold",
        "energy_threshold",
        "coulomb_method",
        "exchange_method",
    }
    _dft_keywords_extra: Set[str] = _scf_keywords_extra.copy()
    _hf_keywords_extra: Set[str] = _scf_keywords_extra.copy()
    _xtb_keywords_extra: Set[str] = {}

    # Energy commands that are currently supported and their available keywords
    _energy_commands: Dict[str, Any] = {
        "dft": _dft_keywords_extra,
        "hf": _hf_keywords_extra,
        "xtb": _xtb_keywords_extra,
    }

    class Config(ProgramHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which(
            "entos", return_bool=True, raise_error=raise_error, raise_msg="Please install via https://www.entos.info/"
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("entos")
        if which_prog not in self.version_cache:
            with popen([which_prog, "--version"]) as exc:
                exc["proc"].wait(timeout=15)
            self.version_cache[which_prog] = safe_version(exc["stdout"].split()[2])

        return self.version_cache[which_prog]

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Run entos
        """
        # Check if entos executable is found
        self.found(raise_error=True)

        # Check entos version
        if parse_version(self.get_version()) < parse_version("0.7.1"):
            raise TypeError(f"entos version {self.get_version()} not supported")

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

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_infiles: Optional[Dict[str, str]] = None,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        scratch_messy: bool = False,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:
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
        exe_success, proc = execute(
            commands,
            infiles=infiles,
            outfiles=outfiles,
            scratch_name=scratch_name,
            scratch_directory=inputs["scratch_directory"],
            scratch_messy=scratch_messy,
            timeout=timeout,
        )

        # Entos does not create an output file and only prints to stdout
        proc["outfiles"]["results.json"] = proc["stdout"]
        return exe_success, proc

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:

        # Write the geom xyz file with unit au
        # TODO Make entos dtype in QCElemental
        xyz_file = input_model.molecule.to_string(dtype="xyz", units="Angstrom")

        # Create input file
        if template is None:

            # Determine the energy_command
            energy_command = self.determine_energy_command(input_model.model.method)

            # Define base options for the energy command (options that can be taken directly from input_model)
            scf_keywords_from_input_model = {
                "ao": input_model.model.basis,
                "charge": input_model.molecule.molecular_charge,
                "spin": float(input_model.molecule.molecular_multiplicity - 1),
            }
            energy_keywords_from_input_model = {
                "dft": {"xc": input_model.model.method.upper(), **scf_keywords_from_input_model},
                "hf": {**scf_keywords_from_input_model},
                "xtb": {"charge": input_model.molecule.molecular_charge,},
            }

            # Resolve keywords (extra options) for the energy command
            caseless_keywords = {k.lower(): v for k, v in input_model.keywords.items()}
            energy_keywords_extra = {}
            for key in caseless_keywords.keys():
                if key in self._energy_commands[energy_command]:
                    energy_keywords_extra[key] = caseless_keywords[key]

            # Additional sub trees
            structure = {"structure": {"file": "geometry.xyz"}}  # Structure sub tree

            # Create the input dictionary for a energy call
            if input_model.driver == "energy":
                input_dict = {
                    energy_command: {
                        **structure,
                        **energy_keywords_from_input_model[energy_command],
                        **energy_keywords_extra,
                    },
                }
            # Create the input dictionary for a gradient call
            elif input_model.driver == "gradient":
                input_dict = {
                    "gradient": {
                        **structure,
                        energy_command: {**energy_keywords_from_input_model[energy_command], **energy_keywords_extra,},
                    },
                }
            elif input_model.driver == "hessian":
                input_dict = {
                    "hessian": {
                        **structure,
                        energy_command: {**energy_keywords_from_input_model[energy_command], **energy_keywords_extra},
                    },
                }
            else:
                raise NotImplementedError(f"Driver {input_model.driver} not implemented for entos.")

            # Write input file
            input_file = [
                f"entos_policy(version = '{self.get_version()}')",
                "json_results := ",
            ] + self.write_input_recursive(input_dict)
            input_file = "\n".join(input_file)

        # Use the template input file if present
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
            "infiles": {"dispatch.in": input_file, "geometry.xyz": xyz_file},
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
        }

    def write_input_recursive(self, d: Dict[str, Any]) -> List:
        input_file = []
        for key, value in d.items():
            if isinstance(value, dict):
                input_file.append(key + "(")
                rec_input = self.write_input_recursive(value)
                indented_line = map(lambda x: "  " + x, rec_input)
                input_file.extend(indented_line)
                input_file.append(")")
            else:
                if isinstance(value, str):
                    input_file.append("{0} = '{1}'".format(key, value))
                elif isinstance(value, bool):
                    input_file.append("{0} = {1}".format(key, str(value).lower()))
                else:
                    input_file.append("{0} = {1}".format(key, value))
        return input_file

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":

        scf_map = {"energy": "scf_total_energy", "n_iter": "scf_iterations"}
        scf_extras = {
            "converged": "scf_converged",
            "ao_basis": {"basis": {"n_functions", "shells"}},
            "density": "scf_density",
            "orbitals": "scf_orbitals",
            "fock": "fock",
        }
        dft_map = scf_map.copy()

        hf_map = scf_map.copy()

        xtb_map = scf_map.copy()

        energy_command_map = {"dft": dft_map, "hf": hf_map, "xtb": xtb_map}

        # Determine the energy_command
        energy_command = self.determine_energy_command(input_model.model.method)

        gradient_map = {"gradient": "gradient"}
        gradient_map.update({"energy": "scf_total_energy"})
        # TODO Uncomment once entos adds scf_map to gradient json results
        # gradient_map.update(energy_command_map[energy_command])

        hessian_map = {"hessian": "hessian"}
        hessian_map.update(energy_command_map[energy_command])

        # Initialize properties dictionary
        properties = {}

        # Determine whether to use the energy map or the gradient map
        if input_model.driver == "energy":
            entos_map = energy_command_map[energy_command]
        elif input_model.driver == "gradient":
            entos_map = gradient_map
        elif input_model.driver == "hessian":
            entos_map = hessian_map
        else:
            raise NotImplementedError(f"Driver {input_model.driver} not implemented for entos.")

        # Parse the results.json output from entos
        load_results = json.loads(outfiles["results.json"])
        entos_results = load_results["json_results"]
        for key in entos_map.keys():
            if key in entos_results:
                properties[entos_map[key]] = entos_results[key]

        # Determine the correct return_result
        output_data = {}
        if input_model.driver == "energy":
            if "scf_total_energy" in properties:
                output_data["return_result"] = properties["scf_total_energy"]
            else:
                raise KeyError(f"Could not find {input_model.model} total energy")
        elif input_model.driver == "gradient":
            if "gradient" in properties:
                output_data["return_result"] = properties.pop("gradient")
            else:
                raise KeyError("Gradient not found.")
        elif input_model.driver == "hessian":
            if "hessian" in properties:
                output_data["return_result"] = properties.pop("hessian")
            else:
                raise KeyError("Hessian not found.")
        else:
            raise NotImplementedError(f"Driver {input_model.driver} not implemented for entos.")

        output_data["properties"] = properties
        output_data["schema_name"] = "qcschema_output"
        output_data["success"] = True

        return AtomicResult(**{**input_model.dict(), **output_data})

    # Determine the energy_command
    def determine_energy_command(self, method):
        """
        Determine the energy command in entos
        """

        if method.upper() in self._dft_functionals:
            energy_command = "dft"
        elif method.upper() == "HF":
            energy_command = "hf"
        elif method.upper() == "XTB":
            energy_command = "xtb"
        else:
            energy_command = method.lower()

        # Check method is supported
        energy_commands = [key for key in self._energy_commands]
        if energy_command not in energy_commands:
            raise InputError(f"Energy method, {method}, not implemented for entos.")

        return energy_command
