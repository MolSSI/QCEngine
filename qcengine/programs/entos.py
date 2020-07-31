"""
The entos QCEngine Harness
"""

import json
import string
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Set, Tuple

import numpy as np
from qcelemental.models import AtomicResult, BasisSet
from qcelemental.util import parse_version, safe_version, which

from ..exceptions import InputError, UnknownError
from ..util import execute, popen
from .model import ProgramHarness
from .util import (
    cca_ao_order_spherical,
    get_ao_conversion,
    reorder_column_ao_indices,
    reorder_row_and_column_ao_indices,
)

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


def entos_ao_order_spherical(max_angular_momentum: int) -> Dict[int, List[int]]:
    ao_order = {}
    for ang_mom in range(max_angular_momentum):
        ao_order[ang_mom] = [x for x in range(ang_mom, -1 * ang_mom - 1, -1)]
    return ao_order


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

    # This map order converts entos ordering to CCA ordering
    # Entos spherical basis ordering for each angular momentum. Follows reverse order of CCA.
    _entos_to_cca_ao_order = {"spherical": get_ao_conversion(cca_ao_order_spherical(10), entos_ao_order_spherical(10))}

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
                "xtb": {"charge": input_model.molecule.molecular_charge},
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
                    }
                }
            # Create the input dictionary for a gradient call
            elif input_model.driver == "gradient":
                input_dict = {
                    "gradient": {
                        **structure,
                        energy_command: {**energy_keywords_from_input_model[energy_command], **energy_keywords_extra},
                    }
                }
            elif input_model.driver == "hessian":
                input_dict = {
                    "hessian": {
                        **structure,
                        energy_command: {**energy_keywords_from_input_model[energy_command], **energy_keywords_extra},
                    }
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
        dft_map = scf_map.copy()
        hf_map = scf_map.copy()
        xtb_map = scf_map.copy()

        energy_command_map = {"dft": dft_map, "hf": hf_map, "xtb": xtb_map}
        extras_map = {"converged": "scf_converged"}
        wavefunction_map = {
            "restricted": {
                "orbitals": "scf_orbitals_a",
                "density": "scf_density_a",
                "fock": "scf_fock_a",
                "eigenvalues": "scf_eigenvalues_a",
                "occupations": "scf_occupations_a",
            },
            "unrestricted": {
                "orbitals_alpha": "scf_orbitals_a",
                "orbitals_beta": "scf_orbitals_b",
                "density_alpha": "scf_density_a",
                "density_beta": "scf_density_b",
                "fock_alpha": "scf_fock_a",
                "fock_beta": "scf_fock_b",
                "eigenvalues_alpha": "scf_eigenvalues_a",
                "eigenvalues_beta": "scf_eigenvalues_b",
                "occupations_alpha": "scf_occupations_a",
                "occupations_beta": "scf_occupations_b",
            },
        }

        # Determine the energy_command
        energy_command = self.determine_energy_command(input_model.model.method)

        gradient_map = {"gradient": "gradient"}
        gradient_map.update({"energy": "scf_total_energy"})
        # TODO Uncomment once entos adds scf_map to gradient json results
        # gradient_map.update(energy_command_map[energy_command])

        hessian_map = {"hessian": "hessian"}
        hessian_map.update(energy_command_map[energy_command])

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
        properties = {}
        load_results = json.loads(outfiles["results.json"])
        entos_results = load_results["json_results"]
        for key in entos_map.keys():
            if key in entos_results:
                properties[entos_map[key]] = entos_results[key]

        # Parse calcinfo_* properties from the results.json
        if "ao_basis" in entos_results.keys():
            properties["calcinfo_nbasis"] = entos_results["ao_basis"]["__Basis"]["n_functions"]
        if "structure" in entos_results.keys():
            properties["calcinfo_natom"] = len(entos_results["structure"]["__Atoms"]["atoms"])

        # Parse wavefunction quantities from entos_results
        wavefunction = {}
        if input_model.protocols.wavefunction != "none":

            # First parse basis set information
            if "ao_basis" in entos_results.keys():
                atom_map = [item[0] for item in entos_results["structure"]["__Atoms"]["atoms"]]

                # Each item in electron_shells is a dictionary containing info for one basis function
                electron_shells_by_center = {}
                for basis_item in entos_results["ao_basis"]["__Basis"]["electron_shells"]:
                    center_index = basis_item["center_index"]

                    electron_shell_info = {
                        "angular_momentum": [basis_item["angular_momentum"]],
                        "harmonic_type": basis_item["function_type"].split("_")[-1],
                        "exponents": basis_item["exponents"],
                        "coefficients": basis_item["coefficients"],
                    }
                    if center_index not in electron_shells_by_center:
                        electron_shells_by_center[center_index] = [electron_shell_info]
                    else:
                        electron_shells_by_center[center_index].append(electron_shell_info)

                # Construct center_data from electron_shells_by_center
                # Note: Duplicate atoms will over write each other
                center_data = {}
                for i in range(len(electron_shells_by_center)):
                    basis_center_info = {"electron_shells": electron_shells_by_center[i]}
                    center_data[atom_map[i]] = basis_center_info

                # Construct BasisSet
                basis_info = {
                    "name": input_model.model.basis,
                    # "description": "", # None provided by entos
                    "center_data": center_data,
                    "atom_map": atom_map,
                    "nbf": entos_results["ao_basis"]["__Basis"]["n_functions"],
                }
                basis_set = BasisSet(**basis_info)
                wavefunction["basis"] = basis_set
            else:
                raise KeyError(
                    f"Basis set information not found so wavefunction protocol {input_model.protocols.wavefunction} is not available."
                )

            # Now parse wavefunction information
            n_channels = entos_results["n_channels"]
            if n_channels == 1:
                wavefunction["restricted"] = True
                for key in wavefunction_map["restricted"].keys():
                    if key in entos_results:
                        if "orbitals" in key:
                            orbitals_transposed = reorder_column_ao_indices(
                                np.array(entos_results[key]), basis_set, self._entos_to_cca_ao_order
                            )
                            wavefunction[wavefunction_map["restricted"][key]] = orbitals_transposed.transpose()
                        elif "density" in key or "fock" in key:
                            wavefunction[wavefunction_map["restricted"][key]] = reorder_row_and_column_ao_indices(
                                entos_results[key], basis_set, self._entos_to_cca_ao_order
                            )
                        else:
                            wavefunction[wavefunction_map["restricted"][key]] = entos_results[key]
            # TODO Add a test in QCEngineRecords
            elif n_channels == 2:
                wavefunction["restricted"] = False
                for key in wavefunction_map["unrestricted"].keys():
                    if key in entos_results:
                        if "orbitals" in key:
                            orbitals_transposed = reorder_column_ao_indices(
                                np.array(entos_results[key]), basis_set, self._entos_to_cca_ao_order
                            )
                            wavefunction[wavefunction_map["restricted"][key]] = orbitals_transposed.transpose()
                        elif "density" in key or "fock" in key:
                            wavefunction[wavefunction_map["restricted"][key]] = reorder_row_and_column_ao_indices(
                                entos_results[key], basis_set, self._entos_to_cca_ao_order
                            )
                        else:
                            wavefunction[wavefunction_map["restricted"][key]] = entos_results[key]

        # Parse results for the extras_map from results.json
        extras = {}
        for key in extras_map.keys():
            if key in entos_results:
                extras[extras_map[key]] = entos_results[key]

        # Initialize output_data by copying over input_model.dict()
        output_data = input_model.dict()

        # Determine the correct return_result
        if input_model.driver == "energy":
            if "scf_total_energy" in properties:
                output_data["return_result"] = properties["scf_total_energy"]
            else:
                raise KeyError(f"Could not find {input_model.model} total energy")
        elif input_model.driver == "gradient" or input_model.driver == "hessian":
            if input_model.driver in properties:
                output_data["return_result"] = properties.pop(input_model.driver)
            else:
                raise KeyError(f"{input_model.driver} not found.")
        else:
            raise NotImplementedError(f"Driver {input_model.driver} not implemented for entos.")

        output_data["properties"] = properties
        if input_model.protocols.wavefunction != "none":
            output_data["wavefunction"] = wavefunction
        output_data["extras"].update(extras)
        output_data["schema_name"] = "qcschema_output"
        output_data["success"] = True

        return AtomicResult(**output_data)

    # Determine the energy_command
    def determine_energy_command(self, method: str) -> str:
        """
        Determine the energy command in entos
        """

        if method.upper() in self._dft_functionals:
            energy_command = "dft"
        else:
            energy_command = method.lower()

        # Check method is supported
        energy_commands = [key for key in self._energy_commands]
        if energy_command not in energy_commands:
            raise InputError(f"The energy method {method} is not implemented in QCEngine for entos.")

        return energy_command
