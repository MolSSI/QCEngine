"""
Calls the Orca executable.
"""

import string
import io
import numpy as np
from typing import Any, Dict, List, Optional, Set, Tuple

import qcelemental
from qcelemental.util import parse_version, which
from qcelemental.models import AtomicResult

from ..exceptions import InputError, UnknownError
from ..util import execute
from .model import ProgramHarness


class OrcaHarness(ProgramHarness):
    _defaults: Dict[str, Any] = {
        "name": "Orca",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    # Set of implemented dft functionals in Orca according to the manual.
    # fmt: off
    _dft_functionals: Set[str] = {
        "HFS" ,"LDA" , "LSD" , "VWN", "VWN5" , "VWN3" , "PWLDA" , "BP86" ,
        "BLYP" , "OLYP" , "GLYP" , "XLYP" , "PW91" , "mPWPW" , "mPWLYP" ,
        "PBE" , "RPBE" , "REVPBE" , "PWP" , "B1LYP" , "B3LYP" , "B3LYP/G" ,
        "O3LYP" , "X3LYP" , "B1P" , "B3P" , "B3PW" , "PW1PW" , "PBE0" ,
        "PW6B95" , "BHANDHLYP" , "TPSS" , "TPSS0" , "M06L" , "M062X" , "wB97"
        , "wB97X" , "wB97X-D3" , "CAM-B3LYP" , "LC-BLYP" , "B2PLYP" ,
        "RI-B2PLYP" , "B2PLYP" , "B2PLYP-D" , "B2PLYP-D3" , "RI-B2PLYP" ,
        "mPW2PLYP-D" , "B2GP-PLYP" , "B2K-PLYP" , "B2T-PLYP" , "PWPB95" ,
        "RI-PWPB95" , "D3BJ" , "D3ZERO" , "D2"
        }

    # fmt: on

    # NOTE: Unrestricted SCF methods must be specified by using keyword reference
    _hf_methods: Set[str] = {"HF", "RHF"}
    _restricted_post_hf_methods: Set[str] = {"MP2", "CCSD", "CCSD(T)"}  # RMP2, RCCSD, RCCSD(T)}
    # TODO Add keyword to specify unrestricted for WF method
    # _unrestricted_post_hf_methods: Set[str] = {"UMP2", "UCCSD", "UCCSD(T)"}
    _post_hf_methods: Set[str] = {*_restricted_post_hf_methods}

    class Config(ProgramHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        return which(
            "orca",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via https://orcaforum.kofo.mpg.de/app.php/portal.",
        )

    # TODO Consider changing this to use molpro --version instead of performing a full execute
    def get_version(self) -> str:
        self.found(raise_error=True)
        which_prog = which("orca")
        if which_prog not in self.version_cache:
            success, output = execute(
                [which_prog, "version.inp", "-d", ".", "-W", "."], infiles={"version.inp": ""}, outfiles=["version.out"]
            )

        _output = output["stdout"].split()

        for index, line in enumerate(_output):
            if "Version" in line:
                found = True
                if found:
                    version = _output[index + 1]
                    break

        return version

    def compute(self, input_data: "AtomicInput", config: "JobConfig") -> "AtomicResult":
        """
        Run Orca
        """
        # Check if Orca executable is found
        self.found(raise_error=True)

        # Check Orca version
        if parse_version(self.get_version()) < parse_version("4.2.1"):
            raise TypeError("Orca version '{}' not supported".format(self.get_version()))

        # Setup the job
        job_inputs = self.build_input(input_data, config)

        # Run Orca
        exe_success, proc = self.execute(job_inputs)

        # Determine whether the calculation succeeded
        if exe_success:
            # If execution succeeded, collect results
            result = self.parse_output(proc, input_data)
            return result
        else:
            # Return UnknownError for error propagation
            return UnknownError(proc["stderr"])

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_infiles: Optional[Dict[str, str]] = None,
        extra_outfiles: Optional[List[str]] = None,
        as_binary: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        scratch_messy: bool = False,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict[str, Any]]:
        """
        For option documentation go look at qcengine/util.execute
        """
        infiles = inputs["infiles"]
        if extra_infiles is not None:
            infiles.update(extra_infiles)

        # Collect all output files and update with extra_outfiles
        outfiles = ["dispatch.engrad"]

        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        # Replace commands with extra_commands if present
        commands = inputs["commands"]
        if extra_commands is not None:
            commands = extra_commands

        # Run the Orca program
        exe_success, proc = execute(
            commands,
            infiles=infiles,
            outfiles=outfiles,
            as_binary=as_binary,
            scratch_name=scratch_name,
            scratch_directory=inputs["scratch_directory"],
            scratch_messy=scratch_messy,
            timeout=timeout,
        )
        return exe_success, proc

    def build_input(
        self, input_model: "AtomicInput", config: "JobConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:

        if template is None:
            input_file = []

            # Resolving keywords
            # caseless_keywords = {k.lower(): v for k, v in input_model.keywords.items()}

            # Write the geom
            xyz_block = input_model.molecule.to_string(dtype="orca", units="Bohr")

            # Determine what SCF type (restricted vs. unrestricted)
            hf_type = "RHF"
            # dft_type = "RKS"
            # if unrestricted:
            #     hf_type = "UHF"
            #     dft_type = "UKS"

            # Write energy call
            energy_call = []

            if input_model.model.method.upper() in self._post_hf_methods:  # post SCF case
                energy_call.append("")
                input_file.append("! {} {}".format(input_model.model.method, input_model.model.basis))
                input_file.append(xyz_block)

            elif input_model.model.method.upper() in self._dft_functionals or self._hf_methods:  # DFT case
                input_file.append("! SP {} {}".format(input_model.model.method, input_model.model.basis))
                input_file.append(xyz_block)

            else:
                raise InputError(f"Method {input_model.model.method} not implemented for Orca.")

            # Write appropriate driver call
            if input_model.driver == "energy":
                input_file.extend(energy_call)
            elif input_model.driver == "gradient":
                input_file[0] = "{} {}".format(input_file[0], "engrad")
            else:
                raise InputError(f"Driver {input_model.driver} not implemented for Orca.")

            parallel = "%pal nproc {} end".format(config.ncores)
            input_file.append(parallel)

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
            "commands": [which("orca"), "dispatch.mol"],
            "infiles": {"dispatch.mol": input_file},
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
        }

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":
        try:
            import cclib
        except ImportError:
            raise ("You need to install the cclib python module...")

        data = cclib.io.ccread(io.StringIO(outfiles["stdout"]))

        properties = {}
        extras = {}
        extras["ev_to_hartrees"] = qcelemental.constants.conversion_factor("eV", "hartree")

        # Process basis set data
        properties["calcinfo_nbasis"] = data.nbasis
        properties["calcinfo_nmo"] = data.nmo
        properties["calcinfo_natom"] = data.natom
        properties["scf_dipole_moment"] = data.moments[1].tolist()

        # Grab the method from input
        method = input_model.model.method.upper()

        # Determining the final energy
        # Throws an error if the energy isn't found for the method specified from the input_model.
        try:
            final_energy = data.scfenergies[-1] * extras["ev_to_hartrees"]
        except:
            raise KeyError(f"Could not find {method} total energy")

        # Initialize output_data by copying over input_model.dict()
        output_data = input_model.dict()

        # Determining return_result
        if input_model.driver == "energy":
            output_data["return_result"] = final_energy
            extras["CURRENT ENERGY"] = final_energy

        elif input_model.driver == "gradient":
            gradient = self._get_gradient(outfiles["outfiles"]["dispatch.engrad"])
            output_data["return_result"] = gradient
            extras["CURRENT ENERGY"] = final_energy
            extras["CURRENT GRADIENT"] = gradient

        # Final output_data assignments needed for the AtomicResult object

        output_data["properties"] = properties
        output_data["extras"].update(extras)
        output_data["schema_name"] = "qcschema_output"
        output_data["stdout"] = outfiles["stdout"]
        output_data["success"] = True

        return AtomicResult(**output_data)

    def _get_gradient(self, gradient_file):
        """Get gradient from engrad Orca file"""
        copy = False
        found = False
        gradient = []

        for line in gradient_file.splitlines():
            if "gradient" in line:
                found = True
                copy = True
            if found and copy:
                try:
                    gradient.append(float(line))
                except ValueError:
                    pass

        dim = int(len(gradient) / 3)

        return np.array(gradient).reshape(dim, 3)
