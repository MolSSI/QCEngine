"""
Calls the Orca executable.
"""

import string
import cclib
import io
import numpy as np
from typing import Any, Dict, List, Optional, Set, Tuple

from qcelemental.models import AtomicResult
from qcelemental.util import parse_version, which

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

    # Set of implemented dft functionals in Molpro according to dfunc.registry (version 2019.2)
    # fmt: off
    _dft_functionals: Set[str] = {
        "B86MGC", "B86R", "B86", "B88C", "B88", "B95", "B97DF", "B97RDF", "BR", "BRUEG", "BW", "CS1", "CS2", "DIRAC",
        "ECERFPBE", "ECERF", "EXACT", "EXERFPBE", "EXERF", "G96", "HCTH120", "HCTH147", "HCTH93", "HJSWPBEX", "LTA",
        "LYP", "M052XC", "M052XX", "M05C", "M05X", "M062XC", "M062XX", "M06C", "M06HFC", "M06HFX", "M06LC", "M06LX",
        "M06X", "M12C", "MK00B", "MK00", "P86", "PBEC", "PBESOLC", "PBESOLX", "PBEXREV", "PBEX", "PW86", "PW91C",
        "PW91X", "PW92C", "STEST", "TFKE", "TH1", "TH2", "TH3", "TH4", "THGFCFO", "THGFCO", "THGFC", "THGFL", "TPSSC",
        "TPSSX", "VSXC", "VW", "VWN3", "VWN5", "XC-M05-2X", "XC-M05", "XC-M06-2X", "XC-M06-HF", "XC-M06-L", "XC-M06",
        "XC-M08-HX", "XC-M08-SO", "XC-M11-L", "XC-SOGGA11", "XC-SOGGA11-X", "XC-SOGGA", "FRMTST", "LHF", "TLHF",
        "LXBECKE", "ELP", "NULL", "YTEST", "TREF2", "TREF", "TTEST", "GLE", "GREEN", "SRB88", "SRLYP", "LB94", "EI",
        "SAOP", "USER", "INE", "ECERF2", "ECERFINTER", "ECERFLOCAL2", "ECERFLOCAL", "EXERFLOCAL", "FC", "FCFO", "FCO",
        "FL", "XC-M11", "PBEXANAL", "PBECANAL", "PBESOLCANAL", "PBESOLXANAL", "EXSRLDA", "EXSRLPBE", "ECSRLPBE",
        "ECSRLLPBE", "ECSQRTLPBE", "ECMUDIVLPBE", "EXERFPHS", "ECLERFMUPBE", "ECERFERFCPBE", "ECSQRTLDA", "REVPBEX",
        "B", "B-LYP", "BLYP", "B-P", "BP86", "B-VWN", "B3LYP", "B3LYP3", "B3LYP5", "B88X", "B97", "B97R", "BECKE",
        "BH-LYP", "CS", "D", "HFB", "HFS", "LDA", "LSDAC", "LSDC", "KYP88", "MM05", "MM05-2X", "MM06", "MM06-2X",
        "MM06-L", "MM06-HF", "PBE", "PBE0", "PBE0MOL", "PBEREV", "PW91", "S", "S-VWN", "SLATER", "VS99", "VWN",
        "VWN80", 'M05', "M05-2X", "M06", "M06-2X", "M06-L", "M06-HF", "M08-HX", "M08-SO", "M11-L", "TPSS", "TPSSH",
        "M12HFC", "HJSWPBE", "HJSWPBEH", "TCSWPBE", "PBESOL"
    }
    # fmt: on

    # Different open-shell scenarios:
    # - Restricted references
    #       RHF-RMP2 (like RHF-UMP2 using Molpro's naming convention for CC methods)
    #       RHF-UCCSD
    #       RHF-UCCSD(T)
    #       RHF-RCCSD
    #       RHF-RCCSD(T)
    # - Unrestricted references (Only supported up to UMP2, no CC support)
    #       UHF-UMP2
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

        # name_space = {"molpro_uri": "http://www.molpro.net/schema/molpro-output"}
        # which_prog = which("molpro")
        # if which_prog not in self.version_cache:
        #     success, output = execute(
        #         [which_prog, "version.inp", "-d", ".", "-W", "."],
        #         infiles={"version.inp": ""},
        #         outfiles=["version.out", "version.xml"],
        #     )

        #     if success:
        #         tree = ET.ElementTree(ET.fromstring(output["outfiles"]["version.xml"]))
        #         root = tree.getroot()
        #         version_tree = root.find("molpro_uri:job/molpro_uri:platform/molpro_uri:version", name_space)
        #         year = version_tree.attrib["major"]
        #         minor = version_tree.attrib["minor"]
        #         molpro_version = year + "." + minor
        #         self.version_cache[which_prog] = safe_version(molpro_version)

        return "4.2.1"

    def compute(self, input_data: "AtomicInput", config: "JobConfig") -> "AtomicResult":
        """
        Run Orca
        """
        # Check if Molpro executable is found
        self.found(raise_error=True)

        # Check Orca version
        if parse_version(self.get_version()) < parse_version("4.2.1"):
            raise TypeError("Orca version '{}' not supported".format(self.get_version()))

        # Setup the job
        job_inputs = self.build_input(input_data, config)

        # Run Orca
        binary = ["dispatch.gbw"]
        exe_success, proc = self.execute(job_inputs, as_binary=binary)

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
        outfiles = ["dispatch.gbw", "dispatch.engrad"]

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
            caseless_keywords = {k.lower(): v for k, v in input_model.keywords.items()}
            unrestricted = False
            # Following Molpro default that ROHF is done for open-shell calculations unless unrestricted is specified
            if "reference" in caseless_keywords and caseless_keywords["reference"] == "unrestricted":
                unrestricted = True

            # Memory is in megawords per core for Molpro
            # memory_mw_core = int(config.memory * (1024 ** 3) / 8e6 / config.ncores)
            # input_file.append("memory,{},M".format(memory_mw_core))
            # input_file.append("")

            # Write the geom
            xyz_block = input_model.molecule.to_string(dtype="orca", units="Angstrom")
            # input_file.append(xyz_block)

            # Write the basis set
            # input_file.append("basis={")
            # input_file.append(f"default,{input_model.model.basis}")
            # input_file.append("}")
            # input_file.append("")

            # Determine what SCF type (restricted vs. unrestricted)
            hf_type = "RHF"
            dft_type = "RKS"
            if unrestricted:
                hf_type = "UHF"
                dft_type = "UKS"

            # Write energy call
            energy_call = []
            # If post-hf method is called then make sure to write a HF call first
            if input_model.model.method.upper() in self._post_hf_methods:  # post SCF case
                energy_call.append(f"{{{hf_type}}}")
                energy_call.append("")
                energy_call.append(f"{{{input_model.model.method}}}")
            # If DFT call make sure to write {rks,method}
            elif input_model.model.method.upper() in self._dft_functionals:  # DFT case
                input_file.append("! SP {} {}".format(input_model.model.method, input_model.model.basis))
                input_file.append(xyz_block)

                # energy_call.append(f"{{{dft_type},{input_model.model.method}}}")
            elif input_model.model.method.upper() in self._hf_methods:  # HF case
                energy_call.append(f"{{{hf_type}}}")
            else:
                raise InputError(f"Method {input_model.model.method} not implemented for Molpro.")

            # Write appropriate driver call
            if input_model.driver == "energy":
                input_file.extend(energy_call)
            elif input_model.driver == "gradient":
                input_file[0] = "{} {}".format(input_file[0], "engrad")
            else:
                raise InputError(f"Driver {input_model.driver} not implemented for Molpro.")

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

        data = cclib.io.ccread(io.StringIO(outfiles["stdout"]))

        properties = {}
        extras = {}

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
            final_energy = data.scfenergies[-1] * 0.0367493
        except:
            raise KeyError(f"Could not find {method} total energy")

        # Initialize output_data by copying over input_model.dict()
        output_data = input_model.dict()

        # Determining return_result
        if input_model.driver == "energy":
            output_data["return_result"] = final_energy 
            extras["CURRENT ENERGY"] = final_energy
        elif input_model.driver == "gradient":
            gradient = self.get_gradient(outfiles["outfiles"]["dispatch.engrad"])
            output_data["return_result"] = gradient
            extras["CURRENT GRADIENT"] = gradient

        # Final output_data assignments needed for the AtomicResult object

        output_data["properties"] = properties
        output_data["extras"].update(extras)
        output_data["schema_name"] = "qcschema_output"
        output_data["stdout"] = outfiles["stdout"]
        output_data["success"] = True

        import uuid
        myid = str(uuid.uuid4())
        with open("/tmp/orcafiles/{}".format(myid), "wb") as handle:
            handle.write(outfiles["outfiles"]["dispatch.gbw"])
        output_data["extras"]["gbw"] = myid

        return AtomicResult(**output_data)

    def get_gradient(self, gradient_file):
        """Get gradient from engrad Orca file
        """
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
        
        dim = np.sqrt(len(gradient)).astype(int)

        return np.array(gradient).reshape(dim, dim)