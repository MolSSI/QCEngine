"""
Calls the Q-Chem executable.
"""

import os
import tempfile
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from qcelemental.models import AtomicResult
from qcelemental.util import parse_version, safe_version, which

from ..exceptions import InputError, UnknownError
from ..util import disk_files, execute, popen, temporary_directory
from .model import ProgramHarness


class QChemHarness(ProgramHarness):
    _defaults: Dict[str, Any] = {
        "name": "QChem",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    def found(self, raise_error: bool = False) -> bool:
        return which(
            "qchem",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install by visiting the Q-Chem website.",
        )

    def _get_qc_path(self, config: Optional["JobConfig"] = None):
        paths = os.environ.copy()
        paths["QCSCRATCH"] = tempfile.gettempdir()
        if config and config.scratch_directory:
            paths["QCSCRATCH"] = config.scratch_directory

        # Nothing else to pass in
        if "QC" in os.environ:
            return paths

        # Assume QC path is
        qchem_path = which("qchem")
        if qchem_path:
            paths["QC"] = os.path.dirname(os.path.dirname(qchem_path))

        return paths

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("qchem")
        if which_prog not in self.version_cache:
            with popen([which_prog, "-h"], popen_kwargs={"env": self._get_qc_path()}) as exc:
                exc["proc"].wait(timeout=15)

            if "QC not defined" in exc["stdout"]:
                return safe_version("0.0.0")

            self.version_cache[which_prog] = safe_version(exc["stdout"].splitlines()[0].split()[-1])

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "JobConfig") -> "AtomicResult":
        """
        Run qchem
        """
        # Check if qchem executable is found
        self.found(raise_error=True)

        # Check qchem version
        qceng_ver = "5.2"
        if parse_version(self.get_version()) < parse_version(qceng_ver):
            raise TypeError(f"Q-Chem version <{qceng_ver} not supported (found version {self.get_version()})")

        # Setup the job
        job_inputs = self.build_input(input_model, config)

        # Run qchem
        exe_success, proc = self.execute(job_inputs)

        # Determine whether the calculation succeeded
        if exe_success:
            # If execution succeeded, collect results
            result = self.parse_output(proc["outfiles"], input_model)
            return result
        else:
            outfile = proc["outfiles"]["dispatch.out"]
            if "fatal error occurred in module qparser" in outfile:
                raise InputError(proc["outfiles"]["dispatch.out"])
            else:
                # Return UnknownError for error propagation
                raise UnknownError(proc["outfiles"]["dispatch.out"])

    def execute(
        self,
        inputs: Dict[str, Any],
        *,
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

        binary_files = [os.path.join("savepath", x) for x in ["99.0", "131.0", "132.0"]]

        # Collect all output files and extend with with extra_outfiles
        outfiles = ["dispatch.out"]
        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        # Replace commands with extra_commands if present
        commands = inputs["commands"] + ["savepath"]
        if extra_commands is not None:
            commands = extra_commands

        envs = self._get_qc_path()

        with temporary_directory(parent=inputs["scratch_directory"], suffix="_qchem_scratch") as tmpdir:
            envs["QCSCRATCH"] = tmpdir
            bdict = {x: None for x in binary_files}

            with disk_files({}, bdict, cwd=tmpdir, as_binary=binary_files):
                exe_success, proc = execute(
                    commands,
                    infiles=infiles,
                    outfiles=outfiles,
                    scratch_name=scratch_name,
                    scratch_directory=tmpdir,
                    scratch_messy=scratch_messy,
                    timeout=timeout,
                    environment=envs,
                )

            proc["outfiles"].update({os.path.split(k)[-1]: v for k, v in bdict.items()})

        if (proc["outfiles"]["dispatch.out"] is None) or (
            "Thank you very much for using Q-Chem" not in proc["outfiles"]["dispatch.out"]
        ):
            exe_success = False

        # QChem does not create an output file and only prints to stdout
        return exe_success, proc

    def build_input(
        self, input_model: "AtomicInput", config: "JobConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:

        # Check some bounds on what cannot be parsed
        if "ccsd" in input_model.model.method.lower() or "ccd" in input_model.model.method.lower():
            raise InputError("Cannot handle CC* methods currently.")

        # Build keywords
        keywords = {k.upper(): v for k, v in input_model.keywords.items()}
        keywords["INPUT_BOHR"] = "TRUE"
        keywords["MEM_TOTAL"] = str(int(config.memory * 1024))  # In MB

        if input_model.driver == "energy":
            keywords["JOBTYPE"] = "sp"
        elif input_model.driver == "gradient":
            keywords["JOBTYPE"] = "force"
        elif input_model.driver == "hessian":
            keywords["JOBTYPE"] = "freq"
        else:
            raise InputError(f"Driver of type {input_model.driver} is not yet supported.")

        if input_model.molecule.fix_com or input_model.molecule.fix_orientation:
            keywords["SYM_IGNORE"] = "TRUE"

        keywords["METHOD"] = input_model.model.method
        if input_model.model.basis:
            keywords["BASIS"] = input_model.model.basis

        # Begin the input file
        input_file = []
        input_file.append(
            f"""$comment
Automatically generated Q-Chem input file by QCEngine
$end
            """
        )

        # Add Molecule, TODO: Add to QCElemental
        mol = input_model.molecule
        input_file.append("$molecule")
        input_file.append(f"""{int(mol.molecular_charge)} {mol.molecular_multiplicity}""")

        for real, sym, geom in zip(mol.real, mol.symbols, mol.geometry):
            if real is False:
                raise InputError("Cannot handle ghost atoms yet.")
            input_file.append(f"{sym} {geom[0]:14.8f} {geom[1]:14.8f} {geom[2]:14.8f}")

        input_file.append("$end\n")

        # Write out the keywords
        input_file.append("$rem")
        for k, v in keywords.items():
            input_file.append(f"{k:20s}  {v}")
        input_file.append("$end\n")

        ret = {
            "infiles": {"dispatch.in": "\n".join(input_file)},
            "commands": [which("qchem"), "-nt", str(config.ncores), "dispatch.in", "dispatch.out"],
            "scratch_directory": config.scratch_directory,
        }

        return ret

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":

        output_data = {}
        outfiles["dispatch.out"]

        bdata = {}
        for k, v in outfiles.items():
            if k == "dispatch.out":
                continue
            if v is None:
                continue
            bdata[k] = np.frombuffer(v)

        if input_model.driver == "energy":
            output_data["return_result"] = bdata["99.0"][-1]
        elif input_model.driver == "gradient":
            output_data["return_result"] = bdata["131.0"]
        elif input_model.driver == "hessian":
            output_data["return_result"] = bdata["132.0"]
        else:
            raise ValueError(f"Could not parse drive of type {driver}.")

        properties = {
            "nuclear_repulsion_energy": bdata["99.0"][0],
            "scf_total_energy": bdata["99.0"][1],
            "return_energy": bdata["99.0"][-1],
        }

        # Correct CCSD because its odd?
        # if input_model.model.method.lower() == "ccsd":
        #     m1 = re.findall(" CCSD correlation energy.+=.+\d+\.\d+", outfiles["dispatch.out"])
        #     m2 = re.findall(" CCSD total energy.+=.+\d+\.\d+", outfiles["dispatch.out"])

        output_data["properties"] = properties
        output_data["stdout"] = outfiles["dispatch.out"]
        output_data["success"] = True

        return AtomicResult(**{**input_model.dict(), **output_data})
