"""
Calls the Q-Chem executable.
"""

import os
import string
import json
from typing import Any, Dict, List, Optional, Set, Tuple

from qcelemental.models import Result
from qcelemental.util import parse_version, safe_version, which

from ..exceptions import InputError, UnknownError
from ..util import execute, popen, temporary_directory
from .model import ProgramHarness


class QChemHarness(ProgramHarness):
    _defaults: Dict[str, Any] = {
        "name": "qchem",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    def found(self, raise_error: bool = False) -> bool:
        return which('qchem',
                     return_bool=True,
                     raise_error=raise_error,
                     raise_msg='Please install by visiting the Q-Chem website.')

    def _get_qc_path(self, config: Optional['JobConfig'] = None):
        paths = os.environ.copy()
        paths["QCSCRATCH"] = "/tmp"
        if config and config.scratch_directory:
            paths["QCSCRATCH"] = config.scratch_directory

        # Nothing else to pass in
        if "QC" in os.environ:
            return paths

        # Assume QC path is
        qchem_path = which('qchem')
        if qchem_path:
            paths["QC"] = os.path.dirname(os.path.dirname(qchem_path))

        return paths

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which('qchem')
        if which_prog not in self.version_cache:
            with popen([which_prog, '-h'], popen_kwargs={"env": self._get_qc_path()}) as exc:
                exc["proc"].wait(timeout=15)

            if "QC not defined" in exc["stdout"]:
                return safe_version("0.0.0")

            self.version_cache[which_prog] = safe_version(exc["stdout"].splitlines()[0].split()[-1])

        return self.version_cache[which_prog]

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Run qchem
        """
        # Check if qchem executable is found
        self.found(raise_error=True)

        # Check qchem version
        if parse_version(self.get_version()) < parse_version("5.2"):
            raise TypeError("Q-Chem version '{}' not supported".format(self.get_version()))

        # Setup the job
        job_inputs = self.build_input(input_data, config)

        # Run qchem
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
        outfiles = ["dispatch.out"]
        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        # Replace commands with extra_commands if present
        commands = inputs["commands"]
        if extra_commands is not None:
            commands = extra_commands

        envs = self._get_qc_path()

        with temporary_directory(parent=inputs["scratch_directory"], suffix="_qchem_scratch") as tmpdir:
            envs["QCSCRATCH"] = os.path.join(tmpdir, "scratch").replace("//", "/")
            exe_success, proc = execute(commands,
                                        infiles=infiles,
                                        outfiles=outfiles,
                                        scratch_name=scratch_name,
                                        scratch_directory=tmpdir,
                                        scratch_messy=scratch_messy,
                                        timeout=timeout,
                                        environment=envs)

        # QChem does not create an output file and only prints to stdout
        return exe_success, proc

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:

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
            keywords["NO_REORIENT"] = "TRUE"

        keywords["METHOD"] = input_model.model.method
        if input_model.model.basis:
            keywords["BASIS"] = input_model.model.basis

        # Begin the input file
        input_file = []
#         input_file.append(f"""$comment
# Automatically generated Q-Chem input file by QCEngine
# $end
#             """)

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
            "infiles": {
                "dispatch.in": "\n".join(input_file)
            },
            "commands": [which("qchem"), "-nt", str(config.ncores), "dispatch.in", "dispatch.out"],
            "scratch_directory": config.scratch_directory,
        }

        return ret

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':

        output_data = {}

        outfile = outfiles["dispatch.out"]
        print(outfile)



        output_data["return_result"] = 5
        output_data["properties"] = {}
        output_data['stdout'] = outfiles["dispatch.out"]
        output_data['success'] = True


        return Result(**{**input_model.dict(), **output_data})


