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
from ..util import execute, popen
from .model import ProgramHarness


class QChemHarness(ProgramHarness):
    _defaults: Dict[str, Any] = {
        "name": "entos",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    def found(self, raise_error: bool = False) -> bool:
        return which('entos', return_bool=True, raise_error=raise_error, raise_msg='Please install by visiting the Q-Chem website.')

    def _get_qc_path(self, config: Optional['JobConfig']=None):
        paths = {"QCSCRATCH": "/tmp"}
        if config and config.scratch_directory:
            paths["QCSCRATCH"] = config.scratch_directory

        # Nothing else to pass in
        if "QC" in os.environ:
            return paths

        # Assume QC path is
        entos_path = which('entos')
        if entos_path:
            paths["QC"] = os.path.dirname(os.path.dirname(entos_path))

        return paths

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which('entos')
        if which_prog not in self.version_cache:
            with popen([which_prog, '-h'], popen_kwargs={"env": self._get_qc_path()}) as exc:
                exc["proc"].wait(timeout=15)
            self.version_cache[which_prog] = safe_version(exc["stdout"].splitlines()[0].split()[-1])

        return self.version_cache[which_prog]

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Run entos
        """
        # Check if entos executable is found
        self.found(raise_error=True)

        # Check entos version
        if parse_version(self.get_version()) < parse_version("5.2"):
            raise TypeError("Q-Chem version '{}' not supported".format(self.get_version()))

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
        outfiles = ["dispatch.out", "results.json"]
        if extra_outfiles is not None:
            outfiles.extend(extra_outfiles)

        # Replace commands with extra_commands if present
        commands = inputs["commands"]
        if extra_commands is not None:
            commands = extra_commands

        # Run the entos program
        exe_success, proc = execute(commands,
                                    infiles=infiles,
                                    outfiles=outfiles,
                                    scratch_name=scratch_name,
                                    scratch_directory=inputs["scratch_directory"],
                                    scratch_messy=scratch_messy,
                                    timeout=timeout)

        # QChem does not create an output file and only prints to stdout
        proc["outfiles"]["results.json"] = proc["stdout"]
        return exe_success, proc

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:

        # Write the geom xyz file with unit au
        xyz_file = input_model.molecule.to_string(dtype='xyz', units='Angstrom')

        # Build keywords
        keywords = {k.upper(): v for k, v in input_model.keywords.items()}
        keywords["INPUT_BOHR"] = "TRUE"

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

        input_file = []
        input_file.append(f"""
$comment
Automatically generated QCEngine input file
$end
            """)
        input_file.append("$molecule")
        input_file.append(input_model.molecule.to_string("psi4"))
        input_file.append("$end")


    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':

        pass