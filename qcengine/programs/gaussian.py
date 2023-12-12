"""
Calls the GAUSSIAN executable
"""

import os
import re
import tempfile
import warnings
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple
import cclib
from cclib.method import Nuclear

import numpy as np
from qcelemental import constants
from qcelemental.models import AtomicInput, AtomicResult, Molecule, Provenance
from qcelemental.molparse import regex
from qcelemental.util import parse_version, safe_version, which

from qcengine.config import TaskConfig, get_config

from ..exceptions import InputError, UnknownError
from ..util import disk_files, execute, temporary_directory
from .model import ProgramHarness

class GaussianHarness(ProgramHarness):

    _defaults = {
        "name": "Gaussian",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass


    def found(self, raise_error: bool = False) -> bool:
        return which(
            "g09",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install Gaussian. Check it's in your PATH with `which g09`."
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        # Get the node configuration
        #config = get_config()

        which_prog = which("g09")
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "v.com"], {"v.com": ""})

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: TaskConfig) -> "AtomicResult":
        """
        Run Gaussian
        """
        # Check if Gaussian executable is found
        self.found(raise_error=True)

        # Setup the job
        job_inputs = self.build_input(input_model, config)

        # Run Gaussian
        exe_success, proc = self.execute(job_inputs)

        # Determine whether the calculation succeeded
        if exe_success:
            # If execution succeeded, collect results
            result = self.parse_output(proc, input_model)
            
            return result
 
        else:
            proc['outfiles']['stderr'] = proc['outfiles']['output.log']
            outfile = proc['outfiles']['output.log']
            
            if 'Error termination via ' in outfile:
                raise InputError(proc['outfiles']['output.log'])
            else:
                # Return UnknownError for error propagation
                raise UnknownError(proc["outfiles"]["output.log"])

    def build_input(
        self, input_model: AtomicInput, config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        gaussian_ret = {
            "infiles": {},
            "commands": [which("g09"),  "input.inp", "output.log"],
            "scratch_directory": config.scratch_directory
            }

        input_file = '''%mem=20MW
#P HF/6-31G(d) scf=tight

test1 HF/6-31G(d) sp formaldehyde

0 1
C1
O2  1  r2
H3  1  r3  2  a3
H4  1  r4  2  a4  3  d4

r2=1.20
r3=1.0
r4=1.0
a3=120.
a4=120.
d4=180.
'''
        gaussian_ret['infiles']['input.inp'] = input_file

        return gaussian_ret

    def execute(self,
                inputs,
                extra_outfiles: Optional[Dict[str, str]] = None,
                extra_commands: Optional[List[str]] = None,
                scratch_name = None,
                timeout: Optional[int] = None
               ):

        success, dexe = execute(
            inputs["commands"],
            inputs["infiles"],
            outfiles = ['output.log'],
            scratch_messy = True
            )
        
        if (dexe['outfiles']['output.log'] is None) or (
           'Error termination via' in dexe['outfiles']['output.log']):
           print ('THERE IS AN ERROR!')
           
           success = False
        
        return success, dexe
 
    def parse_output(self, outfiles: Dict[str, str], input_model: AtomicInput) -> AtomicResult:
        output_data = {}
        stdout = outfiles.pop('stdout')
        stderr = outfiles.pop('stderr')

        method = input_model.model.method.lower()
        #method = method[4:] if method.startswith("") else method

        try:
            qcvars, gaussianmol, module = harvest(
                input_model.molecule, method, stdout, **outfiles
            )
        except:
            pass

        #properties = {
        #    "nuclear_repulsion_energy": bdata['99.0'][0],
        #    "scf_total_energy": bdata["99.0"][1],
        #    "return_energy": bdata["99.0"][-1],
        #}

        qcvars = {}

        #props, prov = self._parse_logfile_common(outtext, input_model.dict())
        #output_data["provenance"] = prov
        #output_data["properties"] = properties
        #output_data["properties"].update(props)
        output_data["stdout"] = stdout
        output_data["success"] = True

        merged_data = {**input_model.dict(), **output_data}
        #merged_data["extras"]["qcvars"] = qcvars
        print(merged_data)

        return AtomicResult(**merged_data)
