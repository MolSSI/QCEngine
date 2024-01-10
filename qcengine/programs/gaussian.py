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

        # Build keywords
        keywords = {k.upper(): v for k, v in input_model.keywords.items()}
        keywords = {'scf_damp': 'true',
                    'scf_diis': 'false'
                    }
        if input_model.driver == "energy":
            keywords["JOBTYPE"] = "sp"
        elif input_model.driver == "gradient":
            keywords["JOBTYPE"] = "force"
        elif input_model.driver == "hessian":
            keywords["JOBTYPE"] = "freq"
        else:
            raise InputError(f"Driver {input_model.driver} not implemented for Gaussian.")

        if input_model.molecule.fix_com or input_model.molecule.fix_orientation:
            keywords["SYM_IGNORE"] = "TRUE"
        
        # Begin input file
        input_file = []
        input_file.append('%mem={}MW'.format(int(config.memory * 1024 / 100))
        input_file.append("#P {}/{}".format(input_model.model.method, input_model.model.basis) + '\n')
        input_file.append("write your comment here\n")
  
        # Create a mol object
        mol = input_model.molecule
        input_file.append(f'{int(mol.molecular_charge)} {mol.molecular_multiplicity}')

        # Write the geometry
        for real, sym, geom in zip(mol.real, mol.symbols, mol.geometry):
            if real is False:
                raise InputError('Cannot handle ghost atoms yet.')
            input_file.append(f'{sym} {geom[0]:14.8f} {geom[1]:14.8f} {geom[2]:14.8f}')
        input_file.append("\n")

        gaussian_ret = {
            'infiles': {'input.inp': '\n'.join(input_file)},
            'commands': [which("g09"),  'input.inp', 'output.log'],
            'scratch_directory': config.scratch_directory
            }

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
        properties = {}
        
        tmp_output_path = outfiles['scratch_directory']
        tmp_output_file = os.path.join(tmp_output_path, 'output.log')
        data = cclib.io.ccread(tmp_output_file)
        
        last_occupied_energy = data.moenergies[0][data.homos[0]]
        output_data['HOMO ENERGY'] = last_occupied_energy
        #print (F'HOMO ENERGY: {last_occupied_energy:2.6f} eV')
        
        scf_energy = data.scfenergies[0]
        output_data['SCF ENERGY'] = scf_energy
        print (F'SCF ENERGY: {scf_energy:3.6f} eV')
        
        #if input_model.driver == 'energy':
           #output_data['return_result'] = 
        #print (os.system('ccget --list ' + tmp_output_file)) #data available in the output for parsing
        
        #if input_model.driver == 'energy':
        #   print (cclib.__version__)
        #   print (output_data)
        #print (input_model)
        
        #provenance = Provenance(creator="Gaussian", version=self.get_version(), routine="g09").dict()
        
        #properties = {
        #    'nuclear_repulsion_energy': Nuclear(data).repulsion_energy(),
        #    'scf_total_energy': data.scfenergies[0],
        #    'return_energy': data.scfenergies[0]
        #    }

        outtext = ''
        cnt = 0
        f = open(os.path.join(tmp_output_path, 'output.log'), 'r')
        outtext = f.readlines()
        for num, line in enumerate(outtext, 1):
           if 'Version=' in line:
              version_line = line.split('Version=')[-1]
              version_line = version_line.strip()
              version_line = version_line.split('\\')[0]
              print ('version : ', version_line)
              #cnt = num + 1
           #if cnt == num:
              #version_line += line.split("\\")[0]
              #print ('version is: ', version_line)
              

        provenance = Provenance(creator="Gaussian 09", version=self.get_version(), routine='g09').dict()
        print ('we are in mobj')
        mobj = re.search(r"^Job cpu time:*seconds.$", outfiles['output.log'])
        print ('we are in mobj')
        if (mobj):
           print ('mobj true')

        stdout = outfiles.pop('stdout')
        stderr = outfiles.pop('stderr')
        #print("\nPRINT STDOUT: \n", stdout)
        
        method = input_model.model.method.lower()
        #method = method[4:] if method.startswith("") else method

        #output_data["stdout"] = stdout
        output_data["success"] = True

        merged_data = {**input_model.dict(), **output_data}

        print('\nPRINT MERGED DATA: \n', merged_data)

        return AtomicResult(**merged_data)
