"""
Calls the GAUSSIAN executable
"""

import os
import re
import tempfile
import warnings
import pprint
import copy

from collections import defaultdict

from typing import Any, Dict, List, Optional, Tuple

import cclib
from cclib.method import Nuclear

from qcelemental import constants
from qcelemental.models import AtomicInput, AtomicResult, Molecule, Provenance
from qcelemental.models import OptimizationInput, OptimizationResult, BasisSet
from qcelemental.molparse import regex
from qcelemental.util import parse_version, safe_version, which, which_import

from ..exceptions import InputError, UnknownError
from ..util import disk_files, execute, temporary_directory
from .util import error_stamp
from .model import ProgramHarness

from qcengine.config import TaskConfig, get_config
from qcengine.procedures.model import ProcedureHarness

class GaussianHarness(ProgramHarness):

    _defaults = {
        "name": "gaussian",
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
        '''
        It checks to see if GaussianHarness is ready for operation
        
        Parameters
        ----------
        raise_error: bool
           Passed on to control negative return between False and ModuleNotFoundError raised.

        Returns
        -------
        bool
           If both gaussian and its dependency cclib are found, returns True
           If gaussian or cclib are missing, returns False.
           If raise_error is True and gaussian or cclib are missing, the error message for the missing one is raised.
        '''
        qc = which(
            "g16",
            return_bool = True,
            raise_error = False,
            raise_msg = "Please install Gaussian. Check it's in your PATH with `which g16`."
        )

        if not qc:
            qc = which(
            "g09",
            return_bool = True,
            raise_error = raise_error,
            raise_msg = "Please install Gaussian. Check it's in your PATH with `which g09` or `which g16`."
        )
        
        dep = which_import(
                "cclib",
                return_bool = True,
                raise_error = raise_error,
                raise_msg = "For gaussian harness, please install cclib by typing in `conda install -c conda-forge cclib`."
                )
        
        return qc & dep

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("g16")
        if which_prog is None:
            which_prog = which('g09')
        
        self.version_cache[which_prog] = safe_version(which_prog.split('/')[-1])
        
        return self.version_cache[which_prog]

    def compute(self, input_model: 'AtomicInput', config: 'TaskConfig') -> 'AtomicResult':
        """
        Run Gaussian program
        """
        
        # Check if Gaussian executable is found
        self.found(raise_error=True)
        
        # Setup the job
        job_inputs = self.build_input(input_model, config)
        
        # Run Gaussian
        success, _exe = self.execute(job_inputs)

        stdin = job_inputs['infiles']['input.inp']

        if isinstance(input_model.model.basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use string basis name.")

        # Determine whether the calculation succeeded
        if success:
            # If execution succeeded, collect results
            return self.parse_output(_exe, input_model)
        else:
            _exe['outfiles']['stderr'] = _exe['outfiles']['output.log']
            outfile = _exe['outfiles']['output.log']
            
            if 'Error termination via ' in outfile:
                raise InputError(error_stamp(stdin, outfile))
         
            else:
                # Return UnknownError for error propagation
                raise UnknownError(error_stamp(stdin, outfile))

    def build_input(
            self, input_model: AtomicInput,
            config: TaskConfig,
            template: Optional[str] = None) -> Dict[str, Any]:

        if template is None:
            input_file = []
            caseless_keywords = {k.lower(): v for k, v in input_model.keywords.items()}

        # Build keywords
        keywords = {k.upper(): v for k, v in input_model.keywords.items()}

        gaussian_kw = []

        if input_model.driver == "energy":
            gaussian_kw.append("sp")
        elif input_model.driver == "gradient":
            gaussian_kw.append("force")
        elif input_model.driver == "hessian":
            gaussian_kw.append("freq")
        else:
            raise InputError(f"Driver {input_model.driver} not implemented for Gaussian.")

        #if input_model.molecule.fix_com or input_model.molecule.fix_orientation:
        #    keywords["SYM_IGNORE"] = "TRUE"
        if 'SCF_CONVERGENCE' in keywords:
           gaussian_kw.append('SCF=' + keywords["SCF_CONVERGENCE"])
        if 'POPULATION' in keywords:
           gaussian_kw.append('Pop=' + keywords['POPULATION'])   
           
        keywords = {'scf_damp': 'true',
                    'scf_diis': 'false'}

        save_fchk = False

        if input_model.protocols.native_files == 'all':
            save_fchk = True
            gaussian_kw.append('formcheck')

        # Begin input file
        input_file = []
        input_file.append('%mem={}MB'.format(int(config.memory * 1024))) # In MB
        input_file.append('#P {}/{}'.format(input_model.model.method, input_model.model.basis) + ' ' + ' '.join(gaussian_kw) + '\n')
        input_file.append('write your comment here\n')
  
        # Handle the geometry
        molcmd, moldata = input_model.molecule.to_string(dtype = 'gaussian', units = 'Angstrom', return_data = True)
        input_file.append(molcmd.lstrip())

#        print ('*' * 100)
#        print ('\n'.join(input_file))
#        print ('*' * 100)

        gaussian_ret = {
            'infiles': {'input.inp': '\n'.join(input_file)},
            'commands': [which("g09"),  'input.inp', 'output.log'],
            'scratch_directory': config.scratch_directory,
            'scratch_messy': config.scratch_messy
        }
        
        return gaussian_ret

    def execute(self,
                inputs,
                extra_infiles: Optional[Dict[str, str]] = None,
                extra_outfiles: Optional[Dict[str, str]] = None,
                extra_commands: Optional[List[str]] = None,
                scratch_name = None,
                timeout: Optional[int] = None
               ):

        success, dexe = execute(
            inputs['commands'],
            inputs['infiles'],
            outfiles = ['output.log'],
            scratch_directory = inputs['scratch_directory'],
            scratch_messy = True #inputs['scratch_messy']
            )
        
        return success, dexe
 
    def parse_output(self, outfiles: Dict[str, str], input_model: 'AtomicInput') -> 'AtomicResult':

        tmp_output_path = outfiles['scratch_directory']
        print ('tmp_output_path=',tmp_output_path)
        tmp_output_file = os.path.join(tmp_output_path, 'output.log')
        print ('tmp_output_file=',tmp_output_file)
        data = cclib.io.ccread(tmp_output_file)
        print ('DATA=', data)
        scf_energy = data.scfenergies[0] / constants.conversion_factor("hartree", "eV") # Change from the eV unit to the Hartree unit
        print ('SCF ENERGY: ', scf_energy)
        
        output_data = {
        	"schema_version": 1,
        	"molecule": input_model.molecule,
        	'extras': input_model.extras,
            "native_files": {k: v for k, v in outfiles.items() if v is not None},
            'properties': '',
            'provenance': Provenance(creator="gaussian", version=self.get_version(), routine='g09').dict(),
            'return_result': scf_energy,
            "stderr": outfiles['outfiles']['output.log'],,
            "stdout": outfiles['outfiles']['output.log'],
            "success": True
        		}

        #output_data = {}
        #properties = {}
        #cclib_vars = {}
        
        #tmp_output_path = outfiles['scratch_directory']
        #tmp_output_file = os.path.join(tmp_output_path, 'output.log')
        #data = cclib.io.ccread(tmp_output_file)
        cclib_vars = data.getattributes(True)
        
        #last_occupied_energy = data.moenergies[0][data.homos[0]]
        #output_data['HOMO ENERGY'] = last_occupied_energy
        
        #scf_energy = data.scfenergies[0] / constants.conversion_factor("hartree", "eV") # Change from the eV unit to the Hartree unit
        #output_data['SCF ENERGY'] = scf_energy
        
        #if input_model.driver == 'energy':
        #    output_data['return_result'] = scf_energy
        #if input_model.driver == 'gradient':
        #    output_data['return_result'] = data.grads
        #    #output_data['return_gradient'] = data.grads
        
        #print (os.system('ccget --list ' + tmp_output_file)) #data available in the output for parsing

        #if input_model.driver == 'energy':
        #   print (cclib.__version__)        
        #   print (output_data)
        #print (input_model)

        properties = {
            'nuclear_repulsion_energy': Nuclear(data).repulsion_energy()
        }


        if input_model.model.method.lower() in ['hf', 'scf']:
            scf_energy = data.scfenergies[0] / constants.conversion_factor("hartree", "eV")
            properties['scf_total_energy'] = scf_energy
            properties['return_energy'] = scf_energy
            output_data['return_result'] = scf_energy
        
        if input_model.model.method.lower().startswith('mp'):
            scf_energy = data.scfenergies[0] / constants.conversion_factor("hartree", "eV") # Change from the eV unit to the Hartree unit
            mp2_energy = data.mpenergies[0] / constants.conversion_factor("hartree", "eV")
            properties['scf_total_energy'] = scf_energy
            properties['mp2_total_energy'] = mp2_energy[0]
            properties['return_energy'] = mp2_energy[0]
            output_data['return_result'] = mp2_energy[0]

        if input_model.model.method.lower().startswith('cc'):
            scf_energy = data.scfenergies[0] / constants.conversion_factor("hartree", "eV")
            cc_energy = data.ccenergies[0] / constants.conversion_factor("hartree", "eV")
            properties['scf_total_energy'] = scf_energy
            properties['ccsd_prt_pr_total_energy'] = cc_energy
            properties['return_energy'] = cc_energy
            output_data['return_result'] = cc_energy           

        properties['calcinfo_nbasis'] = data.nbasis
        properties['calcinfo_nmo'] = data.nmo
        properties['calcinfo_natom'] = data.natom
        
        output_data['properties'] = properties
        #output_data['stdout'] = outfiles['outfiles']['output.log']
        #output_data['success'] = True
        ##print ('output_data: ', output_data)

        #stdout = outfiles.pop('stdout')
        #stderr = outfiles.pop('stderr')
        #print("\nPRINT STDOUT: \n", stdout)
        

        #method = input_model.model.method.lower()
        #method = method[4:] if method.startswith("") else method

        # filter unwanted data
        #to_remove = ['atomnos', 'atomcoords', 'natom']
        #output_data['extras'] = {'cclib': {k:v for k, v in cclib_vars.items() if k not in to_remove}}

        # HACK - scf_values can have NaN
        # remove for now
        #output_data['extras']['cclib'].pop('scfvalues')

        if (input_model.protocols.native_files == 'all'):
            output_data['native_files'] = {}
            tmp_input_file = os.path.join(tmp_output_path, 'input.inp')
            
            with open(tmp_input_file, 'rt') as f:
                output_data['native_files']['input.inp'] = f.read()
                
            # formcheck keyword always creates the Test.FChk file
            tmp_check_file = os.path.join(tmp_output_path, 'Test.FChk')
            
            with open(tmp_check_file, 'rt') as f:
                output_data['native_files']['gaussian.fchk'] = f.read()
                
        merged_data = {**input_model.dict(), **output_data}

        return AtomicResult(**merged_data)

class GaussianDriverProcedure(ProcedureHarness):
    '''
    Optimize a geometry
    '''

    _defaults = {'name': 'GaussianDriver', 'procedure': 'optimization'}

    class Config(ProcedureHarness.Config):
        pass

    def found(self, raise_error: bool = False) -> bool:
        gaussian_harness = GaussianHarness()
        return gaussian_harness(raise_error)
