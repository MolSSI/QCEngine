import re, os

import cclib

from typing import Union, Dict, Any, List, Tuple

from qcelemental.models import (OptimizationInput, OptimizationResult,
                                AtomicInput, AtomicResult,
                                Molecule, Provenance)

from qcengine.exceptions import UnknownError, InputError
from qcengine.procedures.model import ProcedureHarness
from qcengine.config import TaskConfig
from qcengine.programs.gaussian import GaussianHarness
from qcengine.programs.util import PreservingDict, error_stamp

class GaussianDriverProcedure(ProcedureHarness):
    '''
    This makes a structure for geometry optimization in Gaussian.
    '''
    
    _defaults = {'name': 'GaussianDriver', 'procedure': 'optimization'}
    
    class Config(ProcedureHarness.Config):
        pass
    
    def found(self, raise_error: bool = False) -> bool:
        ga_harness = GaussianHarness()
        return ga_harness.found(raise_error)
        
    def get_version(self) -> str:
        ga_harness = GaussianHarness()
        return ga_harness.get_version()
    
    def build_input_model(self, input_data: Union[Dict[str, Any], 'OptimizationInput']) -> OptimizationInput:
        return self._build_model(input_data, OptimizationInput)
    
    
    def compute(self, input_data: OptimizationInput, config: TaskConfig) -> 'BaseModel':
        ga_harness = GaussianHarness()
        self.found(raise_error = True)
        
        keywords = input_data.keywords.copy()
        keywords.update(input_data.input_specification.keywords)
        
        if keywords.get('program', 'gaussian').lower() != 'gaussian' :
            raise InputError('GaussianDriver procedure only works with Gaussian software package.')
        
        input_data.input_specification.extras['is_driver'] = True
        
        # Make an atomic input
        atomic_input = AtomicInput(
            molecule = input_data.initial_molecule,
            driver = 'energy',
            keywords = keywords,
            **input_data.input_specification.dict(exclude = {'driver', 'keywords'}),
        )
        
        # Build an input
        job_inputs = ga_harness.build_input(atomic_input, config, 'opt')
        stdin = job_inputs['infiles']['input.inp']
        #print ('STDIN:', stdin)
                
        # Run the input
        success, _exe = ga_harness.execute(job_inputs)
        
        # Determine whether the calculation succeeded
        if success:
        	# If execution succeeded, collect results
            return self.parse_opt_geom_output(_exe, atomic_input) 
        else:
            _exe['outfiles']['stderr'] = _exe['outfiles']['output.log']
            outfile = _exe['outfiles']['output.log']
            
            if 'Error termination via ' in outfile:
                raise InputError(error_stamp(stdin, outfile))
         
            else:
                # Return UnknownError for error propagation
                raise UnknownError(error_stamp(stdin, outfile))
       
    def parse_opt_geom_output(self, outfiles: Dict[str, str], input_model: 'OptimizationInput') -> 'OptimizationResult':

        ga_harness = GaussianHarness()
        
        tmp_output_path = outfiles['scratch_directory']
        #print ('tmp_output_path=',tmp_output_path)
        tmp_output_file = os.path.join(tmp_output_path, 'output.log')
        #print ('tmp_output_file=',tmp_output_file)
        data = cclib.io.ccread(tmp_output_file)
        #data = parser.parse()
        print ('There are %i atoms and %i MO' %(data.natom, data.nmo))
        print ('\nOptimized geometries at each step:\n')
        print (data.atomcoords)
       
       # Get the stdout from the calculation (required)
        stdout = outfiles.pop("stdout")
        stderr = outfiles.pop("stderr")

        return OptimizationResult(
            initial_molecule = input_model.initial_molecule,
            input_specification = input_model.input_specification,
            final_molecule = '', #final_step.molecule,
            trajectory = atomic_results,
            energies = [float(r.extras["qcvars"]["CURRENT ENERGY"]) for r in atomic_results],
            stdout = stdout,
            stderr = stderr,
            success = True,
            #provenance = Provenance(creator= '', version = self.get_version(), routine = 'gaussian_opt'),
        )
    
 
