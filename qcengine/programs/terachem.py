"""
Calls the TeraChem executable.
"""

from typing import Any, Dict, Optional

from qcelemental.models import ComputeError, FailedOperation, Provenance, Result
import os

from ..units import ureg
from .executor import ProgramExecutor


class TeraChemExecutor(ProgramExecutor):

    _defaults = {
         "name": "TeraChem",
         "scratch": True,
         "thread_safe": False,
         "thread_parallel": True,
         "node_parallel": False,
         "managed_memory": True,
    }

    class Config(ProgramExecutor.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs}) 

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        pass

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig', 
                    template: Optional[str]=None) -> Dict[str, Any]:
        #Write the geom xyz file
        xyz_file = [] 
        s = str(len(input_model.molecule.symbols))  
        xyz_file.append(s+"\n")
        for sym, geom in zip(input_model.molecule.symbols, input_model.molecule.geometry):
            s = "{:<4s} {:>{width}.{prec}f} {:>{width}.{prec}f} {:>{width}.{prec}f}".format(sym, *geom, width=14, prec=10)
            xyz_file.append(s)
  
        xyz_file = "\n".join(xyz_file)

        # Write input file
        input_file = []
        input_file.append("# molecule definition")
        input_file.append( "charge " + str(int(input_model.molecule.molecular_charge)))
        input_file.append( "spinmult " + str(input_model.molecule.molecular_multiplicity))
        input_file.append( "coordinates geometry.xyz")
        
        input_file.append("\n# model")
        input_file.append("basis " + str(input_model.model.basis))
        
        
        input_file.append("\n# driver")
        input_file.append("run " +input_model.driver)
        
        input_file.append("\n# keywords")
        for k, v in input_model.keywords.items():
            input_file.append("{} {}".format(k, v))
  
        return {
            "commands": ["terachem","example.in"],
            "infiles": {
                "example.in": input_file,
                "geometry.xyz": xyz_file
            },
            "input_result": input_model.copy(deep=True)
        }

        def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':
            output_data = {}
            name_space = {}
            properties = {}
            # Parse the output file
            output_lines = open(outfiles["example.out"]).readlines()
            output_data["energy"] = 0
            output_data["gradients"] = []
            output_data["spin_S2"] = 1 # calculated S(S+1)
            for idx,line in enumerate(output_lines):
                if "FINAL ENERGY" in line:
                    output_data["energy"] = float(line.strip('\n').split()[2])
                elif "Gradient units are Hartree/Bohr" in line:
                    #Gradient is stored as (dE/dx1,dE/dy1,dE/dz1,dE/dx2,dE/dy2,...)
                    for i in range(idx+3,idx+3+natom):
                       grad = output_lines[i].strip('\n').split() 
                       for x in grad:
                           output_data["gradients"].append( float(x) )
                elif "SPIN S-SQUARED" in line:
                    output_data["spin_S2"] = float(line.strip('\n').split()[2])

            # Parse files in scratch folder
            output_data["atomic_charge"] = []
            atomic_charge_lines =  open(outfiles["charge.xls"]).readlines()  
            for line in atomic_charge_lines:
                output_data["atomic_charge"].append(line.strip('\n').split()[-1]) 
            
            output_data["properties"] = properties

            output_data['schema_name'] = 'qcschema_output'
            # TODO Should only return True if TeraChem calculation terminated properly
            output_data['success'] = True

            return Result(**{**input_model.dict(), **output_data})
