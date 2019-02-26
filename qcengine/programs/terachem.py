"""
Calls the TeraChem executable.
"""

from typing import Any, Dict, Optional

from qcelemental.models import Result

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
        properties = {}

        # Parse the output file, collect properties and gradient
        output_lines = open(outfiles["example.out"]).readlines()
        gradients = []
        for idx,line in enumerate(output_lines):
            if "FINAL ENERGY" in line:
                properties["scf_total_energy"] = float(line.strip('\n').split()[2])
                last_scf_line = output_lines[idx-2]
                properties["scf_iterations"] = int(last_scf_line.split()[0])
                if "XC Energy" in output_lines:
                    properties["scf_xc_energy"] = float(last_scf_line.split()[4])
            elif "DIPOLE MOMENT" in line:
                newline = line.replace(',','').replace('}','').replace('{','')
                properties["scf_dipole_moment"] = [ float(x) for x in newline.split()[2:5] ]
            elif "Nuclear repulsion energy" in line:
                properties["nuclear_repulsion_energy"] = float(line.split()[-2])
            elif "Gradient units are Hartree/Bohr" in line:
                #Gradient is stored as (dE/dx1,dE/dy1,dE/dz1,dE/dx2,dE/dy2,...)
                for i in range(idx+3,idx+3+natom):
                   grad = output_lines[i].strip('\n').split() 
                   for x in grad:
                       gradients.append( float(x) )
                       
        if len(gradients) > 0:
            output_data["return_result"] = gradients

        # Commented out the properties currently not supported by QCSchema
        #properites["spin_S2"] = 1 # calculated S(S+1)
        #   elif "SPIN S-SQUARED" in line:
        #       properties["spin_S2"] = float(line.strip('\n').split()[2])
        # Parse files in scratch folder
        #properties["atomic_charge"] = []
        #atomic_charge_lines =  open(outfiles["charge.xls"]).readlines()  
        #for line in atomic_charge_lines:
        #    properties["atomic_charge"].append(line.strip('\n').split()[-1]) 
        
        if "return_result" not in output_data:
            if "scf_total_energy" in properties:
                output_data["return_result"] = properties["scf_total_energy"]
            else:
                raise KeyError("Could not find SCF total energy")

        output_data["properties"] = properties

        output_data['schema_name'] = 'qcschema_output'
        # TODO Should only return True if TeraChem calculation terminated properly
        output_data['success'] = True

        return Result(**{**input_model.dict(), **output_data})
