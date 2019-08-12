"""Compute quantum chemistry using Iowa State's GAMESS executable."""

import pprint
from decimal import Decimal
from typing import Any, Dict, Optional

import qcelemental as qcel
from qcelemental.models import Result
from qcelemental.util import which, safe_version, unnp

from ...util import execute
from ..model import ProgramHarness
from .harvester import harvest

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class GAMESSHarness(ProgramHarness):

    _defaults = {
        "name": "GAMESS",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": True,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] ={}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool=False) -> bool:
        return which('rungms', return_bool=True, raise_error=raise_error, raise_msg='Please install via https://www.msg.chem.iastate.edu/GAMESS/GAMESS.html')

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which('rungms')
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "v.inp"], {"v.inp": ""})

            if success:
                for line in output["stdout"].splitlines():
                    if 'GAMESS VERSION' in line:
                        branch = ' '.join(line.strip(' *\t').split()[3:])
                self.version_cache[which_prog] = safe_version(branch)

        return self.version_cache[which_prog]

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        self.found(raise_error=True)

#        print("QCSCHEMA_INPUT")
#        print(input_data)

        #job_inputs = self.build_input(input_data, config)
        job_inputs = self.fake_input(input_data, config)
#        print('JOB_INPUTS')
#        pp.pprint(job_inputs)
        success, dexe = self.execute(job_inputs)
#        print('SUCCESS', success)
#        pp.pprint(dexe)

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_data)

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:
        pass

    def fake_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:

        # Note decr MEMORY=100000 to get
        # ***** ERROR: MEMORY REQUEST EXCEEDS AVAILABLE MEMORY
        # to test gms fail
        infile = \
""" $CONTRL SCFTYP=ROHF MULT=3 RUNTYP=GRADIENT COORD=CART $END
 $SYSTEM TIMLIM=1 MEMORY=800000 $END
 $SCF    DIRSCF=.TRUE. $END
 $BASIS  GBASIS=STO NGAUSS=2 $END
 $GUESS  GUESS=HUCKEL $END
 $DATA
Methylene...3-B-1 state...ROHF/STO-2G
Cnv  2

Hydrogen   1.0    0.82884     0.7079   0.0
Carbon     6.0
Hydrogen   1.0   -0.82884     0.7079   0.0
 $END
"""
        # edits to rungms
        # set SCR=./
        # set USERSCR=./
        # set GMSPATH=/home/psilocaluser/gits/gamess

        return {
            "commands": [which("rungms"), "gamess"],  # rungms JOB VERNO NCPUS >& JOB.log &
            "infiles": {
                "gamess.inp": infile
            },
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
        }


    def execute(self, inputs, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None):

        success, dexe = execute(inputs["commands"],
                                inputs["infiles"],
                                [],
                                scratch_messy=False,
                                scratch_directory=inputs["scratch_directory"],
                        )
        return success, dexe

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':

#        print('PARSE')
#        pp.pprint(outfiles)

        # gamessmol, if it exists, is dinky, just a clue to geometry of gamess results
        qcvars, gamessgrad, gamessmol = harvest(input_model.molecule, outfiles["stdout"]) #**gamessfiles)

        if gamessgrad is not None:
            qcvars['CURRENT GRADIENT'] = gamessgrad

        qcvars = unnp(qcvars, flat=True)

        output_data = {
            'schema_name': 'qcschema_output',
            'molecule' : gamessmol,
            'schema_version': 1,
            'extras': {},
            'properties': {
                'nuclear_repulsion_energy': gamessmol.nuclear_repulsion_energy(),
            },
            'return_result': qcvars[f'CURRENT {input_model.driver.upper()}'],
            'stdout': outfiles["stdout"],
        }

#        # Absorb results into psi4 data structures
#        for key in qcvars.keys():
#            core.set_variable(key.upper(), float(qcvars[key]))

#        if qcdbmolecule is None and gamessmol is not None:
#            molecule = geometry(gamessmol.create_psi4_string_from_molecule(), name='blank_molecule_psi4_yo')
#            molecule.update_geometry()
#            # This case arises when no Molecule going into calc (cfour {} block) but want
#            #   to know the orientation at which grad, properties, etc. are returned (c4mol).
#            #   c4mol is dinky, w/o chg, mult, dummies and retains name
#            #   blank_molecule_psi4_yo so as to not interfere with future cfour {} blocks

        # got to even out who needs plump/flat/Decimal/float/ndarray/list
        output_data['extras']['qcvars'] = {k.upper(): float(v) if isinstance(v, Decimal) else v for k, v in qcel.util.unnp(qcvars, flat=True).items()}

#        # Quit if gamess threw error
#        if core.get_variable('GAMESS ERROR CODE'):
#            raise ValidationError("""gamess exited abnormally.""")
#
#        P4GAMESS_INFO.clear()
#        P4GAMESS_INFO.update(internal_p4gamess_info)
#
#        optstash.restore()

        output_data['success'] = True

        return Result(**{**input_model.dict(), **output_data})

#            return Result(**output_data)
#        return FailedOperation(success=output_data.pop("success", False), error=output_data.pop("error"), input_data=output_data)
