"""Compute quantum chemistry using Iowa State's GAMESS executable."""

import copy
import pprint
import re
from decimal import Decimal
from typing import Any, Dict, Optional

import numpy as np
import qcelemental as qcel
from qcelemental.models import FailedOperation, Molecule, Result
from qcelemental.util import which, safe_version

from ...util import execute
from ..executor import ProgramExecutor
#from ...extras import provenance_stamp
from ..util import PreservingDict

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class GAMESSExecutor(ProgramExecutor):

    _defaults = {
        "name": "GAMESS",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": True,
        "managed_memory": True,
    }

    class Config(ProgramExecutor.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @staticmethod
    def found(raise_error=False) -> bool:
        is_found = which('rungms', return_bool=True)

        if not is_found and raise_error:
            raise ImportError("Could not find `rungms` for GAMESS in the shell path.")
        else:
            return is_found

#    def get_version(self) -> str:
#        self.found(raise_error=True)
#
#        # Note: anything below v3.2.1 will return the help menu here. but that's fine as version compare evals to False.
#        command = [which('dftd3'), '-version']
#        import subprocess
#        proc = subprocess.run(command, stdout=subprocess.PIPE)
#        candidate_version = proc.stdout.decode('utf-8').strip()
#
#        return safe_version(candidate_version)


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
                #"gamess.inp": infile,
                "gamess.inp": input_model.extras['gamess.inp'],
            },
            "scratch_location": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
        }


    def execute(self, inputs, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None):

        success, dexe = execute(inputs["commands"],
                                inputs["infiles"],
                                [],
                                scratch_messy=True,
                                scratch_location=inputs["scratch_location"],
                        )
        return success, dexe

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':

#        print('PARSE')
#        pp.pprint(outfiles)

        # gamessmol, if it exists, is dinky, just a clue to geometry of gamess results
        qcvars, gamessgrad, gamessmol = harvest(input_model.molecule, outfiles["stdout"]) #**gamessfiles)

        output_data = {
            'schema_name': 'qcschema_output',
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

        if gamessgrad is not None:
            qcvars['CURRENT GRADIENT'] = gamessgrad

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

        # TODO Should only return True if Molpro calculation terminated properly
        output_data['success'] = True

        return Result(**{**input_model.dict(), **output_data})

#            return Result(**output_data)
#        return FailedOperation(success=output_data.pop("success", False), error=output_data.pop("error"), input_data=output_data)


def harvest(p4Mol, gamessout, **largs):
    """Parses all the pieces of output from gamess: the stdout in
    *gamessout* Scratch files are not yet considered at this moment.
    """
    outqcvar, outMol, outGrad = harvest_output(gamessout)

    if outMol:
        if p4Mol:
            if abs(outMol.nuclear_repulsion_energy() - p4Mol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValueError("""gamess outfile (NRE: %f) inconsistent with Psi4 input (NRE: %f).""" % \
                            (outMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy()))
    else:
        raise ValueError("""No coordinate information extracted from gamess output.""")

    return outqcvar, outGrad, outMol


def harvest_output(outtext):
    """Function to separate portions of a gamess output file *outtext*,
        divided by "Step".
        """
    pass_qcvar = []
    pass_coord = []
    pass_grad = []

    for outpass in re.split(
        r'^\s+'+r'--------'+
        r'NSERCH:'+r'([1-9][0-9][0-9][0-9]*)'+r'\s*'+
        r'^\s+'+r'--------', outtext, re.MULTILINE):

        qcvar, gamesscoord, gamessgrad = harvest_outfile_pass(outpass)
        pass_qcvar.append(qcvar)
        pass_coord.append(gamesscoord)
        pass_grad.append(gamessgrad)

    retindx = -1 if pass_coord[-1] else -2
    return pass_qcvar[retindx], pass_coord[retindx], pass_grad[retindx]


def harvest_outfile_pass(outtext):
    """Function to read gamess output file *outtext* and parse important
    quantum chemical information from it in
    """
    qcvar = PreservingDict()
    qcvar_coord = None
    qcvar_grad = None

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    # If calculation fail to converge
    mobj = re.search(
        r'^\s+' + r'(?:GAMESS TERMINATED ABNORMALLY)'+r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('GAMESS TERMINATED ABNORMALLY')

    # If calculation converged
    else:
        mobj = re.search(
            r'^\s+' + r'(?:            TOTAL ENERGY)' + r'\s+=\s*' + NUMBER +r's*$'
            ,outtext, re.MULTILINE)
        if mobj:
            print('matched gamess_RHF energy')
            qcvar['HF TOTAL ENERGY'] = mobj.group(1)

        # Process NRE
        mobj = re.search(
            r'^\s+' + r'(?:   NUCLEAR REPULSION ENERGY)' + r'\s+=\s*' + NUMBER +r'\s*$'
            ,outtext, re.MULTILINE)
        if mobj:
            print('matched NRE')
            qcvar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)

        # Process MP2
        mobj = re.search(
            r'^\s+' + r'E\(0\)' + r'=\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'E\(1\)' + r'=\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'E\(2\)' + r'=\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'E\(MP2\)' + r'=\s+' + NUMBER + r'\s*$'
            ,outtext, re.MULTILINE)
        if mobj:
            print('matched mp2')
            qcvar['MP2 CORRELATION ENERGY'] = mobj.group(3)
            qcvar['MP2 TOTAL ENERGY'] = mobj.group(4)

        # Process CCSD(T)
        mobj = re.search(
            r'^\s+' + r'SUMMARY OF RESULTS' + 
            r'\s+' + r'\n'+
            r'^\s+' + r'REFERENCE ENERGY:' + r'\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'MBPT\(2\) ENERGY:' + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'CCSD    ENERGY:'   + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'CCSD\[T\] ENERGY:' + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*'+
            r'^\s+' + r'CCSD\(T\) ENERGY:' + r'\s+' + NUMBER + r'\s*' + r'CORR.E=\s+' + r'\s+' + NUMBER + r'\s*$'
            ,outtext, re.MULTILINE)
        if mobj:
            print('matched ccsd(t)')
            qcvar['SCF TOTAL ENERGY'] = mobj.group(1)
            qcvar['CCSD CORRELATION ENERGY'] = mobj.group(5)
            qcvar['CCSD TOTAL ENERGY'] = mobj.group(4)
            qcvar['(T) CORRECTION ENERGY'] = Decimal(mobj.group(8)) - Decimal(mobj.group(4)) 
            #qcvar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(5)) - qcvar['SCF TOTAL ENERGY'] 
            qcvar['CCSD(T) CORRELATION ENERGY'] = mobj.group(9)
            qcvar['CCSD(T) TOTAL ENERGY'] = mobj.group(8)

        # Process FCI
        mobj = re.search(
            r'^\s+' + r'ALDET CI PROPERTIES...FOR THE WAVEFUNCTION OF STATE' + r'\s+' + r'(?P<state>\d+)' + r'\s*' +
            r'^\s+' + r'USING THE EXPECTATION VALUE DENSITY' + r'\s*' +
            r'(?:.*?)' +
            r'^\s+' + r'TOTAL ENERGY =' + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
        if mobj:
            print('matched fci')
            qcvar['FCI TOTAL ENERGY'] = mobj.group(2)
            qcvar['CI TOTAL ENERGY'] = mobj.group(2)

        # Process EOM-CCSD
        mobj = re.search(
            r'^\s+' + r'---- SUMMARY OF EOM-CCSD CALCULATIONS ----' + 
            r'\s+' + r'\n' + 
            r'^\s+' + r'EXCITATION      EXCITATION      TOTAL' + r'\s*' +
            r'^\s+' + r'SYMMETRY     ENERGY \(H\)      ENERGY \(EV\)     ENERGY \(H\)          ITERATIONS' + r'\s*' +
            r'^\s+' + r'A' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + r'CONVERGED\s+' +
            r'^\s+' + r'A' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + r'CONVERGED\s+' +
            r'^\s+' + r'A' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + r'CONVERGED\s+' +
            r'^\s+' + r'A' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + r'CONVERGED\s+'
            ,outtext, re.MULTILINE)
        if mobj:
            print('matched eom-ccsd')
            qcvar['EOM-CCSD ROOT 1 EXCIATION ENERGY'] =  mobj.group(1)
            qcvar['EOM-CCSD ROOT 2 EXCIATION ENERGY'] =  mobj.group(4)
            qcvar['EOM-CCSD ROOT 3 EXCIATION ENERGY'] =  mobj.group(7)
            qcvar['EOM-CCSD ROOT 4 EXCIATION ENERGY'] =  mobj.group(10)
            qcvar['EOM-CCSD ROOT 1 TOTAL ENERGY'] =  mobj.group(3)
            qcvar['EOM-CCSD ROOT 2 TOTAL ENERGY'] =  mobj.group(6)
            qcvar['EOM-CCSD ROOT 3 TOTAL ENERGY'] =  mobj.group(9)
            qcvar['EOM-CCSD ROOT 4 TOTAL ENERGY'] =  mobj.group(12)

        # Process DFT
        # mobj = re.search(
        #     r'^\s+' + r'DFT EXCHANGE + CORRELATION ENERGY' + r'=\s+' + NUMBER + r'\s*$'
        #     ,outtext, re.MULTILINE)
        mobj = re.search(
            r'^\s+' + r'DFT EXCHANGE \+ CORRELATION ENERGY' + r'\s+=\s*' + NUMBER +r'\s*$'
            ,outtext, re.MULTILINE)
        if mobj:
            print('matched dft')
            qcvar['DFT TOTAL ENERGY'] = mobj.group(1)

        # Process Geometry
        mobj = re.search(
                r'^\s+' + r'ATOM      ATOMIC                      COORDINATES \(BOHR\)' + r'\s*' +
                r'^\s+' + r'CHARGE         X                   Y                   Z'+ r'\s*' +
                r'((?:\s+([A-Z][a-z]*)+\s+\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)'+r'\s*$'
                ,outtext,re.MULTILINE | re.IGNORECASE)
        if mobj:
                print('matched geom')
                molxyz = '%d bohr\n\n' % len(mobj.group(1).splitlines())
                for line in mobj.group(1).splitlines():
                        lline = line.split()
                        molxyz += '%s %16s %16s %16s\n' % (int(float(lline[-4])), lline[-3], lline[-2], lline[-1])
                qcvar_coord = Molecule(validate=False, **qcel.molparse.to_schema(qcel.molparse.from_string(molxyz, dtype='xyz+', fix_com=True, fix_orientation=True)["qm"], dtype=2))

        # Process Gradient
        mobj = re.search(
                r'^\s+' + r'GRADIENT OF THE ENERGY' + r'\s*'+
                r'^\s+' + r'----------------------' + r'\s*'+
                r'\s+' + r'\n'+
                r'^\s+' + r'UNITS ARE HARTREE/BOHR    E\'X               E\'Y               E\'Z' + r'\s*'+
                r'((?:\s+([1-9][0-9]*)+\s+([A-Z][a-x]*)+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)'+r'\s*$'
                ,outtext, re.MULTILINE)
        if mobj:
                print('matched gradient - after')
                atoms = []
                qcvar_grad = []
                for line in mobj.group(1).splitlines():
                        lline = line.split()
                        if lline == []:
                            pass
                        else:
                            print('printing gradient')
                            atoms.append(lline[1])
                            qcvar_grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])
                qcvar_grad = np.array(qcvar_grad)

    # Process CURRENT Energies
    if 'HF TOTAL ENERGY' in qcvar:
        qcvar['CURRENT REFERENCE ENERGY'] = qcvar['HF TOTAL ENERGY']
        qcvar['CURRENT ENERGY'] = qcvar['HF TOTAL ENERGY']

    if 'MP2 TOTAL ENERGY' in qcvar and 'MP2 CORRELATION ENERGY' in qcvar:
            qcvar['CURRENT CORRELATION ENERGY'] = qcvar['MP2 CORRELATION ENERGY']
            qcvar['CURRENT ENERGY'] = qcvar['MP2 TOTAL ENERGY']

    if 'CCSD(T) TOTAL ENERGY' in qcvar and 'CCSD(T) CORRELATION ENERGY' in qcvar:
        qcvar['CURRENT CORRELATION ENERGY'] = qcvar['CCSD(T) CORRELATION ENERGY']
        qcvar['CURRENT ENERGY'] = qcvar['CCSD(T) TOTAL ENERGY']

    if 'DFT TOTAL ENERGY' in qcvar:
        qcvar['CURRENT REFERENCE ENERGY'] = qcvar['DFT TOTAL ENERGY']
        qcvar['CURRENT ENERGY'] = qcvar['DFT TOTAL ENERGY']

    if 'FCI TOTAL ENERGY' in qcvar: # and 'FCI CORRELATION ENERGY' in qcvar:
        qcvar['CURRENT ENERGY'] = qcvar['FCI TOTAL ENERGY']

    return qcvar, qcvar_coord, qcvar_grad
