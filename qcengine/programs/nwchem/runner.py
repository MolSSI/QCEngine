"""
Calls the NWChem executable.
"""
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple

import numpy as np

import qcelemental as qcel
from qcelemental.models import Provenance, Result
from qcelemental.util import safe_version, which

from ..model import ProgramHarness
from ...util import execute
from .harvester import harvest


class NWChemHarness(ProgramHarness):

    _defaults = {
        "name": "NWChem",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": False,  # ATL: OpenMP only >=6.6 and only for Phi; potential for Mac using MKL and Intel compilers
        "node_parallel": True,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] ={}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool=False) -> bool:
        return which('nwchem', return_bool=True, raise_error=raise_error, raise_msg='Please install via http://www.nwchem-sw.org/index.php/Download')

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which('nwchem')
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "v.nw"], {"v.nw": ""})

            if success:
                for line in output["stdout"].splitlines():
                    if 'nwchem branch' in line:
                        branch = line.strip().split()[-1]
                    if 'nwchem revision' in line:
                        revision = line.strip().split()[-1]
                self.version_cache[which_prog] = safe_version(branch + '+' + revision)

        return self.version_cache[which_prog]

    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs NWChem in executable mode
        """
        self.found(raise_error=True)

        job_inputs = self.fake_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_model)

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:
        pass

    def fake_input(self, input_model: 'ResultInput', config: 'JobConfig',
                   template: Optional[str] = None) -> Dict[str, Any]:

        return {
            "command": [which("nwchem")],
            "infiles": input_model.extras['infiles'],
            "scratch_directory": config.scratch_directory,
            "input_result": input_model.copy(deep=True),
        }

    def execute(self,
                inputs: Dict[str, Any],
                *,
                extra_outfiles=None,
                extra_commands=None,
                scratch_name=None,
                timeout=None) -> Tuple[bool, Dict]:

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            ["job.movecs", "job.hess", "job.db", "job.zmat"],
            scratch_messy=False,
            scratch_directory=inputs["scratch_directory"],
        )
        return success, dexe

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':

        stdout = outfiles.pop("stdout")

        # nwmol, if it exists, is dinky, just a clue to geometry of nwchem results
#        qcvars, c4hess, c4grad, c4mol, version, errorTMP = harvest(input_model.molecule, stdout, **outfiles)
        #ORIGpsivar, nwhess, nwgrad, nwmol, version, errorTMP = harvester.harvest(qmol, nwchemrec['stdout'], **nwfiles)
        qcvars, nwhess, nwgrad, nwmol, version, errorTMP = harvest(input_model.molecule, stdout, **outfiles)

        if nwgrad is not None:
            qcvars['CURRENT GRADIENT'] = nwgrad

        if nwhess is not None:
            qcvars['CURRENT HESSIAN'] = nwhess

        retres = qcvars[f'CURRENT {input_model.driver.upper()}']
        if isinstance(retres, Decimal):
            retres = float(retres)
        elif isinstance(retres, np.ndarray):
            retres = retres.ravel().tolist()

        output_data = {
            'schema_name': 'qcschema_output',
            'schema_version': 1,
            'extras': {
                'outfiles': outfiles,
            },
            'properties': {},
            'provenance': Provenance(creator="NWChem", version=self.get_version(), routine="nwchem"),
            'return_result': retres,
            'stdout': stdout,
        }

        # got to even out who needs plump/flat/Decimal/float/ndarray/list
        # Decimal --> str preserves precision
        output_data['extras']['qcvars'] = {
            k.upper(): str(v) if isinstance(v, Decimal) else v
            for k, v in qcel.util.unnp(qcvars, flat=True).items()
        }

        output_data['success'] = True
        return Result(**{**input_model.dict(), **output_data})
