"""Compute quantum chemistry using Mainz-Austin-Budapest-Gainesville's CFOUR executable."""

import pprint
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple

import numpy as np

import qcelemental as qcel
from qcelemental.models import Provenance, Result
from qcelemental.util import which, safe_version

from ...util import execute
from ..model import ProgramHarness
from .harvester import harvest

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class CFOURHarness(ProgramHarness):

    _defaults = {
        "name": "CFOUR",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which('xcfour',
                     return_bool=True,
                     raise_error=raise_error,
                     raise_msg='Please install via http://cfour.de/')

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which('xcfour')
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "ZMAT"], {"ZMAT": "\nHe\n\n"})

            if success:
                for line in output["stdout"].splitlines():
                    if 'Version' in line:
                        branch = ' '.join(line.strip().split()[1:])
                self.version_cache[which_prog] = safe_version(branch)

        return self.version_cache[which_prog]

    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
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
            "command": [which("xcfour")],
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
            ["GRD", "FCMFINAL", "DIPOL"],
            scratch_messy=False,
            scratch_directory=inputs["scratch_directory"],
        )
        return success, dexe

    def parse_output(self, outfiles: Dict[str, str], input_model: 'ResultInput') -> 'Result':

        stdout = outfiles.pop("stdout")

        # c4mol, if it exists, is dinky, just a clue to geometry of cfour results
        qcvars, c4hess, c4grad, c4mol, version, errorTMP = harvest(input_model.molecule, stdout, **outfiles)

        if c4grad is not None:
            qcvars['CURRENT GRADIENT'] = c4grad

        if c4hess is not None:
            qcvars['CURRENT HESSIAN'] = c4hess

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
            'provenance': Provenance(creator="CFOUR", version=self.get_version(), routine="xcfour"),
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
