"""Compute dispersion correction using Greenwell & Beran's MP2D executable."""

import copy
import os
import pprint
import re
import sys
import traceback
from decimal import Decimal

import numpy as np
import qcelemental as qcel
from qcelemental.models import FailedOperation, Result

from .dftd3 import dashparam
from .dftd3.runner import module_driver  # nasty but temporary and better than duplicating fn
from ..util import which
from .executor import ProgramExecutor
from ..extras import provenance_stamp

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class MP2DExecutor(ProgramExecutor):

    _defaults = {
        "name": "MP2D",
        "scratch": True,
        "thread_safe": True,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": False,
    }

    class Config(ProgramExecutor.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @staticmethod
    def found(raise_error=False) -> bool:
        is_found = which('mp2d', return_bool=True)

        if not is_found and raise_error:
            raise ImportError("Could not find MP2D in the shell path.")
        else:
            return is_found

    def get_version(self) -> str:
        self.found(raise_error=True)

        # Note: no version at present. Need to get Chandler to set one up
        #command = [which('mp2d'), '-version']
        #import subprocess
        #proc = subprocess.run(command, stdout=subprocess.PIPE)
        #candidate_version = proc.stdout.decode('utf-8').strip()
        candidate_version = '0.1'

        from pkg_resources import safe_version
        return safe_version(candidate_version)

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        self.found(raise_error=True)

        # Set up the job
        input_data = input_data.copy().dict()
        input_data["success"] = False

        output_data = run_json(input_data)

        if output_data["success"]:
            return Result(**output_data)
        return FailedOperation(
            success=output_data.pop("success", False), error=output_data.pop("error"), input_data=output_data)


def run_json(jobrec):
    """
    An implementation of the QC JSON Schema (molssi-qc-schema.readthedocs.io/en/latest/index.html#) for MP2D (from Psi4).

    Parameters
    ----------
    mjobrec : qcelemental.models.ResultInput
    jobrec : JSON
        Please see molssi-qc-schema.readthedocs.io/en/latest/spec_components.html for further details.

    """
    jobrec['provenance'] = provenance_stamp(sys._getframe().f_code.co_name + '.' + __name__)

    # strip engine hint
    mtd = jobrec['model']['method']
    if mtd.startswith('mp2d-'):
        jobrec['model']['method'] = mtd[5:]

    if jobrec['driver'].derivative_int() > 1:
        jobrec['success'] = False
        jobrec['error'] = {
            'error_type': 'ValueError',
            'error_message': """MP2D produces max gradient, not {jobrec['driver']}"""
        }
        raise ValueError("""MP2D produces max gradient, not {jobrec['driver']}""")

    try:
        mp2d_driver(jobrec)
    except Exception as exc:
        jobrec['success'] = False
        jobrec['error'] = {
            'error_type': type(exc).__name__,
            'error_message': ''.join(traceback.format_exception(*sys.exc_info())),
        }
        raise exc

    jobrec['success'] = True

    ene = jobrec["extras"]['qcvars']['DISPERSION CORRECTION ENERGY']
    jobrec["extras"]["qcvars"]["CURRENT ENERGY"] = ene
    jobrec['properties'] = {"return_energy": ene}

    if jobrec['driver'] == 'energy':
        jobrec["return_result"] = ene
    elif jobrec['driver'] == 'gradient':
        grad = copy.deepcopy(jobrec["extras"]['qcvars']['DISPERSION CORRECTION GRADIENT'])
        jobrec["extras"]['qcvars']['CURRENT GRADIENT'] = grad
        jobrec["return_result"] = grad

    jobrec["molecule"]["real"] = list(jobrec["molecule"]["real"])

    return jobrec


def mp2d_driver(jobrec, verbose=1):
    """Drive the jobrec@i (input) -> mp2drec@i -> mp2drec@io -> jobrec@io (returned) process."""

    return module_driver(
        jobrec=jobrec, module_label='mp2d', plant=mp2d_plant, harvest=mp2d_harvest, verbose=verbose)


def mp2d_plant(jobrec):
    """Transform the QC input specifications `jobrec` into the command
    and files for MP2D: jobrec@i -> mp2drec@i.

    Parameters
    ----------
    jobrec : dict
        Nested dictionary with input specifications for MP2D in generic
        QC terms.

    Returns
    -------
    mp2drec : dict
        Nested dictionary with input specification for MP2D in
        program-specific commands and files.

    """
    # temp until actual options object
    jobrec['extras']['info'] = dashparam.from_arrays(
        name_hint=jobrec['model']['method'],
        level_hint=jobrec['keywords'].get('level_hint', None),
        param_tweaks=jobrec['keywords'].get('params_tweaks', None),
        dashcoeff_supplement=jobrec['keywords'].get('dashcoeff_supplement', None))

    # this is what the mp2d program needs, not what the job needs
    # * form Params.txt string that governs dispersion calc
    # * form mp2d_geometry.psi4 string that supplies geometry to dispersion calc
    # * form command and arguments

    modulerec = {}
    modulerec['infiles'] = {}
    modulerec['infiles']['Params.txt'] = mp2d_coeff_formatter(jobrec['extras']['info']['dashlevel'],
                                                              jobrec['extras']['info']['dashparams'])

    # Need 'real' field later and that's only guaranteed for molrec
    molrec = qcel.molparse.from_schema(jobrec['molecule'])
    xyz = qcel.molparse.to_string(molrec, dtype='xyz', units='Angstrom', ghost_format='')
    modulerec['infiles']['mp2d_geometry.psi4'] = '\n'.join(['molecule {', '', *xyz.split('\n')[2:-1], 'units angstrom', '}'])

    modulerec['outfiles'] = [
        'mp2d_gradient',
    ]
    modulerec['env'] = {
            'HOME': os.environ.get('HOME'),
            'PATH': os.pathsep.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(os.pathsep) if x != '']) + \
                    os.pathsep + os.environ.get('PATH'),
            'MP2D_PARAMPATH': '.',
            'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH'),
        }
    modulerec['blocking_files'] = None

    jobrec['molecule']['real'] = molrec['real']

    command = ['mp2d', 'mp2d_geometry.psi4']
    if jobrec['driver'] == 'gradient':
        command.append('-grad')
    if jobrec['extras']['info']['dashlevel'] == 'atmgr':
        command.append('-abc')
    modulerec['command'] = command

    return modulerec


def mp2d_harvest(jobrec, modulerec):
    """Process raw results from read-only `mp2drec` into Datum
    fields in returned `jobrec`: jobrec@i, mp2drec@io -> jobrec@io.

    Parameters
    ----------
    jobrec : dict
        Nested dictionary with input specifications for MP2D in generic
        QC terms.
    mp2drec : dict
        Nested dictionary with input specification and output collection
        from MP2D in program-specific commands, files, & output capture.

    Returns
    -------
    jobrec : dict
        Nested dictionary with input specification and output collection
        from MP2D in generic QC terms.

    """
    # amalgamate output
    text = modulerec['stdout']
    text += '\n  <<<  MP2D Results  >>>\n'

    for fl, contents in modulerec['outfiles'].items():
        if contents is not None:
            text += f'\n  MP2D scratch file {fl} has been read.\n'
            text += contents

    # parse energy output (could go further and break into UCHF, CKS)
    real = np.array(jobrec['molecule']['real'])
    full_nat = real.shape[0]
    real_nat = np.sum(real)

    for ln in modulerec['stdout'].splitlines():
        if re.search('MP2D dispersion correction v', ln):
            version = ln.replace('MP2D dispersion correction', '').replace('-', '').strip().lower()
        elif re.match('   MP2D dispersion correction Eh', ln):
            ene = Decimal(ln.split()[4])
        elif re.match('normal termination of mp2d', ln):
            break
    else:
        if not ((real_nat == 1) and (jobrec['driver'] == 'gradient')):
            raise ValueError('Unsuccessful run. Possibly -D variant not available in dftd3 version.')

    # parse gradient output
    if modulerec['outfiles']['mp2d_gradient'] is not None:
        srealgrad = modulerec['outfiles']['mp2d_gradient']
        realgrad = np.fromstring(srealgrad, count=3 * real_nat, sep=' ').reshape((-1, 3))

    if jobrec['driver'] == 'gradient':
        ireal = np.argwhere(real).reshape((-1))
        fullgrad = np.zeros((full_nat, 3))
        try:
            fullgrad[ireal, :] = realgrad
        except NameError as exc:
            raise NameError('Unsuccessful gradient collection.') from exc

    qcvkey = jobrec['extras']['info']['fctldash'].upper()

    calcinfo = []
    calcinfo.append(qcel.Datum('DISPERSION CORRECTION ENERGY', 'Eh', ene))
    calcinfo.append(qcel.Datum('2-BODY DISPERSION CORRECTION ENERGY', 'Eh', ene))
    if qcvkey:
        calcinfo.append(qcel.Datum(f'{qcvkey} DISPERSION CORRECTION ENERGY', 'Eh', ene))

    if jobrec['driver'] == 'gradient':
        calcinfo.append(qcel.Datum('DISPERSION CORRECTION GRADIENT', 'Eh/a0', fullgrad))
        calcinfo.append(qcel.Datum('2-BODY DISPERSION CORRECTION GRADIENT', 'Eh/a0', fullgrad))
        if qcvkey:
            calcinfo.append(qcel.Datum(f'{qcvkey} DISPERSION CORRECTION GRADIENT', 'Eh/a0', fullgrad))

    calcinfo1 = {info.label: info for info in calcinfo}
    text += qcel.datum.print_variables(calcinfo1)
    calcinfo = {info.label: info.data for info in calcinfo}
    calcinfo = qcel.util.unnp(calcinfo, flat=True)

    jobrec['stdout'] = text
    jobrec['extras']['qcvars'] = calcinfo

    prov = {}
    prov['creator'] = 'mp2d'
    prov['routine'] = sys._getframe().f_code.co_name
    prov['version'] = version
    jobrec['provenance'] = prov

    return jobrec


def mp2d_coeff_formatter(dashlvl, dashcoeff):
    """Return strings for MP2D program parameter file.

    Parameters
    ----------
    dashlvl : {'dmp2'}
        Level of dispersion correction.
    dashcoeff : dict
        Dictionary fully specifying non-fixed parameters for `dashlvl` to drive MP2D.

    Returns
    -------
    str
        Suitable for `Params.txt` file.

    """
    dashformatter = """a_one: {:.6f} a_two: {:.6f} rcut: {:.6f} width: {:.6f} s_8: {:.6f}\n"""

    dashlvl = dashlvl.lower()
    if dashlvl == 'dmp2':
        return dashformatter.format(dashcoeff['a1'], dashcoeff['a2'], dashcoeff['rcut'], dashcoeff['w'], dashcoeff['s8'])
    else:
        raise ValueError(
            """-D correction level %s is not available. Choose among %s.""" % (dashlvl, dashcoeff.keys()))
