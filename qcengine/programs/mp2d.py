"""Compute dispersion correction using Greenwell & Beran's MP2D executable."""

import copy
import json
import os
import pathlib
import pprint
import re
import socket
import sys
import traceback
from decimal import Decimal
from typing import Any, Dict, Optional

import numpy as np
import qcelemental as qcel
from qcelemental.models import FailedOperation, Result

from .dftd3 import dashparam
from ..util import execute, which
from .executor import ProgramExecutor
from ..extras import parse_dertype, provenance_stamp

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
    def found() -> bool:
        return which('mp2d', return_bool=True)

    def get_version(self) -> str:
        if not self.found():
            raise ImportError("Could not find MP2D in the shell path.")
        # Note: anything below v3.2.1 will return the help menu here. but that's fine as version compare evals to False.
        command = [which('mp2d'), '-version']
        import subprocess
        proc = subprocess.run(command, stdout=subprocess.PIPE)
        candidate_version = proc.stdout.decode('utf-8').strip()

        from pkg_resources import safe_version
        return safe_version(candidate_version)

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':

        if not which('mp2d', return_bool=True):
            raise ImportError("Could not find mp2d in the envvar path.")

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

    if jobrec['driver'].derint() > 1:
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


def module_driver(jobrec, module_label, plant, harvest, verbose=4):
    """Drive the jobrec@i (input) -> modulerec@i -> modulerec@io -> jobrec@io (returned) process."""

    if verbose > 2:
        print(f'[1] {module_label.upper()} JOBREC PRE-PLANT (j@i) <<<')
        pp.pprint(jobrec)
        print('>>>')

    modulerec = plant(jobrec)

    # test json roundtrip
    jmodulerec = json.dumps(modulerec)
    modulerec = json.loads(jmodulerec)

    if verbose > 3:
        print(f'[2] {module_label.upper()}REC PRE-SUBPROCESS (m@i) <<<')
        pp.pprint(modulerec)
        print('>>>\n')

    rc, dexe = execute(modulerec.get('command'),
                       modulerec.get('infiles'),
                       modulerec.pop('outfiles'),
                       **{
        'scratch_messy': True,
        'environment': modulerec.get('env'),
        'blocking_files': modulerec.get('blocking_files')
    })
    modulerec.update(dexe)  # updates modulerec

    if verbose > 3:
        print(f'[3] {module_label.upper()}REC POST-SUBPROCESS (m@io) <<<')
        pp.pprint(modulerec)
        print('>>>\n')

    harvest(jobrec, modulerec)  # updates jobrec

    if verbose > 1:
        print(f'[4] {module_label.upper()} JOBREC POST-HARVEST (j@io) <<<')
        pp.pprint(jobrec)
        print('>>>')

    return jobrec


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
    try:
        jobrec['driver']
        jobrec['model']['method']
        jobrec['keywords']
        jobrec['molecule']
    except KeyError as exc:
        raise KeyError('Required field ({}) missing among ({})'.format(str(exc), list(jobrec.keys()))) from exc

    # temp until actual options object
    jobrec['extras']['info'] = dashparam.from_arrays(
        name_hint=jobrec['model']['method'],
        level_hint=jobrec['keywords'].get('level_hint', None),
        param_tweaks=jobrec['keywords'].get('params_tweaks', None),
        dashcoeff_supplement=jobrec['keywords'].get('dashcoeff_supplement', None))

    # this is what the mp2d program needs, not what the job needs
    # * form mp2d_parameters string that governs dispersion calc
    # * form mp2d_geometry string that supplies geometry to dispersion calc
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
            'NONSENSE': None
        }
    modulerec['blocking_files'] = []

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
    try:
        jobrec['molecule']['real']
        jobrec['driver']
        jobrec['provenance']
        jobrec['extras']['info']['fctldash']
    except KeyError as exc:
        raise KeyError('Required field ({}) missing among ({})'.format(str(exc), list(jobrec.keys()))) from exc

    try:
        modulerec['stdout']
    except KeyError as exc:
        raise KeyError('Required field ({}) missing among ({})'.format(str(exc), list(modulerec.keys()))) from exc

    # amalgamate output
    text = modulerec['stdout']
    text += '\n  <<<  MP2D Results  >>>\n'

    for fl, contents in modulerec['outfiles'].items():
        if contents is not None:
            text += f'\n  MP2D scratch file {fl} has been read.\n'
            text += contents

    # parse energy output (could go further and break into E6, E8, E10 and Cn coeff)
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
    # * DFTD3 crashes on one-atom gradients. Avoid the error (above) and just force the correct result (below).
    if modulerec['outfiles']['mp2d_gradient'] is not None:
        srealgrad = modulerec['outfiles']['mp2d_gradient']
        realgrad = np.fromstring(srealgrad, count=3 * real_nat, sep=' ').reshape((-1, 3))

    if jobrec['driver'] == 'gradient':
        ireal = np.argwhere(real).reshape((-1))
        fullgrad = np.zeros((full_nat, 3))
        try:
            fullgrad[ireal, :] = realgrad
        except NameError as exc:
            raise Dftd3Error('Unsuccessful gradient collection.') from exc

    qcvkey = jobrec['extras']['info']['fctldash'].upper()

    # OLD WAY
    calcinfo = []
    if jobrec['extras']['info']['dashlevel'] == 'atmgr':
        calcinfo.append(qcel.Datum('DISPERSION CORRECTION ENERGY', 'Eh', atm))
        calcinfo.append(qcel.Datum('3-BODY DISPERSION CORRECTION ENERGY', 'Eh', atm))
        calcinfo.append(qcel.Datum('AXILROD-TELLER-MUTO 3-BODY DISPERSION CORRECTION ENERGY', 'Eh', atm))

        if jobrec['driver'] == 'gradient':
            calcinfo.append(qcel.Datum('DISPERSION CORRECTION GRADIENT', 'Eh/a0', fullgrad))
            calcinfo.append(qcel.Datum('3-BODY DISPERSION CORRECTION GRADIENT', 'Eh/a0', fullgrad))
            calcinfo.append(qcel.Datum('AXILROD-TELLER-MUTO 3-BODY DISPERSION CORRECTION GRADIENT', 'Eh/a0', fullgrad))

    else:
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

    # NEW WAY
    #module_vars = {}
    #module_vars['DISPERSION CORRECTION ENERGY'] =  ene
    #module_vars['{} DISPERSION CORRECTION ENERGY'.format(qcvkey)] =  ene
    #if jobrec['driver'] == 'gradient':
    #    module_vars['DISPERSION CORRECTION GRADIENT'] = fullgrad
    #    module_vars['{} DISPERSION CORRECTION GRADIENT'.format(qcvkey)] = fullgrad
    #
    #module_vars = PreservingDict(module_vars)
    #qcvars.build_out(module_vars)
    #calcinfo = qcvars.certify(module_vars)
    #text += print_variables(calcinfo)

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
