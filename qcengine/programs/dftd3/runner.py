"""Compute dispersion correction using Grimme's DFTD3 executable."""

import copy
import json
import os
import pathlib
import pprint
import re
import socket
import sys
from decimal import Decimal

import numpy as np
import qcelemental as qcel
from qcelemental.models import FailedOperation, Result

#from ..pdict import PreservingDict
from . import dashparam
from ...util import execute, which
from ..executor import ProgramExecutor
from ...extras import parse_dertype, provenance_stamp

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class DFTD3Executor(ProgramExecutor):

    _defaults = {
        "name": "DFTD3",
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

    def found(self) -> bool:
        return which('dftd3', return_bool=True)

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':

        if not which('dftd3', return_bool=True):
            raise ImportError("Could not find dftd3 in the envvar path.")

        # Setup the job
        input_data = input_data.copy().dict()
        input_data["success"] = False

        output_data = run_json(input_data)

        if output_data["success"]:
            return Result(**output_data)
        return FailedOperation(
            success=output_data.pop("success", False), error=output_data.pop("error"), input_data=output_data)


def run_json(jobrec):
    """
    An implementation of the QC JSON Schema (molssi-qc-schema.readthedocs.io/en/latest/index.html#) implementation in Psi4.

    Parameters
    ----------
    jobrec : JSON
        Please see molssi-qc-schema.readthedocs.io/en/latest/spec_components.html for further details.

    """
    # pp.pprint(jobrec)

    # This is currently a forced override
    if jobrec["schema_name"] not in ["qc_schema_input", "qcschema_input"]:
        raise KeyError(f"""Schema name of '{jobrec["schema_name"]}' not understood.""")

    if jobrec["schema_version"] != 1:
        raise KeyError(f"""Schema version of '{jobrec["schema_version"]}' not understood.""")

    jobrec['provenance'] = provenance_stamp(sys._getframe().f_code.co_name + '.' + __name__)

    # strip engine hint
    mtd = jobrec['model']['method']
    if mtd.startswith('d3-'):
        jobrec['model']['method'] = mtd[3:]

    # jobrec['model'] = {
    #     'method': name,
    #     'basis': '(auto)',
    # }
    # _, jobrec['driver'] = parse_dertype(kwargs['ptype'], max_derivative=1)

    # jobrec['options'] = opts
    # jobrec['options'] = copy.deepcopy(options)
    # # Set options
    # for k, v in json_data["keywords"].items():
    #     core.set_global_option(k, v)

    kw = jobrec["keywords"].copy()
    try:
        dftd3_driver(jobrec)
    except Exception as exc:
        import traceback
        jobrec['success'] = False
        jobrec['error'] = {
            'error_type': type(exc).__name__,
            'error_message': ''.join(traceback.format_exception(*sys.exc_info())),
        }
        raise exc

    jobrec['success'] = True
    #    for k, v in jobrec["extras"]["qcvars"].items():
    #        v = v.data
    #        if isinstance(v, np.ndarray):
    #            v = v.ravel().tolist()
    #        elif isinstance(v, Decimal):
    #            v = float(v)
    #
    #        jobrec["extras"]["qcvars"][k] = v

    jobrec["extras"]["qcvars"]["CURRENT ENERGY"] = jobrec["extras"]['qcvars']['DISPERSION CORRECTION ENERGY']
    jobrec['properties'] = {"return_energy": jobrec["extras"]['qcvars']['CURRENT ENERGY']}

    if jobrec['driver'] == 'energy':
        jobrec["return_result"] = jobrec["properties"]["return_energy"]
    elif jobrec['driver'] == 'gradient':
        jobrec["extras"]['qcvars']['CURRENT GRADIENT'] = copy.deepcopy(
            jobrec["extras"]['qcvars']['DISPERSION CORRECTION GRADIENT'])
        jobrec["return_result"] = jobrec["extras"]["qcvars"]["CURRENT GRADIENT"]

    jobrec["molecule"]["real"] = list(jobrec["molecule"]["real"])
    #    jobrec["extras"] = {"qcvars": jobrec.pop("qcvars"),
    #                       "info": jobrec.pop("keywords")}
    jobrec["keywords"] = kw

    return jobrec


def run_dftd3(name, molecule, options, **kwargs):
    """QCDriver signature for computing `name` on `molecule` with `options` with engine `DFTD3`."""

    # * ONLY takes self-sufficient fctl-dash for name
    # * tweakparams are only valid options

    opts = {}

    jobrec = {}
    jobrec['schema_name'] = 'qcschema_input'
    jobrec['schema_version'] = 1  # --> 2
    jobrec['provenance'] = provenance_stamp(sys._getframe().f_code.co_name + '.' + __name__)

    # strip engine hint
    if name.startswith('d3-'):
        name = name[3:]

    jobrec.update(molecule.to_schema(dtype=1))  # --> jobrec['molecule'] =
    jobrec['model'] = {
        'method': name,
        'basis': '(auto)',
    }
    _, jobrec['driver'] = parse_dertype(kwargs['ptype'], max_derivative=1)
    jobrec['keywords'] = opts
    #jobrec['options'] = copy.deepcopy(options)

    try:
        dftd3_driver(jobrec)
    except Exception as exc:
        import traceback
        jobrec['success'] = False
        jobrec['error'] = {
            'error_type': type(exc).__name__,
            'error_message': ''.join(traceback.format_exception(*sys.exc_info())),
        }
        raise exc

    jobrec['success'] = True
    jobrec['extras']['qcvars']['CURRENT ENERGY'] = copy.deepcopy(
        jobrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'])
    if jobrec['driver'] == 'gradient':
        jobrec['extras']['qcvars']['CURRENT GRADIENT'] = copy.deepcopy(
            jobrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT'])

    return jobrec


def run_dftd3_from_arrays(molrec,
                          name_hint=None,
                          level_hint=None,
                          param_tweaks=None,
                          ptype='energy',
                          dashcoeff_supplement=None,
                          verbose=1):
    """Specialized signature disentangling dispersion level and
    parameters for computing on `molecule` with engine `DFTD3`. See
    `dashparam.from_array` for parameter details.

    """
    jobrec = {}
    jobrec['schema_name'] = 'qc_schema_input'  # --> qcschema_input
    jobrec['schema_version'] = 1  # --> 2
    jobrec['provenance'] = provenance_stamp(sys._getframe().f_code.co_name + '.' + __name__)

    # strip engine hint
    if name_hint.startswith('d3-'):
        name_hint = name_hint[3:]

    opts = {}
    opts['level_hint'] = level_hint
    opts['params_tweaks'] = param_tweaks
    opts['dashcoeff_supplement'] = dashcoeff_supplement

    jobrec.update(qcel.molparse.to_schema(molrec, dtype=1))  # --> jobrec['molecule'] =
    jobrec['model'] = {
        'method': name_hint,
        'basis': '(auto)',
    }
    _, jobrec['driver'] = parse_dertype(ptype, max_derivative=1)
    jobrec['keywords'] = opts
    #jobrec['options'] = copy.deepcopy(options)

    try:
        dftd3_driver(jobrec)
    except Exception as exc:
        import traceback
        jobrec['success'] = False
        jobrec['error'] = {
            'error_type': type(exc).__name__,
            'error_message': ''.join(traceback.format_exception(*sys.exc_info())),
        }
        raise exc

    jobrec['success'] = True
    jobrec['extras']['qcvars']['CURRENT ENERGY'] = copy.deepcopy(
        jobrec['extras']['qcvars']['DISPERSION CORRECTION ENERGY'])
    if jobrec['driver'] == 'gradient':
        jobrec['extras']['qcvars']['CURRENT GRADIENT'] = copy.deepcopy(
            jobrec['extras']['qcvars']['DISPERSION CORRECTION GRADIENT'])

    return jobrec


def dftd3_driver(jobrec, verbose=1):
    """Drive the jobrec@i (input) -> dftd3rec@i -> dftd3rec@io -> jobrec@io (returned) process."""

    return module_driver(
        jobrec=jobrec, module_label='dftd3', plant=dftd3_plant, harvest=dftd3_harvest, verbose=verbose)


def module_driver(jobrec, module_label, plant, harvest, verbose=1):
    """Drive the jobrec@i (input) -> modulerec@i -> modulerec@io -> jobrec@io (returned) process.

    Input Fields
    ------------

    Optional Input Fields
    ---------------------

    Output Fields
    -------------

    Optional Output Fields
    ----------------------

    """
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

    #run(modulerec)  # updates modulerec
    #execute(modulerec)  # updates modulerec

    command = modulerec.pop('command')
    infiles = modulerec.pop('infiles')
    outfiles = modulerec.pop('outfiles')
    env = modulerec.pop('env')
    blocking_files = modulerec.pop('blocking_files')

    ans, dans = execute(command, infiles, outfiles, **{
        'scratch_messy': True,
        'environment': env,
        'blocking_files': blocking_files
    })
    modulerec.update(dans)
    modulerec.update({'command': command, 'infiles': infiles, 'env': env})

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


def dftd3_plant(jobrec):
    """Transform the QC input specifications `jobrec` into the command
    and files for DFTD3: jobrec@i -> dftd3rec@i.

    Parameters
    ----------
    jobrec : dict
        Nested dictionary with input specifications for DFTD3 in generic
        QC terms.

    Returns
    -------
    dftd3rec : dict
        Nested dictionary with input specification for DFTD3 in
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
    dftd3rec = dashparam.from_arrays(
        name_hint=jobrec['model']['method'],
        level_hint=jobrec['keywords'].get('level_hint', None),
        param_tweaks=jobrec['keywords'].get('params_tweaks', None),
        dashcoeff_supplement=jobrec['keywords'].get('dashcoeff_supplement', None))
    # sketchy: adding to options during planting season
    jobrec['keywords'].update(dftd3rec)

    # this is what the dftd3 program needs, not what the job needs
    # * form dftd3_parameters string that governs dispersion calc
    # * form dftd3_geometry string that supplies geometry to dispersion calc
    # * form command and arguments

    dftd3rec['infiles'] = {}
    dftd3rec['infiles']['.dftd3par.local'] = dftd3_coeff_formatter(dftd3rec['dashlevel'], dftd3rec['dashparams'])

    # Have to pass outer level, not jobrec['molecule'] b/c qc_schema is in outer
    # Need 'real' field later and that's only guaranteed for molrec
    molrec = qcel.molparse.from_schema(jobrec)
    dftd3rec['infiles']['dftd3_geometry.xyz'] = qcel.molparse.to_string(
        molrec, dtype='xyz', units='Angstrom', ghost_format='')

    dftd3rec['outfiles'] = [
        'dftd3_gradient',
        'dftd3_abc_gradient',
    ]
    dftd3rec['env'] = {
            'HOME': os.environ.get('HOME'),
            'PATH': os.pathsep.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(os.pathsep) if x != '']) + \
                    os.pathsep + os.environ.get('PATH'),
            'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH'),
            'NONSENSE': None
        }
    dftd3rec['blocking_files'] = [os.path.join(pathlib.Path.home(), '.dftd3par.' + socket.gethostname())]

    jobrec['molecule']['real'] = molrec['real']

    command = ['dftd3', 'dftd3_geometry.xyz']
    if jobrec['driver'] == 'gradient':
        command.append('-grad')
    if dftd3rec['dashlevel'] == 'atmgr':
        command.append('-abc')
    dftd3rec['command'] = command

    return dftd3rec


def dftd3_harvest(jobrec, dftd3rec):
    """Process raw results from read-only `dftd3rec` into Datum
    fields in returned `jobrec`: jobrec@i, dftd3rec@io -> jobrec@io.

    Parameters
    ----------
    jobrec : dict
        Nested dictionary with input specifications for DFTD3 in generic
        QC terms.
    dftd3rec : dict
        Nested dictionary with input specification and output collection
        from DFTD3 in program-specific commands, files, & output capture.

    Returns
    -------
    jobrec : dict
        Nested dictionary with input specification and output collection
        from DFTD3 in generic QC terms.

    Notes
    -----
    Central to harvesting is the fact (to the planting, not to the DFTD3
    program) that 2-body and 3-body are run separately. Because of how
    damping functions work (see GH:psi4/psi4#1407), some 2-body damping
    schemes can give wrong answers for 3-body. And because 3-body is
    set to run with some dummy values, the 2-body values are no good.

    """
    try:
        jobrec['molecule']['real']
        jobrec['driver']
        jobrec['provenance']
        jobrec['keywords']['fctldash']
    except KeyError as exc:
        raise KeyError('Required field ({}) missing among ({})'.format(str(exc), list(jobrec.keys()))) from exc

    try:
        dftd3rec['stdout']
    except KeyError as exc:
        raise KeyError('Required field ({}) missing among ({})'.format(str(exc), list(dftd3rec.keys()))) from exc

    # amalgamate output
    text = dftd3rec['stdout']
    text += '\n  <<<  DFTD3 Results  >>>\n'

    for fl, contents in dftd3rec['outfiles'].items():
        if contents is not None:
            text += f'\n  DFTD3 scratch file {fl} has been read.\n'
            text += contents

    # parse energy output (could go further and break into E6, E8, E10 and Cn coeff)
    real = np.array(jobrec['molecule']['real'])
    full_nat = real.shape[0]
    real_nat = np.sum(real)

    for ln in dftd3rec['stdout'].splitlines():
        if re.search('DFTD3 V', ln):
            version = ln.replace('DFTD3', '').replace('|', '').strip().lower()
        elif re.match(' Edisp /kcal,au', ln):
            ene = Decimal(ln.split()[3])
        elif re.match(r" E6\(ABC\) \"   :", ln):  # c. v3.2.0
            raise ValueError("Cannot process ATM results from DFTD3 prior to v3.2.1.")
        elif re.match(r""" E6\(ABC\) /kcal,au:""", ln):
            atm = Decimal(ln.split()[-1])
        elif re.match(' normal termination of dftd3', ln):
            break
    else:
        if not ((real_nat == 1) and (jobrec['driver'] == 'gradient')):
            raise ValueError('Unsuccessful run. Possibly -D variant not available in dftd3 version.')

    # parse gradient output
    # * DFTD3 crashes on one-atom gradients. Avoid the error (above) and just force the correct result (below).
    if dftd3rec['outfiles']['dftd3_gradient'] is not None:
        srealgrad = dftd3rec['outfiles']['dftd3_gradient'].replace('D', 'E')
        realgrad = np.fromstring(srealgrad, count=3 * real_nat, sep=' ').reshape((-1, 3))
    elif real_nat == 1:
        realgrad = np.zeros((1, 3))

    if dftd3rec['outfiles']['dftd3_abc_gradient'] is not None:
        srealgrad = dftd3rec['outfiles']['dftd3_abc_gradient'].replace('D', 'E')
        realgradabc = np.fromstring(srealgrad, count=3 * real_nat, sep=' ').reshape((-1, 3))
    elif real_nat == 1:
        realgradabc = np.zeros((1, 3))

    if jobrec['driver'] == 'gradient':
        ireal = np.argwhere(real).reshape((-1))
        fullgrad = np.zeros((full_nat, 3))
        rg = realgradabc if (dftd3rec['dashlevel'] == 'atmgr') else realgrad
        try:
            fullgrad[ireal, :] = rg
        except NameError as exc:
            raise Dftd3Error('Unsuccessful gradient collection.') from exc

    qcvkey = jobrec['keywords']['fctldash'].upper()

    # OLD WAY
    calcinfo = []
    if dftd3rec['dashlevel'] == 'atmgr':
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
    jobrec['extras'] = {}
    jobrec['extras']['qcvars'] = calcinfo

    prov = {}
    prov['creator'] = 'dftd3'
    prov['routine'] = sys._getframe().f_code.co_name
    prov['version'] = version
    jobrec['provenance'] = prov

    return jobrec


def dftd3_coeff_formatter(dashlvl, dashcoeff):
    """Return strings for DFTD3 program parameter file.

             s6      rs6      s18     rs8     alpha6      version
             ------------------------------------------------------
    d2:      s6      sr6      s8=0.0  a2=None alpha6      version=2
    d3zero:  s6      sr6      s8      a2=sr8  alpha6      version=3
    d3bj:    s6      a1       s8      a2      alpha6=None version=4
    d3mzero: s6      sr6      s8      beta    alpha6=14.0 version=5
    d3mbj:   s6      a1       s8      a2      alpha6=None version=6
    atmgr:   s6=1.0  sr6=None s8=None a2=None alpha6      version=3 (needs -abc, too)

    Parameters
    ----------
    dashlvl : {'d2', 'd3zero', d3bj', 'd3mzero', 'd3mbj', 'atmgr'}
        Level of dispersion correction.
    dashcoeff : dict
        Dictionary fully specifying non-fixed parameters (table above) for `dashlvl` to drive DFTD3.

    Notes
    -----
    The `atmgr` dashlvl is intended for use only to get the three-body Axilrod-Teller-Muto
    three body dispersion correction. Therefore, dummy parameters are passed for two-body damping
    function, and it will give garbage for two-body component of dispersion correction.

    Returns
    -------
    str
        Suitable for `.dftd3par` file.

    """
    dashformatter = """{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:6}\n"""

    dashlvl = dashlvl.lower()
    if dashlvl == 'd2':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['sr6'], 0.0, 0.0, dashcoeff['alpha6'], 2)
    elif dashlvl == 'd3zero':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['sr6'], dashcoeff['s8'], dashcoeff['sr8'],
                                    dashcoeff['alpha6'], 3)
    elif dashlvl == 'd3bj':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['a1'], dashcoeff['s8'], dashcoeff['a2'], 0.0, 4)
    elif dashlvl == 'd3mzero':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['sr6'], dashcoeff['s8'], dashcoeff['beta'], 14.0, 5)
    elif dashlvl == 'd3mbj':
        return dashformatter.format(dashcoeff['s6'], dashcoeff['a1'], dashcoeff['s8'], dashcoeff['a2'], 0.0, 6)
    elif dashlvl == 'atmgr':
        # need to set first four parameters to something other than None, otherwise DFTD3 gets mad or a bit wrong
        return dashformatter.format(1.0, 0.0, 0.0, 0.0, dashcoeff['alpha6'], 3)
    else:
        raise ValueError(f"""-D correction level {dashlvl} is not available. Choose among {dashcoeff.keys()}.""")


"""
Notes
-----
The DFTD3 executable must be independently compiled and found in :envvar:`PATH` or :envvar:`PSIPATH`.
research site: https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3
Psi4 mode: When `psi4` the python module is importable at `import qcdb`
           time, Psi4 mode is activated, with the following alterations:
           * output goes to output file
           * gradient returned as psi4.core.Matrix, not list o'lists
           * scratch is written to randomly named subdirectory of psi scratch
           * psivar "DISPERSION CORRECTION ENERGY" is set
           * `verbose` triggered when PRINT keywork of SCF module >=3
"""
