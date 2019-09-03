"""
Calls the PylibEFP interface to LibEFP.
"""
import pprint
from typing import Dict

import qcelemental as qcel
from qcelemental.models import Provenance, Result
from qcelemental.util import safe_version, which_import

from ..exceptions import InputError  #, RandomError, ResourceError, UnknownError
from .model import ProgramHarness

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class PylibEFPHarness(ProgramHarness):

    _defaults = {
        "name": "PylibEFP",
        "scratch": False,
        "thread_safe": False,
        "thread_parallel": False,  # can be but not the way Psi usually builds it
        "node_parallel": False,
        "managed_memory": False,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which_import('pylibefp',
                            return_bool=True,
                            raise_error=raise_error,
                            raise_msg='Please install via `conda install pylibefp -c psi4`.')

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import('pylibefp')
        if which_prog not in self.version_cache:
            import pylibefp
            self.version_cache[which_prog] = safe_version(pylibefp.__version__)

        candidate_version = self.version_cache[which_prog]

        if "undef" in candidate_version:
            raise TypeError(
                "Using custom build without tags. Please pull git tags with `git pull origin master --tags`.")

        return candidate_version

    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs PylibEFP in API mode
        """
        self.found(raise_error=True)

        #        if parse_version(self.get_version()) < parse_version("1.2"):
        #            raise ResourceError("Psi4 version '{}' not understood.".format(self.get_version()))

        # Setup the job
        input_data = input_model.dict(encoding="json")
        #        input_data["nthreads"] = config.ncores
        #        input_data["memory"] = int(config.memory * 1024 * 1024 * 1024 * 0.95)  # Memory in bytes
        input_data["success"] = False
        input_data["return_output"] = True

        # initialize EFP fragments
        import pylibefp
        # fix units when parsing efp string
        efpobj = pylibefp.from_dict({**input_model.molecule.extras['efp_molecule']['extras'], 'units':'Bohr'})

        # print efp geom in [A]
        print(efpobj.banner())
        print(efpobj.geometry_summary(units_to_bohr=qcel.constants.bohr2angstroms))

        if input_model.model.method != 'efpefp':
            raise InputError(f"Method not efpefp: {input_model.model.method}")

        # set keywords
        # * make keywords keys case insensitive
        # * label changes the accepted names of the keywords (xr vs. exch)
        # * append changes the defaults upon which they act (off for libefp vs. on for psi)
        opts = {k.lower(): v for k, v in input_model.keywords.items()}
        keywords_label = opts.pop('keywords_label', 'libefp').lower()
        results_label = opts.pop('results_label', 'libefp').lower()
        try:
            efpobj.set_opts(opts, label=keywords_label, append=keywords_label)
        except pylibefp.EFPSyntaxError as e:
            raise InputError(e.message)

        if input_model.driver == 'energy':
            do_gradient = False
        elif input_model.driver == 'gradient':
            do_gradient = True
        else:
            raise InputError

        # compute and report
        efpobj.compute(do_gradient=do_gradient)
        print(efpobj.energy_summary(label=results_label))

        ene = efpobj.get_energy(label=results_label)

        pp.pprint(ene)
        print('<<< get_opts():  ', efpobj.get_opts(), '>>>')
        #print('<<< summary():   ', efpobj.summary(), '>>>')
        print('<<< get_energy():', ene, '>>>')
        print('<<< get_atoms(): ', efpobj.get_atoms(), '>>>')
        print(efpobj.energy_summary())
        print(efpobj.geometry_summary(units_to_bohr=qcel.constants.bohr2angstroms))
        print(efpobj.geometry_summary(units_to_bohr=1.0))

        ###### psi4 proc
        #def run_efp(name, **kwargs):
        #        try:
        #            efpobj = efp_molecule.EFP
        #        except AttributeError:
        #            raise ValidationError("""Method 'efp' not available without EFP fragments in molecule""")
        #
        #        core.set_variable('EFP ELST ENERGY', ene['electrostatic'] + ene['charge_penetration'] + ene['electrostatic_point_charges'])
        #        core.set_variable('EFP IND ENERGY', ene['polarization'])
        #        core.set_variable('EFP DISP ENERGY', ene['dispersion'])
        #        core.set_variable('EFP EXCH ENERGY', ene['exchange_repulsion'])
        #        core.set_variable('EFP TOTAL ENERGY', ene['total'])
        #        core.set_variable('CURRENT ENERGY', ene['total'])
        #
        #        if do_gradient:
        #            core.print_out(efpobj.gradient_summary())
        #
        #            core.set_variable('EFP TORQUE', torq)
        #
        #                output_data = input_data

        if input_model.driver == 'energy':
            retres = ene['total']


#        elif input_model.driver == 'gradient':
#            torq = efpobj.get_gradient()

        output_data = {
            'schema_name': 'qcschema_output',
            'schema_version': 1,
            'extras': {
                'local_properties': ene,
                #     'outfiles': outfiles,
            },
            'properties': {},
            'provenance': Provenance(creator="PylibEFP", version=self.get_version(), routine="pylibefp"),
            'return_result': retres,
            #'stdout': stdout,
        }

        output_data['success'] = True
        return Result(**{**input_model.dict(), **output_data})
