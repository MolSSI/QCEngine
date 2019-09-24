"""
Calls the NWChem executable.
"""
import copy
from decimal import Decimal
from pathlib import Path
import pprint
import re
from subprocess import Popen, PIPE
from typing import Any, Dict, Optional, Tuple

import numpy as np

import qcelemental as qcel
from qcelemental.models import Provenance, Result
from qcelemental.util import safe_version, which

from ...exceptions import InputError
from ...util import execute, temporary_directory
from ..model import ProgramHarness

class TurbomoleHarness(ProgramHarness):
    
    _defaults = {
        "name": "Turbomole",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": False,
        "node_parallel": True,
        "managed_memory": True,
    }

    version_cache: Dict[str, str] = {}
    
    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which('define',
                     return_bool=True,
                     raise_error=raise_error,
                     raise_msg='Please install via http://www.cosmologic.de/turbomole/home.html')

    def get_version(self) -> str:
        which_prog = which('define')
        if which_prog not in self.version_cache:
            # We use basically a dummy stdin as we dont want to pipe any real
            # input into define. We only want to parse the version number from
            # the string.
            stdout = self.execute_define("\n")
            version_re  = re.compile("TURBOMOLE V([\d\.]+)")
            mobj = version_re.search(stdout)
            version = mobj[1]
            self.version_cache[which_prog] = safe_version(version)
        return self.version_cache[which_prog]

    def compute(self, input_model: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs NWChem in executable mode
        """
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        if 'There is an error in the input file' in dexe["stdout"]:
            raise InputError(dexe["stdout"])

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_model)

    def prepare_define_stdin(self, method: str, basis: str, molecule: 'Molecule') -> str:
        charge = molecule.molecular_charge
        mult = molecule.molecular_multiplicity
        unrestricted = False

        def occ_num_mo_data(charge: int, mult: int,
                            unrestricted: Optional[bool] = False) -> str:
            """Handles the 'Occupation Number & Molecular Orbital' section
            of define."""
            unrestricted = unrestricted or (mult != 1)
            unpaired = mult - 1
            charge = int(charge)

            occ_num_mo_data_stdin = f"""eht
            y
            {charge}
            y
            """
            if unrestricted:
                # Somehow Turbomole/define asks us if we want to write
                # natural orbitals... we don't want to.
                occ_num_mo_data_stdin = f"""eht
                y
                {charge}
                n
                u {unpaired}
                *
                n
                """
            return occ_num_mo_data_stdin

        kwargs = {
            "init_guess": occ_num_mo_data(charge, mult, unrestricted),
            "title": "QCEngine Turbomole",
            "xcfunc": "pbe0",
            "grid": "m4",
            "scf_conv": 8,
            "scf_iters": 150,
        }

        stdin = """
        {title}
        a coord
        *
        no
        b
        all def2-SVP
        *
        {init_guess}
        dft
        on
        func
        {xcfunc}
        grid
        {grid}

        scf
        conv
        {scf_conv}
        iter
        {scf_iters}

        *
        """.format(**kwargs)

        return stdin

    def build_input(self, input_model: 'ResultInput', config: 'JobConfig',
                    template: Optional[str] = None) -> Dict[str, Any]:
        turbomolrec = {
            'infiles': {},
            'scratch_directory': config.scratch_directory,
        }

        # Handle molecule
        # TODO: what's up with moldata? Do I need it?
        coord_str, moldata = input_model.molecule.to_string(dtype='turbomole', return_data=True)

        # Prepare stdin for define call
        model = input_model.model
        stdin = self.prepare_define_stdin(model.method, model.basis,
                                          input_model.molecule,
        )
        # TODO: I would probably like to use temporary_directory but the
        # context manager seems to return a str instead of a Path?!
        with temporary_directory(suffix="_define_scratch") as tmpdir:
            tmpdir = Path(tmpdir)
            with open(tmpdir / "coord", "w") as handle:
                handle.write(coord_str)
            stdout = self.execute_define(stdin, cwd=tmpdir)
            # The define scratch will be populated by some files that we want to keep
            to_keep = "basis auxbasis coord control alpha beta mos".split()

            for fn in to_keep:
                full_fn = tmpdir / fn
                if not full_fn.exists():
                    continue
                with open(full_fn) as handle:
                    turbomolrec['infiles'][fn] = handle.read()

        turbomolrec['command'] = ["dscf"]

        return turbomolrec

    def execute_define(self, stdin: str, cwd: Optional["Path"] = None) -> str:
        # TODO: replace this with a call to the default execute provided by QCEngine
        proc = Popen("define", universal_newlines=True,
                     stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=cwd,
        )
        # TODO: add timeout? Unless the disk hangs this should never take long...
        stdout, _ = proc.communicate(stdin)#, timeout=3)
        proc.terminate()
        return stdout

    def execute(self,
                inputs: Dict[str, Any],
                *,
                extra_outfiles=None,
                extra_commands=None,
                scratch_name=None,
                timeout=None) -> Tuple[bool, Dict]:

        import pdb; pdb.set_trace()
        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            # TODO: scratch_messy?
            # scratch_messy=False,
        )
        return success, dexe

    def parse_output(self,
                     outfiles: Dict[str, str],
                     input_model: 'ResultInput') -> 'Result':  # lgtm: [py/similar-function]

        stdout = outfiles.pop("stdout")

        # nwmol, if it exists, is dinky, just a clue to geometry of nwchem results
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
            'provenance': Provenance(creator="Turbomole", version=self.get_version(),
                                     routine="turbomole"),
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
