"""
Calls the Turbomole executable.
"""
import os
import re
from decimal import Decimal
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from qcelemental.models import AtomicResult, Provenance
from qcelemental.util import safe_version, which

from ...util import execute, temporary_directory
from ..model import ProgramHarness
from .define import execute_define, prepare_stdin
from .harvester import harvest
from .methods import KEYWORDS, METHODS


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
        return which(
            "define",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via http://www.cosmologic.de/turbomole/home.html",
        )

    def get_version(self) -> str:
        which_prog = which("define")
        if which_prog not in self.version_cache:
            # We use basically a dummy stdin as we dont want to pipe any real
            # input into define. We only want to parse the version number from
            # the string.
            with temporary_directory(suffix="_define_scratch") as tmpdir:
                tmpdir = Path(tmpdir)
                stdout = execute_define("\n", cwd=tmpdir)
            # Tested with V7.3 and V7.4.0
            version_re = re.compile("TURBOMOLE (?:rev\. )?(V.+?)\s+")
            mobj = version_re.search(stdout)
            version = mobj[1]
            self.version_cache[which_prog] = safe_version(version)
        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        # TODO: handle input errors?! But then define probably already crashed...
        # if 'There is an error in the input file' in dexe["stdout"]:
        # raise InputError(dexe["stdout"])

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_model)

    def build_input(
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        turbomolrec = {"infiles": {}, "outfiles": {"control": "control"}, "scratch_directory": config.scratch_directory}

        # Handle molecule
        # TODO: what's up with moldata? Do I need it?
        coord_str, moldata = input_model.molecule.to_string(dtype="turbomole", return_data=True)

        # Prepare stdin for define call
        model = input_model.model
        # geeopt will hold the for which to calculate the gradient.
        # 'x' corresponds to the ground state, 'a 1' would be the GS too.
        # 'a1 2' would be the 1st excited state of the irreducible group A1.
        # Right now only GS are supported, so this is hardcoded as 'x'.
        geoopt = "x" if input_model.driver.derivative_int() > 0 else ""
        stdin, subs = prepare_stdin(
            model.method,
            model.basis,
            input_model.keywords,
            input_model.molecule.molecular_charge,
            input_model.molecule.molecular_multiplicity,
            geoopt,
        )
        with temporary_directory(suffix="_define_scratch") as tmpdir:
            tmpdir = Path(tmpdir)
            with open(tmpdir / "coord", "w") as handle:
                handle.write(coord_str)
            stdout = execute_define(stdin, cwd=tmpdir)
            # The define scratch will be populated by some files that we want to keep
            to_keep = "basis auxbasis coord control alpha beta mos".split()

            for fn in to_keep:
                full_fn = tmpdir / fn
                if not full_fn.exists():
                    continue
                with open(full_fn) as handle:
                    turbomolrec["infiles"][fn] = handle.read()

        env = os.environ.copy()
        env["PARA_ARCH"] = "SMP"
        env["PARNODES"] = str(config.ncores)
        env["SMPCPUS"] = str(config.ncores)
        # TODO: set memory

        turbomolrec["environment"] = env

        keywords = input_model.keywords
        ri_calculation = any([keywords.get(ri_kw, False) for ri_kw in KEYWORDS["ri"]])

        # Set appropriate commands. We always need a reference wavefunction
        # so the first command will be dscf or ridft.
        commands = ["ridft"] if ri_calculation else ["dscf"]
        # ricc2 will also calculate the gradient
        if model.method in METHODS["ricc2"]:
            commands.append("ricc2")
        # Gradient calculation for DFT/HF
        elif input_model.driver.derivative_int() == 1:
            grad_command = "rdgrad" if ri_calculation else "grad"
            commands.append(grad_command)
        elif input_model.driver.derivative_int() == 2:
            freq_command = "aoforce"
            commands.append(freq_command)
            # Add
            #   $noproj
            #   $nprhessian file=nprhessian
            # to control file.
            turbomolrec["outfiles"]["hessian"] = "nprhessian"

        if input_model.driver.derivative_int() == 1:
            turbomolrec["outfiles"]["gradient"] = "gradient"
        command = ["; ".join(commands)]
        turbomolrec["command"] = command

        # TODO: check if the chosen commands are available with which()?

        return turbomolrec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            inputs["outfiles"],
            shell=True,
            # TODO: scratch_messy?
            # scratch_messy=False,
        )
        return success, dexe

    def parse_output(
        self, outfiles: Dict[str, str], input_model: "AtomicInput"
    ) -> "AtomicResult":  # lgtm: [py/similar-function]

        stdout = outfiles.pop("stdout")

        # nwmol, if it exists, is dinky, just a clue to geometry of nwchem results
        qcvars, gradient, hessian = harvest(input_model.molecule, stdout, **outfiles)

        if gradient is not None:
            qcvars["CURRENT GRADIENT"] = gradient

        if hessian is not None:
            qcvars["CURRENT HESSIAN"] = hessian

        retres = qcvars[f"CURRENT {input_model.driver.upper()}"]
        if isinstance(retres, Decimal):
            retres = float(retres)

        output_data = input_model.dict()
        output_data["extras"]["outfiles"] = outfiles
        output_data["properties"] = {}
        output_data["provenance"] = Provenance(creator="Turbomole", version=self.get_version(), routine="turbomole")
        output_data["return_result"] = retres
        output_data["stdout"] = stdout
        output_data["success"] = True

        return AtomicResult(**output_data)
