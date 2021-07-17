"""Compute quantum chemistry using Iowa State's GAMESS executable."""

import copy
import pprint
import sys
import traceback
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple

from qcelemental.models import AtomicInput, AtomicResult, BasisSet, Provenance
from qcelemental.util import safe_version, which

from ...exceptions import InputError, UnknownError
from ...util import execute
from ..model import ProgramHarness
from ..qcvar_identities_resources import build_atomicproperties, build_out
from .germinate import muster_modelchem
from .harvester import harvest
from .keywords import format_keywords

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)


class GAMESSHarness(ProgramHarness):
    """

    Notes
    -----
    Required edits to the ``rungms`` script are as follows::
        set SCR=./                                  # will be managed by QCEngine instead
        set USERSCR=./                              # ditto
        set GMSPATH=/home/psilocaluser/gits/gamess  # full path to installation

    """

    _defaults = {
        "name": "GAMESS",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": True,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "rungms",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via https://www.msg.chem.iastate.edu/GAMESS/GAMESS.html",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("rungms")
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "v.inp"], {"v.inp": ""})

            if success:
                for line in output["stdout"].splitlines():
                    if "GAMESS VERSION" in line:
                        branch = " ".join(line.strip(" *\t").split()[3:])
                self.version_cache[which_prog] = safe_version(branch)

        return self.version_cache[which_prog]

    def compute(self, input_data: AtomicInput, config: "TaskConfig") -> AtomicResult:
        self.found(raise_error=True)

        job_inputs = self.build_input(input_data, config)
        success, dexe = self.execute(job_inputs)

        if "INPUT HAS AT LEAST ONE SPELLING OR LOGIC MISTAKE" in dexe["stdout"]:
            raise InputError(dexe["stdout"])

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_data)

    def build_input(
        self, input_model: AtomicInput, config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        gamessrec = {
            "infiles": {},
            "scratch_directory": config.scratch_directory,
            "scratch_messy": config.scratch_messy,
        }

        opts = copy.deepcopy(input_model.keywords)

        def partition(total: float, fraction_replicated: float, ncores: int) -> Tuple[int, int]:
            """Return mwords and memddi values for ``total`` memory (already in mwords) partitioned so (0, 1]
            ``fraction_replicated`` is replicated over ``ncores`` and remainder distributed.
            """
            replicated_summed = total * fraction_replicated
            replicated_each = int(max(1, replicated_summed / ncores))
            distributed = int(total - replicated_summed)

            return replicated_each, distributed

        # Handle molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="gamess", units="Bohr", return_data=True)
        opts.update(moldata["keywords"])

        # Handle calc type and quantum chemical method
        opts.update(muster_modelchem(input_model.model.method, input_model.driver.derivative_int()))

        # Handle basis set
        if isinstance(input_model.model.basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use string basis name.")
        if input_model.model.basis is None:
            raise InputError("None for model.basis is not useable.")

        # * for gamess, usually insufficient b/c either ngauss or ispher needed
        opts["basis__gbasis"] = input_model.model.basis

        # Handle memory
        # * [GiB] --> [M QW]
        # * docs on mwords: "This is given in units of 1,000,000 words (as opposed to 1024*1024 words)"
        # * docs: "the memory required on each processor core for a run using p cores is therefore MEMDDI/p + MWORDS."
        # * int() rounds down
        mwords_total = int(config.memory * (1024 ** 3) / 8e6)

        for mem_frac_replicated in (1, 0.5, 0.1, 0.75):
            mwords, memddi = partition(mwords_total, mem_frac_replicated, config.ncores)
            trial_opts = copy.deepcopy(opts)
            trial_opts["contrl__exetyp"] = "check"
            trial_opts["system__parall"] = not (config.ncores == 1)
            trial_opts["system__mwords"] = mwords
            trial_opts["system__memddi"] = memddi
            trial_gamessrec = {
                "infiles": {"trial_gamess.inp": format_keywords(trial_opts) + molcmd},
                "command": [which("rungms"), "trial_gamess"],
                "scratch_messy": False,
                "scratch_directory": config.scratch_directory,
            }
            success, dexe = self.execute(trial_gamessrec)

            # this would be a lot cleaner if there was a unique or list of memory error strings
            if "INPUT HAS AT LEAST ONE SPELLING OR LOGIC MISTAKE" in dexe["stdout"]:
                raise InputError(dexe["stdout"])
            elif "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-" in dexe["stdout"]:
                pass
            else:
                opts["system__mwords"] = mwords
                opts["system__memddi"] = memddi
                break

        # TODO: "semget errno=ENOSPC -- check system limit for sysv semaphores." `ipcs -l`, can fix if too many gamess tests run at once by config.ncores = 4

        # ERROR: ROHF'S CCTYP MUST BE CCSD OR CR-CCL, WITH SERIAL EXECUTION
        if opts.get("contrl__cctyp", "").lower() == "ccsd" and opts.get("contrl__scftyp", "").lower() == "rohf":
            config.ncores = 1

        # Handle conversion from schema (flat key/value) keywords into local format
        optcmd = format_keywords(opts)

        gamessrec["infiles"]["gamess.inp"] = optcmd + molcmd
        gamessrec["command"] = [
            which("rungms"),
            "gamess",
            "00",
            str(config.ncores),
        ]  # rungms JOB VERNO NCPUS >& JOB.log &

        return gamessrec

        # Note decr MEMORY=100000 to get
        # ***** ERROR: MEMORY REQUEST EXCEEDS AVAILABLE MEMORY
        # to test gms fail

        # $CONTRL SCFTYP=ROHF MULT=3 RUNTYP=GRADIENT COORD=CART $END
        # $SYSTEM TIMLIM=1 MEMORY=800000 $END
        # $SCF    DIRSCF=.TRUE. $END
        # $BASIS  GBASIS=STO NGAUSS=2 $END
        # $GUESS  GUESS=HUCKEL $END
        # $DATA
        # Methylene...3-B-1 state...ROHF/STO-2G
        # Cnv  2
        #
        # Hydrogen   1.0    0.82884     0.7079   0.0
        # Carbon     6.0
        # Hydrogen   1.0   -0.82884     0.7079   0.0
        # $END

    def execute(self, inputs, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None):

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            ["gamess.dat"],
            scratch_messy=inputs["scratch_messy"],
            scratch_directory=inputs["scratch_directory"],
        )
        return success, dexe

    def parse_output(self, outfiles: Dict[str, str], input_model: AtomicInput) -> AtomicResult:

        # Get the stdout from the calculation (required)
        stdout = outfiles.pop("stdout")
        stderr = outfiles.pop("stderr")

        # gamessmol, if it exists, is dinky, just a clue to geometry of gamess results
        try:
            qcvars, gamesshess, gamessgrad, gamessmol, module = harvest(input_model.molecule, input_model.model.method, stdout, **outfiles)
        except Exception as e:
            raise UnknownError(
                "STDOUT:\n"
                + stdout
                + "\nSTDERR:\n"
                + stderr
                + "\nTRACEBACK:\n"
                + "".join(traceback.format_exception(*sys.exc_info()))
            )

        try:
            if gamessgrad is not None:
                qcvars[f"{input_model.model.method.upper()[4:]} TOTAL GRADIENT"] = gamessgrad
                qcvars["CURRENT GRADIENT"] = gamessgrad

            if gamesshess is not None:
                qcvars[f"{input_model.model.method.upper()[4:]} TOTAL HESSIAN"] = gamesshess
                qcvars["CURRENT HESSIAN"] = gamesshess

            if input_model.driver.upper() == "PROPERTIES":
                retres = qcvars[f"CURRENT ENERGY"]
            else:
                retres = qcvars[f"CURRENT {input_model.driver.upper()}"]
        except KeyError as e:
            raise UnknownError(
                "STDOUT:\n"
                + stdout
                + "\nSTDERR:\n"
                + stderr
                + "\nTRACEBACK:\n"
                + "".join(traceback.format_exception(*sys.exc_info()))
            )

        build_out(qcvars)
        atprop = build_atomicproperties(qcvars)

        provenance = Provenance(creator="GAMESS", version=self.get_version(), routine="rungms").dict()
        if module is not None:
            provenance["module"] = module

        output_data = {
            "schema_version": 1,
            "molecule": gamessmol,
            "extras": {"outfiles": outfiles, **input_model.extras},
            "properties": atprop,
            "provenance": provenance,
            "return_result": retres,
            "stderr": stderr,
            "stdout": stdout,
            "success": True,
        }

        # got to even out who needs plump/flat/Decimal/float/ndarray/list
        # * formerly unnp(qcvars, flat=True).items()
        output_data["extras"]["qcvars"] = {
            k.upper(): str(v) if isinstance(v, Decimal) else v for k, v in qcvars.items()
        }

        return AtomicResult(**{**input_model.dict(), **output_data})
