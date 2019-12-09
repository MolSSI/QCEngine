"""Compute quantum chemistry using Iowa State's GAMESS executable."""

import copy
import pprint
from decimal import Decimal
from typing import Any, Dict, Optional

import qcelemental as qcel
from qcelemental.models import AtomicResult
from qcelemental.util import safe_version, unnp, which

from ...exceptions import InputError
from ...util import execute
from ..model import ProgramHarness
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

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
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
        self, input_model: "AtomicInput", config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        gamessrec = {"infiles": {}, "scratch_directory": config.scratch_directory}

        opts = copy.deepcopy(input_model.keywords)

        # Handle memory
        # for gamess, [GiB] --> [MW]
        opts["system__mwords"] = int(config.memory * (1024 ** 3) / 8e6)

        # Handle molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="gamess", units="Bohr", return_data=True)
        opts.update(moldata["keywords"])

        # Handle calc type and quantum chemical method
        opts.update(muster_modelchem(input_model.model.method, input_model.driver.derivative_int()))

        # Handle basis set
        # * for gamess, usually insufficient b/c either ngauss or ispher needed
        opts["basis__gbasis"] = input_model.model.basis

        # Handle conversion from schema (flat key/value) keywords into local format
        optcmd = format_keywords(opts)

        gamessrec["infiles"]["gamess.inp"] = optcmd + molcmd
        gamessrec["command"] = [which("rungms"), "gamess"]  # rungms JOB VERNO NCPUS >& JOB.log &

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
            inputs["command"], inputs["infiles"], [], scratch_messy=False, scratch_directory=inputs["scratch_directory"]
        )
        return success, dexe

    def parse_output(self, outfiles: Dict[str, str], input_model: "AtomicInput") -> "AtomicResult":

        # gamessmol, if it exists, is dinky, just a clue to geometry of gamess results
        qcvars, gamessgrad, gamessmol = harvest(input_model.molecule, outfiles["stdout"])

        if gamessgrad is not None:
            qcvars["CURRENT GRADIENT"] = gamessgrad

        qcvars = unnp(qcvars, flat=True)

        output_data = {
            "schema_name": "qcschema_output",
            "molecule": gamessmol,
            "schema_version": 1,
            "extras": {},
            "properties": {"nuclear_repulsion_energy": gamessmol.nuclear_repulsion_energy()},
            "return_result": qcvars[f"CURRENT {input_model.driver.upper()}"],
            "stdout": outfiles["stdout"],
        }

        # got to even out who needs plump/flat/Decimal/float/ndarray/list
        output_data["extras"]["qcvars"] = {
            k.upper(): float(v) if isinstance(v, Decimal) else v for k, v in qcel.util.unnp(qcvars, flat=True).items()
        }

        # copy qcvars into schema where possible
        qcvars_to_properties = {
            "DFT XC ENERGY": "scf_xc_energy",
            "ONE-ELECTRON ENERGY": "scf_one_electron_energy",
            "TWO-ELECTRON ENERGY": "scf_two_electron_energy",
            "SCF TOTAL ENERGY": "scf_total_energy",
            "MP2 CORRELATION ENERGY": "mp2_correlation_energy",
            "MP2 TOTAL ENERGY": "mp2_total_energy",
            "CCSD CORRELATION ENERGY": "ccsd_correlation_energy",
            "CCSD TOTAL ENERGY": "ccsd_total_energy",
            "CCSD(T) CORRELATION ENERGY": "ccsd_prt_pr_correlation_energy",
            "CCSD(T) TOTAL ENERGY": "ccsd_prt_pr_total_energy",
        }
        for qcvar in qcvars:
            if qcvar in qcvars_to_properties:
                output_data["properties"][qcvars_to_properties[qcvar]] = qcvars[qcvar]
        if {"SCF DIPOLE X", "SCF DIPOLE Y", "SCF DIPOLE Z"} & set(qcvars.keys()):
            conv = Decimal(qcel.constants.conversion_factor("debye", "e * bohr"))
            output_data["properties"]["scf_dipole_moment"] = [
                qcvars["SCF DIPOLE X"] * conv,
                qcvars["SCF DIPOLE Y"] * conv,
                qcvars["SCF DIPOLE Z"] * conv,
            ]
        output_data["success"] = True

        return AtomicResult(**{**input_model.dict(), **output_data})
