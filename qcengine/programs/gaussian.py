"""
Calls the GAUSSIAN excutable
"""

import os
import re
import tempfile
import warnings
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from qcelemental import constants
from qcelemental.models import AtomicInput, AtomicResult, Molecule, Provenance
from qcelemental.molparse import regex
from qcelemental.util import parse_version, safe_version, which

from qcengine.config import TaskConfig, get_config

from ..exceptions import InputError, UnknownError
from ..util import disk_files, execute, temporary_directory
from .model import ProgramHarness

class GaussianHarness(ProgramHarness):

    _defaults = {
        "name": "Gaussian",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass


    def found(self, raise_error: bool = False) -> bool:
        return which(
            "g09",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install Gaussian. Check it's in your PATH with `which g09`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        # Get the node configuration
        #config = get_config()

        which_prog = which("g09")
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "v.com"], {"v.com": ""})

            #mobj = re.search(r"Gaussian\s+([\d.]+)\s+for", exc["stdout"])
            #if not mobj:
            #    mobj = re.search(r"Gaussian version:\s+([\d.]+)\s+", exc["stdout"])

            #if mobj:
            #    self.version_cache[which_prog] = safe_version(mobj.group(1))

            # if "QC not defined" in exc["stdout"]:
            #else:
            #    return safe_version("09")

        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: TaskConfig) -> "AtomicResult":
        """
        Run Gaussian
        """
        # Check if Gaussian executable is found
        self.found(raise_error=True)

        # Setup the job
        job_inputs = self.build_input(input_model, config)

        # Run Gaussian
        exe_success, proc = self.execute(job_inputs)

        # Determine whether the calculation succeeded
        print("SUCCESS", exe_success)
        if exe_success:
            # If execution succeeded, collect results
            return self.parse_output(proc["outfiles"], input_model)


    def build_input(
        self, input_model: AtomicInput, config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        gaussian_ret = {
            "infiles": {},
            "commands": [which("g09"),  "input.inp"],
            "scratch_directory": config.scratch_directory
            }

        input_file = '''%mem=20MW
#P HF/6-31G(d) scf=tight

test1 HF/6-31G(d) sp formaldehyde

0 1
C1
O2  1  r2
H3  1  r3  2  a3
H4  1  r4  2  a4  3  d4

r2=1.20
r3=1.0
r4=1.0
a3=120.
a4=120.
d4=180.

        '''
        gaussian_ret['infiles']["input.inp"] = input_file

        return gaussian_ret

    def execute(self, inputs, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None):

        success, dexe = execute(
        inputs["commands"],
        inputs["infiles"],
        )
        print(dexe['proc'].stdout.read())
        print(dexe['proc'].stderr.read())
        return success, dexe

    def parse_output(self, outfiles: Dict[str, str], input_model: AtomicInput) -> AtomicResult:
        print("IN PARSE OUTPUT")

        #output_data = {}
        stdout = outfiles.pop('stdout')
        stderr = outfiles.pop('stderr')

        method = input_model.model.method.lower()
        #method = method[4:] if method.startswith("") else method

        try:
            qcvars, gaussianmol, module = harvest(
                input_model.molecule, method, stdout, **outfiles
            )
        except:
            pass

        #properties = {
        #    "nuclear_repulsion_energy": bdata['99.0'][0],
        #    "scf_total_energy": bdata["99.0"][1],
        #    "return_energy": bdata["99.0"][-1],
        #}

        qcvars = {}

        #props, prov = self._parse_logfile_common(outtext, input_model.dict())
        #output_data["provenance"] = prov
        #output_data["properties"] = properties
        #output_data["properties"].update(props)
        output_data["stdout"] = stdout
        output_data["success"] = True

        merged_data = {**input_model.dict(), **output_data}
        #merged_data["extras"]["qcvars"] = qcvars
        print(merged_data)

        return AtomicResult(**merged_data)

    # DELETE THIS FUNC
    def _parse_logfile_common(self, outtext: str, input_dict: Dict[str, Any]):
        """
        Parses fields from log file that are not parsed from QCSCRATCH in parse_output
        """

        properties = {}
        provenance = Provenance(creator="QChem", version=self.get_version(), routine="qchem").dict()
        mobj = re.search(r"This is a multi-thread run using ([0-9]+) threads", outtext)
        if mobj:
            provenance["nthreads"] = int(mobj.group(1))

        mobj = re.search(r"Total job time:\s*" + NUMBER + r"s\(wall\)", outtext)
        if mobj:
            provenance["wall_time"] = float(mobj.group(1))

        mobj = re.search(r"Archival summary:\s*\n[0-9]+\\[0-9+]\\([\w\.\-]+)\\", outtext)
        if mobj:
            provenance["hostname"] = mobj.group(1)

        mobj = re.search(r"\n\s*There are\s+(\d+) alpha and\s+(\d+) beta electrons\s*\n", outtext)
        if mobj:
            properties["calcinfo_nalpha"] = int(mobj.group(1))
            properties["calcinfo_nbeta"] = int(mobj.group(2))

        mobj = re.search(r"\n\s*There are\s+\d+ shells and\s+(\d+) basis functions\s*\n", outtext)
        if mobj:
            properties["calcinfo_nbasis"] = int(mobj.group(1))

        mobj = re.search(r"\n\s*RI-MP2 CORRELATION ENERGY\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj:
            properties["mp2_correlation_energy"] = float(mobj.group(1))

        mobj = re.search(r"\n\s*RI-MP2 SINGLES ENERGY\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj:
            properties["mp2_singles_energy"] = float(mobj.group(1))

        mobj_aaaa = re.search(r"\n\s*RI-MP2 ENERGY \(aa\|aa\)\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        mobj_bbbb = re.search(r"\n\s*RI-MP2 ENERGY \(bb\|bb\)\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj_aaaa and mobj_bbbb:
            properties["mp2_same_spin_correlation_energy"] = float(mobj_aaaa.group(1)) + float(mobj_bbbb.group(1))

        mobj_aabb = re.search(r"\n\s*RI-MP2 ENERGY \(aa\|bb\)\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        mobj_bbaa = re.search(r"\n\s*RI-MP2 ENERGY \(bb\|aa\)\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj_aaaa and mobj_bbbb:
            properties["mp2_opposite_spin_correlation_energy"] = float(mobj_aabb.group(1)) + float(mobj_bbaa.group(1))

        properties["calcinfo_natom"] = len(input_dict["molecule"]["symbols"])

        mobj = re.search(r"\n\s*(\d+)\s+" + NUMBER + r"\s+" + NUMBER + r"\s+Convergence criterion met\s*\n", outtext)
        if mobj:
            properties["scf_iterations"] = int(mobj.group(1))

        mobj = re.search(
            r"\n\s+Dipole Moment \(Debye\)\s*\n\s+X\s+" + NUMBER + r"\s+Y\s+" + NUMBER + r"\s+Z\s+" + NUMBER + r"\s*\n",
            outtext,
        )
        if mobj:
            cf = constants.conversion_factor("debye", "e * bohr")
            properties["scf_dipole_moment"] = [float(mobj.group(i)) * cf for i in range(1, 4)]

        return properties, provenance

    def parse_logfile(self, outfiles: Dict[str, str]) -> AtomicResult:
        """
        Parses a log file.
        """
        warnings.warn(
            "parse_logfile will result in precision loss for some fields due to trunctation in " "Q-Chem output files."
        )

        outtext = outfiles["dispatch.log"]

        mobj = re.search(r"(?:User input\:|Running Job?)\s+\d+\s+of\s+(\d+)", outtext)
        if mobj:
            if int(mobj.group(1)) > 1:
                raise InputError("Multi-job Q-Chem log files not supported.")

        input_dict = {}
        mobj = re.search(r"\n-{20,}\nUser input:\n-{20,}\n(.+)\n-{20,}", outtext, re.DOTALL)
        if mobj:
            inputtext = mobj.group(1)

            rem_match = re.search(r"\$rem\s*\n([^\$]+)\n\s*\$end", inputtext, re.DOTALL | re.IGNORECASE)
            if rem_match:
                rem_text = rem_match.group(1)
                lines = rem_text.split("\n")
                keywords = {}
                for line in lines:
                    s = re.sub(r"(^|[^\\])!.*", "", line).split()
                    if len(s) == 0:
                        continue
                    keywords[s[0].lower()] = s[1].lower()
                input_dict["model"] = {}
                input_dict["model"]["method"] = keywords.pop("method").lower()
                input_dict["model"]["basis"] = keywords.pop("basis").lower()
                if "jobtype" in keywords:
                    jobtype = keywords.pop("jobtype")
                else:
                    jobtype = "sp"
                _jobtype_to_driver = {
                    "sp": "energy",
                    "force": "gradient",
                    "freq": "hessian",
                }  # optimization intentionally not supported
                try:
                    input_dict["driver"] = _jobtype_to_driver[jobtype]
                except KeyError:
                    raise KeyError(f"Jobtype {jobtype} not supported in qchem log file parser.")

                for key in keywords:
                    if keywords[key] == "false":
                        keywords[key] = False
                    if keywords[key] == "true":
                        keywords[key] = True
                input_dict["keywords"] = keywords

            molecule_match = re.search(r"\$molecule\s*\n([^\$]+)\n\s*\$end", inputtext, re.DOTALL | re.IGNORECASE)
            if molecule_match:
                molecule_text = molecule_match.group(1)
                if keywords.get("input_bohr", False):
                    molecule_text += "\nunits au"
                molecule = Molecule.from_data(molecule_text, dtype="psi4")
                input_dict["molecule"] = molecule.dict()

            _excluded_rem = {
                "input_bohr",
                "mem_total",
                "mem_static",
            }  # parts of the rem normally written by build_input, and which do not affect results
            for item in _excluded_rem:
                if item in input_dict["keywords"]:
                    input_dict["keywords"].pop(item)

        try:
            qcscr_result = self.parse_output(outfiles, AtomicInput(**input_dict)).dict()
        except KeyError:
            props, prov = self._parse_logfile_common(outtext, input_dict)
            qcscr_result = {"properties": props, "provenance": prov, **input_dict}

        mobj = re.search(r"\n\s*Total\s+energy in the final basis set =\s+" + NUMBER + r"\s*\n", outtext)
        if mobj and qcscr_result["properties"].get("scf_total_energy", None) is None:
            qcscr_result["properties"]["scf_total_energy"] = float(mobj.group(1))

        mobj = re.search(r"\n\s*Nuclear Repulsion Energy =\s+" + NUMBER + r"\s+hartrees\s*\n", outtext)
        if mobj and qcscr_result["properties"].get("nuclear_repulsion_energy", None) is None:
            qcscr_result["properties"]["nuclear_repulsion_energy"] = float(mobj.group(1))

        mobj = re.search(r"\n\s*RI-MP2 TOTAL ENERGY\s+=\s+" + NUMBER + r"\s+au\s*\n", outtext)
        if mobj and qcscr_result["properties"].get("mp2_total_energy", None) is None:
            qcscr_result["properties"]["mp2_total_energy"] = float(mobj.group(1))

        _mp2_methods = {"rimp2"}

        method = input_dict["model"]["method"].lower()
        if qcscr_result["properties"].get("return_energy", None) is None:
            if method in _mp2_methods:
                qcscr_result["properties"]["return_energy"] = qcscr_result["properties"]["mp2_total_energy"]
            elif method in _scf_methods:
                qcscr_result["properties"]["return_energy"] = qcscr_result["properties"]["scf_total_energy"]
            else:
                raise NotImplementedError(f"Method {method} not supported by logfile parser for energy driver.")

        if input_dict["driver"] == "gradient" and qcscr_result.get("return_result", None) is None:

            def read_matrix(text):
                lines = text.split("\n")
                i = 0
                mdict = defaultdict(dict)
                maxcol = 0
                maxrow = 0
                while i < len(lines):
                    cols = [int(idx) for idx in lines[i].split()]
                    maxcol = max(maxcol, *cols)
                    i += 1
                    while i < len(lines):
                        s = lines[i].split()
                        if len(s) <= len(cols):
                            break
                        assert len(s) == len(cols) + 1, s
                        row = int(s[0])
                        maxrow = max(maxrow, row)
                        data = [float(field) for field in s[1:]]
                        for col_idx, col in enumerate(cols):
                            mdict[row - 1][col - 1] = data[col_idx]
                        i += 1

                ret = np.zeros((maxrow, maxcol))
                for row in mdict:
                    for col in mdict[row]:
                        ret[row, col] = mdict[row][col]
                return ret

            if method in _mp2_methods:
                mobj = re.search(
                    r"\n\s+Full Analytical Gradient of MP2 Energy \(in au.\)\s*\n"
                    r"([\s\d\.EDed\+\-]+)\n"
                    r"\s*Gradient time:",
                    outtext,
                )

            elif method in _scf_methods:
                mobj = re.search(
                    r"\n\s+Gradient of SCF Energy\s*\n([\s\d\.EDed\+\-]+)\n\s*Max gradient component =", outtext
                )

            else:
                raise NotImplementedError(f"Method {method} not supported by the logfile parser for gradient driver.")

            if mobj:
                qcscr_result["return_result"] = read_matrix(mobj.group(1)).T

        qcscr_result["success"] = True  # XXX: have a nice day?
        qcscr_result["stdout"] = outtext
        if input_dict["driver"] == "energy" and qcscr_result.get("return_result", None) is None:
            qcscr_result["return_result"] = qcscr_result["properties"]["return_energy"]

        return AtomicResult(**qcscr_result)
