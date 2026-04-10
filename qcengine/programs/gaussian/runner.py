"""Compute quantum chemistry using Gaussian 16 or Gaussian 09."""

import re
from decimal import Decimal
from typing import Any, Dict, Optional, Tuple

import numpy as np
import qcelemental as qcel
from qcelemental.models.v2 import AtomicInput, AtomicResult, BasisSet, Provenance
from qcelemental.util import safe_version, which, which_import

from ...exceptions import InputError, ResourceError, UnknownError
from ...util import execute
from ..model import ProgramHarness
from ..qcvar_identities_resources import build_atomicproperties, build_out
from ..util import error_stamp
from .germinate import muster_modelchem
from .harvester import harvest, is_normal_termination
from .keywords import build_com_file, build_route_line

# rq-cb0f61f2
_GAUSSIAN_EXECUTABLES = ("g16", "g09")

# Regex to extract version from Gaussian startup output.
# rq-cb0f61f2
_VERSION_RE = re.compile(r"Gaussian\s+(\d+),\s+Revision\s+([\w.]+),")


def _find_gaussian() -> Optional[str]:
    """Return the path of the first Gaussian executable found on PATH, or None.

    rq-cb0f61f2 rq-d3b8d8c3
    """
    for exe in _GAUSSIAN_EXECUTABLES:
        path = which(exe, return_bool=False)
        if path:
            return path
    return None


# rq-998da273 rq-436a54d5
class GaussianHarness(ProgramHarness):

    # rq-998da273
    _defaults = {
        "name": "Gaussian",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    # rq-4cf7c460 rq-d3b8d8c3 rq-8e21b7d7 rq-b5e72c6d rq-bd98d417 rq-242ef1d8
    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Return True iff both a Gaussian executable and cclib are available.

        Checks for g16 first, then g09.  Also requires cclib to be importable.
        """
        qc = _find_gaussian() is not None
        if not qc:
            if raise_error:
                raise ModuleNotFoundError(
                    "Gaussian executable (g16 or g09) not found on PATH. "
                    "Install Gaussian and ensure the executable is on your PATH."
                )
            return False

        dep = which_import(
            "cclib",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="For the Gaussian harness, please install cclib: `pip install cclib`.",
        )
        return dep

    # rq-186c5060 rq-850ede4f rq-3e134642 rq-20d9a784
    def get_version(self) -> str:
        """Return a normalised Gaussian version string (e.g. "16.C.02").

        Runs the Gaussian executable with a blank .com input file and parses the
        version from the startup banner.  Result is cached per executable path.
        """
        self.found(raise_error=True)

        exe_path = _find_gaussian()
        if exe_path in self.version_cache:
            return self.version_cache[exe_path]

        # Run Gaussian with an empty .com file to provoke the version banner.
        # Gaussian writes the version line to the log file even when the input
        # is blank (causing a segfault), so collect version_probe.log and
        # accept any exit code by passing a large exit_code threshold.
        success, output = execute(
            [exe_path, "version_probe.com"],
            {"version_probe.com": "\n"},
            ["version_probe.log"],
            exit_code=255,
        )
        log_text = (output.get("outfiles") or {}).get("version_probe.log", "") or ""
        stdout = output.get("stdout", "") or ""
        stderr = output.get("stderr", "") or ""
        combined = log_text + stdout + stderr

        m = _VERSION_RE.search(combined)
        if not m:
            raise UnknownError(
                f"Could not parse Gaussian version from output.\n"
                f"stdout: {stdout}\nstderr: {stderr}"
            )

        version_str = f"{m.group(1)}.{m.group(2)}"
        self.version_cache[exe_path] = version_str
        return version_str

    # rq-99807237 rq-e62578e2 rq-4d683dac rq-1a051bc2
    # rq-c77a2c13 rq-165224fc rq-edb73b83 rq-a4c449da rq-71aceb93
    def compute(self, input_model: AtomicInput, config: "TaskConfig") -> AtomicResult:
        """Top-level compute entry point for the Gaussian harness."""
        self.found(raise_error=True)

        # rq-c77a2c13
        if isinstance(input_model.specification.model.basis, BasisSet):
            raise InputError(
                "QCSchema BasisSet for model.basis is not supported by the Gaussian harness. "
                "Use a string basis name (e.g. 'STO-3G')."
            )

        job_inputs = self.build_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        log_text = dexe["outfiles"].get("gaussian.log") or ""
        stdin = job_inputs["infiles"].get("gaussian.com", "")
        stdout = dexe.get("stdout", "") or ""
        stderr = dexe.get("stderr", "") or ""

        # rq-a4c449da
        if "galloc: could not allocate memory" in log_text:
            raise ResourceError(error_stamp(stdin, log_text, stderr))

        # rq-165224fc
        if (
            "Ernie: is not defined as a Gaussian input file" in log_text
            or "Ernie: unrecognized" in log_text
        ):
            raise InputError(error_stamp(stdin, log_text, stderr))

        # rq-edb73b83
        if "Error in internal coordinate system" in log_text:
            raise InputError(error_stamp(stdin, log_text, stderr))

        # rq-71aceb93 — absence of normal termination with no known pattern
        if not is_normal_termination(log_text):
            raise UnknownError(error_stamp(stdin, log_text, stderr))

        dexe["outfiles"]["stdout"] = stdout
        dexe["outfiles"]["stderr"] = stderr
        dexe["outfiles"]["input"] = stdin
        return self.parse_output(dexe["outfiles"], input_model)

    # rq-9cdf3f7d rq-74e0110c rq-0052a2c7
    def build_input(
        self,
        input_model: AtomicInput,
        config: "TaskConfig",
        template: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Construct the Gaussian .com input file and execution record."""
        # rq-0aa99076
        if isinstance(input_model.specification.model.basis, BasisSet):
            raise InputError(
                "QCSchema BasisSet for model.basis is not supported by the Gaussian harness."
            )
        if input_model.specification.model.basis is None:
            raise InputError("model.basis must be a string basis name for the Gaussian harness.")

        # rq-5053204b
        if not all(input_model.molecule.real):
            raise InputError("The Gaussian harness does not support ghost atoms (real=False).")

        mol = input_model.molecule
        method = input_model.specification.model.method
        basis = input_model.specification.model.basis
        driver = input_model.specification.driver.value  # DriverEnum → plain lowercase string

        # rq-ef62979b
        modelchem = muster_modelchem(method, driver, mol.molecular_multiplicity)

        # rq-74e0110c — build atom block in Ångström with 10+ decimal places
        bohr_to_ang = qcel.constants.conversion_factor("bohr", "angstrom")
        geom_ang = mol.geometry.reshape(-1, 3) * bohr_to_ang
        atom_lines = [
            f"{sym:2s}  {x:20.10f}  {y:20.10f}  {z:20.10f}"
            for sym, (x, y, z) in zip(mol.symbols, geom_ang)
        ]
        atom_block = "\n".join(atom_lines)

        # rq-db92093a rq-4ad16d95
        mem_gb = max(1, int(config.memory))
        link0 = {
            "NProcShared": config.ncores,
            "Mem": f"{mem_gb}GB",
        }

        route_line = build_route_line(
            modelchem["method_string"],
            basis,
            modelchem["job_type"],
            modelchem["extra_keywords"],
            input_model.specification.keywords,
        )

        com_text = build_com_file(
            link0=link0,
            route_line=route_line,
            title="QCEngine Gaussian calculation",
            charge=int(mol.molecular_charge),
            multiplicity=int(mol.molecular_multiplicity),
            atom_block=atom_block,
        )

        exe_path = _find_gaussian()

        return {
            "infiles": {"gaussian.com": com_text},
            "command": [exe_path, "gaussian.com"],
            "scratch_directory": config.scratch_directory,
            "scratch_messy": config.scratch_messy,
        }

    def execute(
        self,
        inputs: Dict[str, Any],
        extra_outfiles: Optional[list] = None,
        extra_commands: Optional[list] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Dict]:
        """Run the Gaussian executable; collect gaussian.log."""
        # rq-cb0f61f2
        # Gaussian's exit code is unreliable; we always inspect the log.
        # Pass exit_code=1 so execute() returns success=True for both 0 and 1.
        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            ["gaussian.log"],
            scratch_messy=inputs["scratch_messy"],
            scratch_directory=inputs["scratch_directory"],
            exit_code=1,
        )
        return success, dexe

    def parse_output(
        self,
        outfiles: Dict[str, str],
        input_model: AtomicInput,
    ) -> AtomicResult:
        """Parse Gaussian output files and return an AtomicResult."""
        stdout = outfiles.pop("stdout", "") or ""
        stderr = outfiles.pop("stderr", "") or ""
        input_text = outfiles.pop("input", "") or ""
        log_text = outfiles.get("gaussian.log", "") or ""

        method = input_model.specification.model.method.lower()
        # Strip any "gaussian-" prefix in case caller added one
        if method.startswith("gaussian-"):
            method = method[len("gaussian-"):]

        try:
            qcvars, grad, hess, out_mol = harvest(input_model.molecule, method, log_text)
        except Exception as e:
            raise UnknownError(error_stamp(input_text, log_text, stderr)) from e

        method_upper = method.upper()

        if grad is not None:
            qcvars[f"{method_upper} TOTAL GRADIENT"] = grad
            qcvars["CURRENT GRADIENT"] = grad

        if hess is not None:
            qcvars[f"{method_upper} TOTAL HESSIAN"] = hess
            qcvars["CURRENT HESSIAN"] = hess

        driver = input_model.specification.driver.value  # DriverEnum → plain lowercase string
        try:
            if driver.upper() == "PROPERTIES":
                retres = float(qcvars["CURRENT ENERGY"])
            else:
                retres = qcvars[f"CURRENT {driver.upper()}"]
        except KeyError:
            raise UnknownError(error_stamp(input_text, log_text, stderr))

        # Serialise Decimal → float/str so JSON encoding works downstream
        qcvars_serialised = {
            k.upper(): (str(v) if isinstance(v, Decimal) else v)
            for k, v in qcvars.items()
        }

        build_out(qcvars)
        atprop = build_atomicproperties(qcvars)

        provenance = Provenance(
            creator="Gaussian",
            version=self.get_version(),
            routine="gaussian",
        ).model_dump()

        output_data = {
            "schema_version": 2,
            "input_data": input_model,
            "molecule": out_mol,
            "extras": {**input_model.specification.extras, "qcvars": qcvars_serialised},
            "native_files": {
                "gaussian.com": input_text,
                "gaussian.log": log_text,
            },
            "properties": atprop,
            "provenance": provenance,
            "return_result": retres,
            "stderr": stderr,
            "stdout": stdout,
            "success": True,
        }

        return AtomicResult(**output_data)
