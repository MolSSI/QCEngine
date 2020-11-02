"""
Calls adcc
"""
import json
import os
import sys
from pathlib import Path
import importlib
from typing import TYPE_CHECKING, Any, Dict, Optional, Tuple
import numpy as np

from qcelemental.models import AtomicResult, AtomicResultProperties, Provenance
from qcelemental.util import deserialize, parse_version, safe_version, which_import

from ..exceptions import InputError, RandomError, ResourceError, UnknownError
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


class AdccHarness(ProgramHarness):

    _defaults = {
        "name": "adcc",
        "scratch": False,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Whether adcc harness is ready for operation.
        Parameters
        ----------
        raise_error: bool
            Passed on to control negative return between False and ModuleNotFoundError raised.
        Returns
        -------
        bool
            If adcc is found, returns True.
            If raise_error is False and adcc is missing, returns False.
            If raise_error is True and adcc is missing, the error message is raised.
        """
        return which_import(
            "adcc",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install adcc -c adcc`.",
        )

    def get_version(self) -> str:
        """Return the currently used version of adcc"""
        self.found(raise_error=True)

        which_prog = which_import("adcc")
        if which_prog not in self.version_cache:
            import adcc

            self.version_cache[which_prog] = safe_version(adcc.__version__)
        return self.version_cache[which_prog]

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs adcc
        """
        self.found(raise_error=True)
        import adcc

        mol = input_model.molecule
        model = input_model.model
        kws = input_model.keywords
        xyz = mol.to_string(dtype="xyz+", units="Bohr")
        xyz = "\n".join(xyz.split("\n")[2:])
        scfres = adcc.backends.run_hf(
            backend=None,  # auto-select available backend
            xyz=xyz,
            multiplicity=mol.molecular_multiplicity,
            charge=mol.molecular_charge,
            basis=model.basis,
            #   conv_tol=model.scf_conv,
            #   max_iter=,
        )

        # handle defaults nicely...
        n_singlets = kws.pop("n_singlets", 3)

        adcc.set_n_threads(config.ncores)
        adcc_state = adcc.run_adc(scfres, method=model.method, n_singlets=n_singlets)

        input_data = input_model.dict(encoding="json")
        output_data = input_data.copy()
        output_data["success"] = adcc_state.converged
        provenance = Provenance(creator="adcc", version=self.get_version()).dict()
        provenance["nthreads"] = adcc.get_n_threads()
        output_data["provenance"] = provenance
        output_data["properties"] = AtomicResultProperties()
        output_data["return_result"] = adcc_state.excitation_energy

        extract_props = input_model.driver == "properties"
        output_data["extras"]["qcvars"] = self._extract_qcvars(adcc_state, extract_props)
        return AtomicResult(**output_data)

    def _extract_qcvars(self, state, extract_props=False):
        qcvars = {}
        name = state.method.name
        NAME = name.upper()
        is_cvs_adc3 = state.method.level >= 3 and state.ground_state.has_core_occupied_space
        mp = state.ground_state
        mp_energy = mp.energy(state.method.level if not is_cvs_adc3 else 2)
        mp_corr = 0.0
        qcvars[f"HF TOTAL ENERGY"] = mp.reference_state.energy_scf
        if state.method.level > 1:
            for level in range(2, state.method.level + 1):
                if level >= 3 and is_cvs_adc3:
                    continue
                energy = mp.energy_correction(level)
                mp_corr += energy
                qcvars[f"MP{level} CORRELATION ENERGY"] = energy
                qcvars[f"MP{level} TOTAL ENERGY"] = mp.energy(level)
        qcvars["EXCITATION KIND"] = state.kind.upper()
        qcvars[f"{NAME} ITERATIONS"] = state.n_iter
        qcvars[NAME + " EXCITATION ENERGIES"] = state.excitation_energy
        qcvars["NUMBER OF EXCITED STATES"] = len(state.excitation_energy)
        if extract_props:
            qcvars["HF DIPOLE"] = mp.reference_state.dipole_moment
            if state.method.level > 1:
                qcvars["MP2 DIPOLE"] = mp.dipole_moment(2)
            # transition properties
            qcvars[f"{NAME} TRANSITION DIPOLES (LEN)"] = state.transition_dipole_moment
            qcvars[f"{NAME} TRANSITION DIPOLES (VEL)"] = state.transition_dipole_moment_velocity
            qcvars[f"{NAME} OSCILLATOR STRENGTHS (LEN)"] = state.oscillator_strength
            qcvars[f"{NAME} OSCILLATOR STRENGTHS (VEL)"] = state.oscillator_strength_velocity
            qcvars[f"{NAME} ROTATIONAL STRENGTHS (VEL)"] = state.rotatory_strength
            # state properties
            qcvars[f"{NAME} STATE DIPOLES"] = state.state_dipole_moment
        return qcvars
