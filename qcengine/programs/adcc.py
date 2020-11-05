"""
Calls adcc
"""
from typing import Dict, TYPE_CHECKING

from qcelemental.util import safe_version, which_import
from qcelemental.models import (AtomicResult, AtomicResultProperties,
                                Provenance)

from .model import ProgramHarness
from ..exceptions import InputError, UnknownError

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig

#


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
            If adcc and psi4 are found, returns True.
            If raise_error is False and adcc or psi4 is missing, returns False.
            If raise_error is True and adcc or psi4 are missing, the error message is raised.
        """
        found_adcc = which_import(
            "adcc",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install adcc -c adcc`.",
        )
        found_psi4 = which_import(
            "psi4",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install psi4 for adcc harness via `conda install psi4 -c psi4`.",
        )
        return found_adcc and found_psi4

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
        import psi4

        mol = input_model.molecule
        model = input_model.model
        conv_tol = input_model.keywords.get("conv_tol", 1e-6)

        psi4_molecule = psi4.core.Molecule.from_schema(dict(mol.dict(), fix_symmetry="c1"))
        psi4.core.clean()
        psi4.core.be_quiet()
        psi4.set_options(
            {
                "basis": model.basis,
                "scf_type": "pk",
                "e_convergence": conv_tol / 100,
                "d_convergence": conv_tol / 10,
                # 'maxiter': max_iter,
                "reference": "RHF" if mol.molecular_multiplicity == 1 else "UHF",
            }
        )
        _, wfn = psi4.energy("HF", return_wfn=True, molecule=psi4_molecule)

        adcc.set_n_threads(config.ncores)
        try:
            adcc_state = adcc.run_adc(wfn, method=model.method, **input_model.keywords)
        except ValueError as e:
            # TODO Once the new adcc is realeased switch to "adcc.InputError as e"
            raise InputError(str(e))
        except Exception as e:
            raise UnknownError(str(e))

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
        mp_corr = 0.0
        qcvars["HF TOTAL ENERGY"] = mp.reference_state.energy_scf
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
