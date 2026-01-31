"""
Calls adcc
"""
from typing import TYPE_CHECKING, Any, ClassVar, Dict

from qcelemental.models.v2 import AtomicResult, BasisSet, Provenance
from qcelemental.util import safe_version, which_import

from ..exceptions import InputError, UnknownError
from .model import ProgramHarness
from .qcvar_identities_resources import build_atomicproperties

if TYPE_CHECKING:
    from qcelemental.models.v2 import AtomicInput

    from ..config import TaskConfig

#


class AdccHarness(ProgramHarness):
    """Interface for adcc project."""

    _defaults: ClassVar[Dict[str, Any]] = {
        "name": "adcc",
        "scratch": False,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

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
            raise_msg="Please install via `conda install adcc -c conda-forge`.",
        )
        found_psi4 = which_import(
            "psi4",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install psi4 for adcc harness via `conda install psi4 -c conda-forge`.",
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
        model = input_model.specification.model
        conv_tol = input_model.specification.keywords.get("conv_tol", 1e-6)

        if input_model.specification.driver not in ["energy", "properties"]:
            raise InputError(f"Driver {input_model.specification.driver} not implemented for ADCC.")

        if isinstance(input_model.specification.model.basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use string basis name.")
        if not input_model.specification.model.basis:
            raise InputError("Model must contain a basis set.")

        psi4_molecule = psi4.core.Molecule.from_schema(dict(mol.model_dump(), fix_symmetry="c1"))
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
        compute_success = False
        try:
            adcc_state = adcc.run_adc(wfn, method=model.method, **input_model.specification.keywords)
            compute_success = adcc_state.converged
        except adcc.InputError as e:
            raise InputError(str(e))
        except Exception as e:
            raise UnknownError(str(e))

        input_data = input_model.model_dump(encoding="json")
        output_data = {"input_data": input_data, "extras": {}, "molecule": mol}
        output_data["success"] = compute_success

        if compute_success:
            output_data["return_result"] = adcc_state.excitation_energy[0]

            extract_props = input_model.specification.driver == "properties"
            qcvars = adcc_state.to_qcvars(recurse=True, properties=extract_props)
            qcvars["CURRENT ENERGY"] = adcc_state.excitation_energy[0]
            atprop = build_atomicproperties(qcvars)
            output_data["extras"]["qcvars"] = qcvars
            output_data["properties"] = atprop

        provenance = Provenance(creator="adcc", version=self.get_version(), routine="adcc").model_dump()
        provenance["nthreads"] = adcc.get_n_threads()
        output_data["provenance"] = provenance

        return AtomicResult(**output_data)
