"""
Calls the TorchANI package.
"""

from typing import TYPE_CHECKING, Any, ClassVar, Dict

from qcelemental.models.v2 import AtomicResult, Provenance
from qcelemental.util import parse_version, safe_version, which_import

from ..exceptions import InputError, ResourceError
from ..units import ureg
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models.v2 import AtomicInput

    from ..config import TaskConfig


class TorchANIHarness(ProgramHarness):
    """Interface for TorchANI project."""

    _CACHE = {}

    _defaults: ClassVar[Dict[str, Any]] = {
        "name": "TorchANI",
        "scratch": False,
        "thread_safe": True,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": False,
    }
    version_cache: Dict[str, str] = {}

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which_import(
            "torchani",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `pip install torchani`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("torchani")
        if which_prog not in self.version_cache:
            import torchani

            self.version_cache[which_prog] = safe_version(torchani.__version__)

        return self.version_cache[which_prog]

    def get_model(self, name: str) -> "torchani.models.BuiltinModels":
        name = name.lower()

        if name in self._CACHE:
            return self._CACHE[name]

        import torch
        import torch.nn as nn
        import torchani

        class EnsembleEnergies(nn.Module):
            """
            Wrap a TorchANI model and return (species, concatenated_member_energies).

            Compatible with both call styles:
              - m((species, coordinates))
              - m(species, coordinates, *extra_args)
            """
            def __init__(self, base_model):
                super().__init__()
                self.base = base_model

                # expose `model.species` for v2.2 for symbol->index conversion
                if hasattr(base_model, "species"):
                    self.species = base_model.species

            @staticmethod
            def _normalize_call_args(first, rest):
                # forward((species, coords), ...)
                if isinstance(first, (tuple, list)) and len(first) == 2:
                    species, coords = first
                    extra = rest
                else:
                    # forward(species, coords, ...)
                    species = first
                    coords = rest[0]
                    extra = rest[1:]
                return species, coords, extra

            @staticmethod
            def _extract_energies(member_result):
                # TorchANI 2.7+: returns object with .energies
                if hasattr(member_result, "energies"):
                    return member_result.energies
                # TorchANI 2.2.x (and older style): returns (species, energies)
                if isinstance(member_result, (tuple, list)) and len(member_result) >= 2:
                    return member_result[1]
                raise TypeError(f"Unrecognized TorchANI model output: {type(member_result)!r}")

            def forward(self, first, *rest, **kwargs):
                species, coords, extra = self._normalize_call_args(first, rest)

                member_energies = []
                n_members = len(self.base)
                for i in range(n_members):
                    member = self.base[i]
                    out = member((species, coords), *extra, **kwargs)
                    e = self._extract_energies(out)
                    member_energies.append(e)

                # keep old flat vector pattern using torch.cat([...])
                outputs = torch.cat(member_energies, dim=0)
                return species, outputs

        ani_models = {
            "ani1x": torchani.models.ANI1x,
            "ani1ccx": torchani.models.ANI1ccx,
        }
        if parse_version(self.get_version()) >= parse_version("2.0"):
            ani_models["ani2x"] = torchani.models.ANI2x

        if name not in ani_models:
            raise InputError(f"TorchANI only accepts methods: {list(ani_models.keys())}")

        base = ani_models[name]()         # the actual TorchANI model
        wrapped = EnsembleEnergies(base)  # your compatibility wrapper
        self._CACHE[name] = wrapped

        return self._CACHE[name]

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs TorchANI in FF typing
        """
        # Check if exists and version
        self.found(raise_error=True)
        if parse_version(self.get_version()) < parse_version("0.9"):
            raise ResourceError("QCEngine's TorchANI wrapper requires version 0.9 or greater.")

        import numpy as np
        import torch
        import torchani

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # Failure flag
        ret_data = {"success": False}

        # helpers functions for version compatibility
        def species_tensor_from_symbols(symbols, model, device):
            """
            Returns a (1, A) LongTensor species tensor compatible with both:
              - torchani 2.7+ (atomic numbers input)
              - torchani 2.2.x (model-index input via model.species)
            """
            symbols = list(symbols)

            # Newer TorchANI (e.g., 2.7+) supports atomic-number conversion in torchani.utils
            to_Z = getattr(torchani.utils, "ChemicalSymbolsToAtomicNumbers", None)
            if to_Z is not None:
                return to_Z()(symbols).to(device).unsqueeze(0)

            # Older TorchANI (e.g., 2.2.4): convert to model indices using model.species
            from torchani.utils import ChemicalSymbolsToInts
            if not hasattr(model, "species"):
                raise RuntimeError(
                    "TorchANI version lacks ChemicalSymbolsToAtomicNumbers, and model has no .species; "
                    "cannot build species tensor."
                )
            return ChemicalSymbolsToInts(model.species)(symbols).to(device).unsqueeze(0)

        def get_ensemble_energies(model_result):
            """
            Return a 1D tensor of ensemble member energies for a single geometry.
            Works with:
              - TorchANI result objects: result.energies
              - old style tuples: (species, energies)
              - your current wrapper tuple output
            """
            if hasattr(model_result, "energies"):
                return model_result.energies
            if isinstance(model_result, (tuple, list)) and len(model_result) >= 2:
                return model_result[1]
            raise TypeError(f"Unrecognized TorchANI model output: {type(model_result)!r}")

        # Build model
        method = input_data.specification.model.method
        model = self.get_model(method).to(device)

        # Validate and build species tensor
        symbols = list(input_data.molecule.symbols)

        known_sym = {"H", "C", "N", "O"}
        if method.lower() == "ani2x":
            known_sym.update({"S", "F", "Cl"})
        unknown_sym = set(symbols) - known_sym
        if unknown_sym:
            raise InputError(f"TorchANI model '{method}' does not support symbols: {unknown_sym}.")

        num_atoms = len(symbols)
        species = species_tensor_from_symbols(symbols, model, device)

        # Build coord array -- single geometry, batch dim = 1
        geom_array = input_data.molecule.geometry.reshape(1, -1, 3) * ureg.conversion_factor("bohr", "angstrom")

        # Keep gradients enabled for gradient/hessian drivers
        coordinates = torch.tensor(geom_array.tolist(), device=device, requires_grad=True)

        # Run model
        result = model((species, coordinates))
        ensemble_energies = get_ensemble_energies(result)

        # Ensure 1D member energies for "no batching" usage.
        # Some paths might return shape (n_members, 1) or (n_members,) — normalize to (n_members,)
        if ensemble_energies.ndim == 2 and ensemble_energies.shape[-1] == 1:
            ensemble_energies = ensemble_energies[:, 0]

        # Scalar energy used for derivatives: mean across ensemble members
        energy = ensemble_energies.mean()

        ensemble_std = ensemble_energies.std(unbiased=False)
        ensemble_scaled_std = ensemble_std / np.sqrt(num_atoms)

        ret_data["properties"] = {"return_energy": energy.item()}

        if input_data.specification.driver == "energy":
            ret_data["return_result"] = ret_data["properties"]["return_energy"]
        elif input_data.specification.driver == "gradient":
            derivative = torch.autograd.grad(energy, coordinates, create_graph=False, retain_graph=False)[0].squeeze()
            ret_data["return_result"] = (
                np.asarray(derivative.cpu() * ureg.conversion_factor("angstrom", "bohr")).ravel().tolist()
            )
        elif input_data.specification.driver == "hessian":
            # torchani.utils.hessian expects coordinates and energies (scalar or (batch,))
            hessian = torchani.utils.hessian(coordinates, energies=energy)
            ret_data["return_result"] = np.asarray(hessian.cpu())
        else:
            raise InputError(
                f"TorchANI can only compute energy, gradient, and hessian driver methods. Found {input_data.specification.driver}."
            )

        #######################################################################
        # Description of the quantities stored in `extras`
        #
        # ensemble_energies:
        #   An energy array of all members (models) in an ensemble of models
        #
        # ensemble_energy_avg:
        #   The average value of energy array which is also recorded with as
        #   `energy` in QCEngine
        #
        # ensemble_energy_std:
        #   The standard deviation of energy array
        #
        # ensemble_per_root_atom_disagreement:
        #   The standard deviation scaled by the square root of N, with N being
        #   the number of atoms in the molecule. This is the quantity used in
        #   the query-by-committee (QBC) process in active learning to infer
        #   the reliability of the models in an ensemble, and produce more data
        #   points in the regions where this quantity is below a certain
        #   threshold (inclusion criteria)
        ret_data["input_data"] = input_data
        ret_data["molecule"] = input_data.molecule
        ret_data["extras"] = {
            # keep 1D array of ensemble member energies (no batching)
            "ensemble_energies": ensemble_energies.detach().cpu().numpy(),
            "ensemble_energy_avg": energy.item(),
            "ensemble_energy_std": ensemble_std.item(),
            "ensemble_per_root_atom_disagreement": ensemble_scaled_std.item(),
        }

        ret_data["provenance"] = Provenance(
            creator="torchani", version="unknown", routine="torchani.builtin.aev_computer"
        )

        ret_data["schema_name"] = "qcschema_atomic_result"
        ret_data["success"] = True

        return AtomicResult(**ret_data)
