from typing import TYPE_CHECKING, Any, ClassVar, Dict

from qcelemental.models.v2 import AtomicResult, Provenance
from qcelemental.util import safe_version, which_import

from qcengine.exceptions import InputError
from qcengine.programs.model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models.v2 import AtomicInput

    from qcengine.config import TaskConfig


class ORBHarness(ProgramHarness):
    """A harness to run ORB models <https://github.com/orbital-materials/orb-models>"""

    _CACHE = {}

    _defaults: ClassVar[Dict[str, Any]] = {
        "name": "ORB",
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
            "orb_models",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `pip install orb-models`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("orb_models")
        if which_prog not in self.version_cache:
            import orb_models

            self.version_cache[which_prog] = safe_version(orb_models.__version__)

        return self.version_cache[which_prog]

    @staticmethod
    def resolve_device(keywords: Dict[str, Any]) -> str:
        """Pick the torch device, honouring an explicit ``device`` keyword over autodetection.

        orb loads models with an unconditional ``.cuda(device)``, so cpu and cuda are the only
        options (no mps)
        """
        import torch

        requested = keywords.get("device", "auto")
        if requested != "auto":
            return requested
        return "cuda" if torch.cuda.is_available() else "cpu"

    def load_model(self, name: str, device: str):
        """Return the (model, adapter) pair orb's loaders hand back."""
        # cache on device too, else a cpu-loaded model is handed back for a gpu request
        key = (name.lower(), device)
        if key in self._CACHE:
            return self._CACHE[key]

        from orb_models.forcefield import pretrained

        if key[0] not in pretrained.ORB_PRETRAINED_MODELS:
            raise InputError(f"ORB model {name} not recognized. Available: {sorted(pretrained.ORB_PRETRAINED_MODELS)}")

        self._CACHE[key] = pretrained.ORB_PRETRAINED_MODELS[key[0]](device=device)
        return self._CACHE[key]

    def compute(self, input_data: "AtomicInput", config: "TaskConfig"):
        self.found(raise_error=True)

        from qcengine.units import ureg

        method = input_data.specification.model.method
        device = self.resolve_device(input_data.specification.keywords)
        model, adapter = self.load_model(name=method, device=device)

        if input_data.specification.driver not in ["energy", "gradient"]:
            raise InputError(
                f"ORB can only compute energy and gradient driver methods. Requested {input_data.specification.driver} not supported."
            )

        from ase import Atoms

        # the omol models are charge/spin conditioned and require both together in atoms.info
        atoms = Atoms(
            numbers=input_data.molecule.atomic_numbers,
            positions=input_data.molecule.geometry * ureg.conversion_factor("bohr", "angstrom"),
            info={
                "charge": input_data.molecule.molecular_charge,
                "spin": input_data.molecule.molecular_multiplicity,
            },
        )
        batch = adapter.batch([adapter.from_ase_atoms(atoms, device=device)])
        out = model.predict(batch)

        # conservative regressors name their outputs on the model (and may suffix a level of
        # theory); direct regressors just key on the head names
        energy = out[getattr(model, "energy_name", "energy")].item()
        forces = out[getattr(model, "grad_forces_name", "forces")].detach().cpu().numpy()

        ret_data = {
            "input_data": input_data,
            "molecule": input_data.molecule,
            "success": False,
            "properties": {
                "return_energy": energy * ureg.conversion_factor("eV", "hartree"),
                "return_gradient": -1.0 * forces * ureg.conversion_factor("eV / angstrom", "hartree / bohr"),
                "calcinfo_natom": len(input_data.molecule.atomic_numbers),
            },
            "extras": {"orb": {"device": device}},
        }
        if input_data.specification.driver == "energy":
            ret_data["return_result"] = ret_data["properties"]["return_energy"]
        else:
            ret_data["return_result"] = ret_data["properties"]["return_gradient"]

        ret_data["provenance"] = Provenance(creator="orb_models", version=self.get_version(), routine="load_model")

        ret_data["success"] = True

        return AtomicResult(**ret_data)
