import os
from typing import TYPE_CHECKING, Any, ClassVar, Dict

from qcelemental.models.v2 import AtomicResult, Provenance
from qcelemental.util import safe_version, which_import

from qcengine.exceptions import InputError, ResourceError
from qcengine.programs.model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models.v2 import AtomicInput

    from qcengine.config import TaskConfig


class FairChemHarness(ProgramHarness):
    """A harness to run Meta's UMA models through fairchem <https://github.com/facebookresearch/fairchem>

    ``model.method`` is a fairchem checkpoint name (``uma-s-1p2p1``, ``uma-m-1p1``, ...), a prefix of one
    (``uma`` resolves to the newest small model), or the path to a ``.pt`` file. Named checkpoints must already
    be in fairchem's cache (``~/.cache/fairchem`` or ``FAIRCHEM_CACHE_DIR``); the harness never downloads.

    Keywords:
        task_name (str, required): UMA head: ``omol``, ``omat``, ``oc20``, ``odac`` or ``omc``.
        device (str): ``auto`` (default), ``cpu`` or ``cuda``.
        inference_settings (str): ``turbo`` (default), ``default`` or ``batch``. The first two compile for one
            fixed composition/charge/spin/task and fall back on change; use ``batch`` for varied systems.
    """

    _CACHE = {}

    _defaults: ClassVar[Dict[str, Any]] = {
        "name": "FairChem",
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
            "fairchem.core",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `pip install fairchem-core`.",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("fairchem.core")
        if which_prog not in self.version_cache:
            import fairchem.core

            self.version_cache[which_prog] = safe_version(fairchem.core.__version__)

        return self.version_cache[which_prog]

    @staticmethod
    def resolve_device(keywords: Dict[str, Any]) -> str:
        """Honor an explicit ``device`` keyword, else cuda if available. fairchem accepts only cpu and cuda."""
        import torch

        requested = keywords.get("device", "auto")
        if requested != "auto":
            return requested
        return "cuda" if torch.cuda.is_available() else "cpu"

    @staticmethod
    def resolve_model(method: str) -> str:
        """Map ``method`` to a checkpoint: a file path, an exact name, or the newest name with that prefix."""
        if os.path.isfile(method):
            return method

        from fairchem.core import pretrained_mlip

        name = method.lower()
        if name in pretrained_mlip.available_models:
            return name
        matches = [m for m in pretrained_mlip.available_models if m.startswith(name)]
        if not matches:
            raise InputError(
                f"FairChem model {method} not recognized. Available: {sorted(pretrained_mlip.available_models)}"
            )
        return max(matches)

    @staticmethod
    def check_weights(name: str) -> None:
        """Raise unless the checkpoint is already in fairchem's Hugging Face cache."""
        from fairchem.core import pretrained_mlip
        from fairchem.core._config import CACHE_DIR
        from huggingface_hub import try_to_load_from_cache

        ckpt = pretrained_mlip._MODEL_CKPTS.checkpoints[name]
        path = "/".join(filter(None, [ckpt.subfolder, ckpt.filename]))
        if not isinstance(try_to_load_from_cache(ckpt.repo_id, path, cache_dir=CACHE_DIR, revision=ckpt.revision), str):
            raise ResourceError(
                f"FairChem weights for {name} are not cached in {CACHE_DIR}. They are gated: accept the license at "
                f"https://huggingface.co/{ckpt.repo_id}, run `hf auth login`, then pull them once with\n"
                f'  python -c "from fairchem.core import pretrained_mlip; '
                f"pretrained_mlip.get_predict_unit('{name}', device='cpu')\""
            )

    def load_model(self, name: str, device: str, inference_settings: str):
        key = (name, device, inference_settings)
        if key in self._CACHE:
            return self._CACHE[key]

        from fairchem.core import pretrained_mlip
        from fairchem.core.units.mlip_unit import load_predict_unit

        if os.path.isfile(name):
            self._CACHE[key] = load_predict_unit(name, inference_settings, device=device)
            return self._CACHE[key]

        self.check_weights(name)
        try:
            self._CACHE[key] = pretrained_mlip.get_predict_unit(
                name, inference_settings=inference_settings, device=device
            )
        except OSError as e:  # reference-energy files missing from the cache while offline
            raise ResourceError(f"Could not load {name} ({type(e).__name__}: {e}).") from e
        return self._CACHE[key]

    def compute(self, input_data: "AtomicInput", config: "TaskConfig"):
        self.found(raise_error=True)

        from qcengine.units import ureg

        spec = input_data.specification
        if spec.driver not in ["energy", "gradient"]:
            raise InputError(
                f"FairChem can only compute energy and gradient driver methods. Requested {spec.driver} not supported."
            )
        task_name = spec.keywords.get("task_name")
        if task_name is None:
            raise InputError("FairChem requires the `task_name` keyword: one of omol, omat, oc20, odac, omc.")

        name = self.resolve_model(spec.model.method)
        device = self.resolve_device(spec.keywords)
        inference_settings = spec.keywords.get("inference_settings", "turbo")
        predictor = self.load_model(name, device, inference_settings)

        from ase import Atoms
        from fairchem.core import FAIRChemCalculator

        try:
            calc = FAIRChemCalculator(predictor, task_name=task_name)
        except ValueError as e:  # unknown task_name
            raise InputError(str(e)) from e

        # the omol head reads integer charge and spin multiplicity from atoms.info
        atoms = Atoms(
            numbers=input_data.molecule.atomic_numbers,
            positions=input_data.molecule.geometry * ureg.conversion_factor("bohr", "angstrom"),
            info={
                "charge": int(input_data.molecule.molecular_charge),
                "spin": int(input_data.molecule.molecular_multiplicity),
            },
        )
        atoms.calc = calc
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()

        ret_data = {
            "input_data": input_data,
            "molecule": input_data.molecule,
            "success": False,
            "properties": {
                "return_energy": energy * ureg.conversion_factor("eV", "hartree"),
                "return_gradient": -1.0 * forces * ureg.conversion_factor("eV / angstrom", "hartree / bohr"),
                "calcinfo_natom": len(input_data.molecule.atomic_numbers),
            },
            "extras": {
                "fairchem": {
                    "model": name,
                    "task_name": task_name,
                    "device": device,
                    "inference_settings": inference_settings,
                }
            },
        }
        if spec.driver == "energy":
            ret_data["return_result"] = ret_data["properties"]["return_energy"]
        else:
            ret_data["return_result"] = ret_data["properties"]["return_gradient"]

        ret_data["provenance"] = Provenance(
            creator="fairchem", version=self.get_version(), routine="FAIRChemCalculator"
        )

        ret_data["success"] = True

        return AtomicResult(**ret_data)
