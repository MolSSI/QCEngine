"""
Calls the TorchANI package.
"""

from typing import Dict

from qcelemental.models import Provenance, Result
from qcelemental.util import parse_version, safe_version, which_import

from .model import ProgramHarness
from ..exceptions import InputError, ResourceError
from ..units import ureg


class TorchANIHarness(ProgramHarness):

    _CACHE = {}

    _defaults = {
        "name": "TorchANI",
        "scratch": False,
        "thread_safe": True,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": False
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool=False) -> bool:
        return which_import('torchani', return_bool=True, raise_error=raise_error, raise_msg='Please install via `pip install torchani`.')

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import('torchani')
        if which_prog not in self.version_cache:
            import torchani
            self.version_cache[which_prog] = safe_version(torchani.__version__)

        return self.version_cache[which_prog]

    def get_model(self, name: str) -> 'torchani.models.BuiltinModels':
        name = name.lower()

        if name in self._CACHE:
            return self._CACHE[name]

        import torch
        import torchani

        if name == "ani1x":
            self._CACHE[name] = torchani.models.ANI1x()

        elif name == "ani1ccx":
            self._CACHE[name] = torchani.models.ANI1ccx()

        else:
            return False

        return self._CACHE[name]

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs TorchANI in FF typing
        """

        # Check if existings and version
        self.found(raise_error=True)
        if parse_version(self.get_version()) < parse_version("0.5"):
            raise ResourceError("QCEngine's TorchANI wrapper requires version 0.5 or greater.")

        import torch
        import numpy as np

        device = torch.device('cpu')

        # Failure flag
        ret_data = {"success": False}

        # Build model
        model = self.get_model(input_data.model.method)
        if model is False:
            raise InputError("TorchANI only accepts the ANI1x or ANI1ccx method.")

        # Build species
        species = "".join(input_data.molecule.symbols)
        unknown_sym = set(species) - {"H", "C", "N", "O"}
        if unknown_sym:
            raise InputError(f"TorchANI model '{input_data.model.method}' does not support symbols: {unknown_sym}.")

        species = model.species_to_tensor(species).to(device).unsqueeze(0)

        # Build coord array
        geom_array = input_data.molecule.geometry.reshape(1, -1, 3) * ureg.conversion_factor("bohr", "angstrom")
        coordinates = torch.tensor(geom_array.tolist(), requires_grad=True, device=device)

        _, energy = model((species, coordinates))
        ret_data["properties"] = {"return_energy": energy.item()}

        if input_data.driver == "energy":
            ret_data["return_result"] = ret_data["properties"]["return_energy"]
        elif input_data.driver == "gradient":
            derivative = torch.autograd.grad(energy.sum(), coordinates)[0].squeeze()
            ret_data["return_result"] = np.asarray(
                derivative * ureg.conversion_factor("angstrom", "bohr")).ravel().tolist()
        else:
            raise InputError(f"TorchANI can only compute energy and gradient driver methods. Found {input_data.driver}.")

        ret_data["provenance"] = Provenance(
            creator="torchani", version="unknown", routine='torchani.builtin.aev_computer')

        ret_data["schema_name"] = "qcschema_output"
        ret_data["success"] = True

        # Form up a dict first, then sent to BaseModel to avoid repeat kwargs which don't override each other
        return Result(**{**input_data.dict(), **ret_data})
