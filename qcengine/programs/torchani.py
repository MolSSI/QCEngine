"""
Calls the Psi4 executable.
"""

from qcelemental.models import ComputeError, FailedOperation, Provenance, Result
from qcelemental.util import parse_version, safe_version, which_import

from .executor import ProgramExecutor
from ..units import ureg


class TorchANIExecutor(ProgramExecutor):

    _CACHE = {}

    _defaults = {
        "name": "TorchANI",
        "scratch": False,
        "thread_safe": True,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": False
    }

    class Config(ProgramExecutor.Config):
        pass

    def __init__(self, **kwargs):
        super().__init__(**{**self._defaults, **kwargs})

    @staticmethod
    def found(raise_error: bool=False) -> bool:
        is_found = which_import("torchani", return_bool=True)

        if not is_found and raise_error:
            raise ModuleNotFoundError("Could not find 'torchani' in the Python path. Please 'pip install torchani'.")

        return is_found


    def get_version(self) -> str:
        self.found(raise_error=True)

        import torchani
        return safe_version(torchani.__version__)

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
            ret_data["error"] = ComputeError(
                error_type="version_error", error_message="QCEngine's TorchANI wrapper requires version 0.5 or greater.")
            return FailedOperation(input_data=input_data.dict(), **ret_data)

        import torch
        import numpy as np

        device = torch.device('cpu')

        # Failure flag
        ret_data = {"success": False}

        # Build model
        model = self.get_model(input_data.model.method)
        if model is False:
            ret_data["error"] = ComputeError(
                error_type="input_error", error_message="run_torchani only accepts the ANI1x or ANI1ccx method.")
            return FailedOperation(input_data=input_data.dict(), **ret_data)

        # Build species
        species = "".join(input_data.molecule.symbols)
        unknown_sym = set(species) - {"H", "C", "N", "O"}
        if unknown_sym:
            ret_data["error"] = ComputeError(
                error_type="input_error",
                error_message="The '{}' model does not support symbols: {}.".format(
                    input_data.model.method, unknown_sym))
            return FailedOperation(input_data=input_data.dict(), **ret_data)

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
            ret_data["error"] = ComputeError(
                error_type="input_error",
                error_message="run_torchani did not understand driver method '{}'.".format(input_data.driver))
            return FailedOperation(input_data=input_data.dict(), **ret_data)

        ret_data["provenance"] = Provenance(
            creator="torchani", version="unknown", routine='torchani.builtin.aev_computer')

        ret_data["schema_name"] = "qcschema_output"
        ret_data["success"] = True

        # Form up a dict first, then sent to BaseModel to avoid repeat kwargs which don't override each other
        return Result(**{**input_data.dict(), **ret_data})
