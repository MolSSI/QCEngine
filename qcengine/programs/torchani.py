"""
Calls the Psi4 executable.
"""

from qcelemental.models import Result, ComputeError, Provenance, FailedOperation
from qcengine.units import ureg

_CACHE = {}


def get_model(name):
    name = name.lower()

    if name in _CACHE:
        return _CACHE[name]

    import torch
    import torchani

    if name == "ani1":
        # Build model
        builtin = torchani.neurochem.Builtins()
        model = torch.nn.Sequential(builtin.aev_computer, builtin.models, builtin.energy_shifter)
        _CACHE[name] = model

        return model

    else:
        return False


def torchani(input_data, config):
    """
    Runs TorchANI in FF typing
    """

    import numpy as np
    try:
        import torch
    except ImportError:
        raise ImportError("Could not find PyTorch in the Python path.")
    try:
        import torchani
    except ImportError:
        raise ImportError("Could not find TorchANI in the Python path.")

    device = torch.device('cpu')
    builtin = torchani.neurochem.Builtins()

    # Failure flag
    ret_data = {"success": False}

    # Build model
    model = get_model(input_data.model.method)
    if model is False:
        ret_data["error"] = ComputeError(
            error_type="input_error", error_message="run_torchani only accepts the ANI1 method.")
        return FailedOperation(input_data=input_data.dict(), **ret_data)

    # Build species
    species = "".join(input_data.molecule.symbols)
    unknown_sym = set(species) - {"H", "C", "N", "O"}
    if unknown_sym:
        ret_data["error"] = ComputeError(
            error_type="input_error",
            error_message="The '{}' model does not support symbols: {}.".format(input_data.model.method, unknown_sym))
        return FailedOperation(input_data=input_data.dict(), **ret_data)

    species = builtin.consts.species_to_tensor(species).to(device).unsqueeze(0)

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

    ret_data["provenance"] = Provenance(creator="torchani", version="unknown", routine='torchani.builtin.aev_computer')

    ret_data["schema_name"] = "qcschema_output"
    ret_data["success"] = True

    # Form up a dict first, then sent to BaseModel to avoid repeat kwargs which don't override each other
    return Result(**{**input_data.dict(), **ret_data})
