"""
Calls the Psi4 executable.
"""

from qcengine import units

_CACHE = {}


def get_model(name):
    if name in _CACHE:
        print()
        print("CACHE")
        return _CACHE[name]

    import torch
    import torchani

    if name.lower() == "ani1":
        # Build model
        builtin = torchani.neurochem.Builtins()
        model = torch.nn.Sequential(builtin.aev_computer, builtin.models, builtin.energy_shifter)

        return model

    else:
        return False


def run_torchani(ret_data):
    """
    Runs TorchANI in FF typing
    """

    import numpy as np
    import torch
    import torchani

    device = torch.device('cpu')
    builtin = torchani.neurochem.Builtins()

    # Failure flag
    ret_data["success"] = False

    # Build model
    model = get_model(ret_data["model"]["method"])
    if model is False:
        ret_data["error_message"] = "run_torchani only accepts the ANI1 method."
        return ret_data

    # Build species
    species = "".join(ret_data["molecule"]["symbols"])
    unknown_sym = set(species) - {"H", "C", "N", "O"}
    if unknown_sym:
        ret_data["error_message"] = "The '{}' model does not support symbols: {}.".format(
            ret_data["model"]["method"], unknown_sym)
        return ret_data

    species = builtin.consts.species_to_tensor(species).to(device).unsqueeze(0)

    # Build coord array
    geom_array = np.array(ret_data["molecule"]["geometry"]).reshape(1, -1, 3) * units.bohr_to_angstrom
    coordinates = torch.tensor(geom_array.tolist(), requires_grad=True, device=device)

    _, energy = model((species, coordinates))
    ret_data["properties"] = {"return_energy": energy.item()}

    if ret_data["driver"] == "energy":
        ret_data["return_result"] = ret_data["properties"]["return_energy"]
    elif ret_data["driver"] == "gradient":
        derivative = torch.autograd.grad(energy.sum(), coordinates)[0].squeeze()
        ret_data["return_result"] = np.asarray(derivative / units.bohr_to_angstrom).ravel().tolist()
    else:
        ret_data["error_message"] = "run_torchani did not understand driver method '{}'.".format(ret_data["driver"])
        return ret_data

    ret_data["provenance"] = {"creator": "torchani", "version": "unknown", "routine": "torchani.builtin.aev_computer"}

    ret_data["schema_name"] = "qc_schema_output"
    ret_data["success"] = True

    return ret_data
