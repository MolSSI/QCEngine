"""
Calls the Psi4 executable.
"""

from qcengine import units


def run_torchani(ret_data):
    """
    Runs TorchANI in FF typing
    """

    import numpy as np
    import torch
    import torchani

    device = torch.device('cpu')

    # Failure flag
    ret_data["success"] = False

    # Build model and species array
    if ret_data["model"]["method"] == "ANI1":

        # Build model
        builtin = torchani.neurochem.Builtins()
        model = torch.nn.Sequential(builtin.aev_computer, builtin.models, builtin.energy_shifter)

        # Build species
        unknown_sym = set(ret_data["molecule"]["symbols"]) - {"H", "C", "N", "O"}
        if unknown_sym:
            raise KeyError(
                "The '{}' model does not support symbols: {}.".format(ret_data["model"]["method"], unknown_sym))

        species = builtin.consts.species_to_tensor("".join(ret_data["molecule"]["symbols"])).to(device).unsqueeze(0)
        species = species.double()
    else:
        ret_data["error_message"] = "run_torchani only accepts the ANI1 method."
        return ret_data

    # Build coord array
    geom_array = np.array(ret_data["molecule"]["geometry"]).reshape(1, -1, 3) * units.bohr_to_angstrom
    coordinates = torch.tensor(geom_array, requires_grad=True, device=device)
    coordinates = coordinates.double()

    _, energy = model((species, coordinates))
    ret_data["properties"] = {"return_energy": energy.item()}

    if ret_data["driver"] == "energy":
        ret_data["return_result"] = ret_data["properties"]["return_energy"]
    elif ret_data["driver"] == "gradient":
        derivative = torch.autograd.grad(energy.sum(), coordinates)[0].squeeze()
        ret_data["return_result"] = (derivative / units.bohr_to_angstrom).ravel().tolist()
    else:
        ret_data["error_message"] = "run_torchani did not understand driver method '{}'.".format(ret_data["driver"])
        return ret_data

    ret_data["provenance"] = {"creator": "torchani", "version": "unknown", "routine": "torchani.builtin.aev_computer"}

    ret_data["schema_name"] = "qc_schema_output"
    ret_data["success"] = True

    return ret_data
