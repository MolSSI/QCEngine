"""
Calls the Psi4 executable.
"""

import time
import sys

from . import config
from . import units



def run_rdkit(json):
    """
    Runs RDKit in FF typing
    """

    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem 

    # Failure flag
    json["success"] = False

    # Build the Molecule
    jmol = json["molecule"]

    # Handle errors
    if ("molecular_charge" in jmol) and (abs(jmol["molecular_charge"]) < 1.e-6):
        json["error"] = "run_rdkit does not currently support charged molecules"
        return json

    if ("connectivity" not in jmol):
        json["error"] = "run_rdkit molecule must have a connectivity graph"
        return json

    # Build out the base molecule
    base_mol = Chem.Mol()
    rwMol = Chem.RWMol(base_mol)
    for sym in jmol["symbols"]:
        rwMol.AddAtom(Chem.Atom(sym))

    # Add in connectivity
    bond_types = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }
    for atom1, atom2, bo in jmol["connectivity"]:
        rwMol.AddBond(atom1, atom2, bond_types[bo])
        
    mol = rwMol.GetMol()

    # Write out the conformer
    natom = len(jmol["symbols"])
    conf = Chem.Conformer(natom)
    for line in range(natom):
        conf.SetAtomPosition(line, (units.bohr_to_angstrom * jmol["geometry"][line * 3],
                                    units.bohr_to_angstrom * jmol["geometry"][line * 3 + 1],
                                    units.bohr_to_angstrom * jmol["geometry"][line * 3 + 2]))

    mol.AddConformer(conf)
    Chem.rdmolops.SanitizeMol(mol)

    if json["model"]["method"] == "UFF":
        ff = AllChem.UFFGetMoleculeForceField(mol)
        all_params = AllChem.UFFHasAllMoleculeParams(mol)
    else:
        json["error"] = "run_rdkit can only accepts UFF methods"
        return json

    if all_params is False:
        json["error"] = "run_rdkit did not match all parameters to molecule"
        return json

    ff.Initialize()

    json["properties"] = {"return_energy": ff.CalcEnergy()}

    if json["driver"] == "energy":
        json["return_result"] = json["properties"]["return_energy"]
    elif json["driver"] == "gradient":
        json["return_result"] = [x / units.bohr_to_angstrom for x in ff.CalcGrad()]
    else:
        json["error"] = "run_rdkit did not understand driver method '{}'.".format(json["driver"])
        return json


    json["provenance"] = {
        "creator": "rdkit",
        "version": rdkit.__version__,
        "routine": "rdkit.Chem.AllChem.UFFGetMoleculeForceField"
        }

    json["schema_name"] = "qc_json_output"
    json["success"] = True
    

    return json
