"""
Calls the Psi4 executable.
"""

from . import units



def run_rdkit(input_data):
    """
    Runs RDKit in FF typing
    """

    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Failure flag
    input_data["success"] = False

    # Build the Molecule
    jmol = input_data["molecule"]

    # Handle errors
    if ("molecular_charge" in jmol) and (abs(jmol["molecular_charge"]) < 1.e-6):
        input_data["error"] = "run_rdkit does not currently support charged molecules"
        return input_data

    if "connectivity" not in jmol:
        input_data["error"] = "run_rdkit molecule must have a connectivity graph"
        return input_data

    # Build out the base molecule
    base_mol = Chem.Mol()
    rw_mol = Chem.RWMol(base_mol)
    for sym in jmol["symbols"]:
        rw_mol.AddAtom(Chem.Atom(sym))

    # Add in connectivity
    bond_types = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE
    }
    for atom1, atom2, bo in jmol["connectivity"]:
        rw_mol.AddBond(atom1, atom2, bond_types[bo])

    mol = rw_mol.GetMol()

    # Write out the conformer
    natom = len(jmol["symbols"])
    conf = Chem.Conformer(natom)
    for line in range(natom):
        conf.SetAtomPosition(line, (units.bohr_to_angstrom * jmol["geometry"][line * 3],
                                    units.bohr_to_angstrom * jmol["geometry"][line * 3 + 1],
                                    units.bohr_to_angstrom * jmol["geometry"][line * 3 + 2]))

    mol.AddConformer(conf)
    Chem.rdmolops.SanitizeMol(mol)

    if input_data["model"]["method"] == "UFF":
        ff = AllChem.UFFGetMoleculeForceField(mol)
        all_params = AllChem.UFFHasAllMoleculeParams(mol)
    else:
        input_data["error"] = "run_rdkit can only accepts UFF methods"
        return input_data

    if all_params is False:
        input_data["error"] = "run_rdkit did not match all parameters to molecule"
        return input_data

    ff.Initialize()

    input_data["properties"] = {"return_energy": ff.CalcEnergy()}

    if input_data["driver"] == "energy":
        input_data["return_result"] = input_data["properties"]["return_energy"]
    elif input_data["driver"] == "gradient":
        input_data["return_result"] = [x / units.bohr_to_angstrom for x in ff.CalcGrad()]
    else:
        input_data["error"] = "run_rdkit did not understand driver method '{}'.".format(input_data["driver"])
        return input_data


    input_data["provenance"] = {
        "creator": "rdkit",
        "version": rdkit.__version__,
        "routine": "rdkit.Chem.AllChem.UFFGetMoleculeForceField"
        }

    input_data["schema_name"] = "qc_input_data_output"
    input_data["success"] = True


    return input_data
