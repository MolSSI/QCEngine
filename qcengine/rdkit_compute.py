"""
Calls the Psi4 executable.
"""

from . import units


def run_rdkit(ret_data):
    """
    Runs RDKit in FF typing
    """

    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Failure flag
    ret_data["success"] = False

    # Build the Molecule
    jmol = ret_data["molecule"]

    # Handle errors
    if ("molecular_charge" in jmol) and (abs(jmol["molecular_charge"]) < 1.e-6):
        ret_data["error_message"] = "run_rdkit does not currently support charged molecules"
        return ret_data

    if "connectivity" not in jmol:
        ret_data["error_message"] = "run_rdkit molecule must have a connectivity graph"
        return ret_data

    # Build out the base molecule
    base_mol = Chem.Mol()
    rw_mol = Chem.RWMol(base_mol)
    for sym in jmol["symbols"]:
        rw_mol.AddAtom(Chem.Atom(sym.title()))

    # Add in connectivity
    bond_types = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
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

    if ret_data["model"]["method"] == "UFF":
        ff = AllChem.UFFGetMoleculeForceField(mol)
        all_params = AllChem.UFFHasAllMoleculeParams(mol)
    else:
        ret_data["error_message"] = "run_rdkit can only accepts UFF methods"
        return ret_data

    if all_params is False:
        ret_data["error_message"] = "run_rdkit did not match all parameters to molecule"
        return ret_data

    ff.Initialize()

    ret_data["properties"] = {"return_energy": ff.CalcEnergy() / units.hartree_to_kj_mol}

    if ret_data["driver"] == "energy":
        ret_data["return_result"] = ret_data["properties"]["return_energy"]
    elif ret_data["driver"] == "gradient":
        coef = 1 / (units.bohr_to_angstrom * units.hartree_to_kj_mol)
        ret_data["return_result"] = [x * coef for x in ff.CalcGrad()]
    else:
        ret_data["error_message"] = "run_rdkit did not understand driver method '{}'.".format(ret_data["driver"])
        return ret_data

    ret_data["provenance"] = {
        "creator": "rdkit",
        "version": rdkit.__version__,
        "routine": "rdkit.Chem.AllChem.UFFGetMoleculeForceField"
    }

    ret_data["schema_name"] = "qc_ret_data_output"
    ret_data["success"] = True

    return ret_data
