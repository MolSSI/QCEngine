"""
Calls the RDKit package.
"""

from typing import TYPE_CHECKING, Dict
from importlib import import_module

from qcelemental.models import AtomicResult, Provenance
from qcelemental.util import safe_version, which_import

from ..exceptions import InputError
from ..units import ureg
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


class RDKitHarness(ProgramHarness):

    _defaults = {
        "name": "RDKit",
        "scratch": False,
        "thread_safe": True,
        "thread_parallel": False,
        "node_parallel": False,
        "managed_memory": False,
    }

    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def return_bondtype(bond):
        from rdkit import Chem
        # flexible bond typing would allow for calculated(float) bond orders
        # to be included in possible calculations (via NBO or something)
        # in the future while allowing informatics processes more comprehensive
        # molecular description.
        bond_types = {
            # see http://rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol
            # class rdkit.Chem.rdchem.BondType {values}
            0: Chem.BondType.UNSPECIFIED,
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            4: Chem.BondType.QUADRUPLE,
            5: Chem.BondType.QUINTUPLE,
            6: Chem.BondType.HEXTUPLE,
            7: Chem.BondType.ONEANDAHALF,
            8: Chem.BondType.TWOANDAHALF,
            9: Chem.BondType.THREEANDAHALF,
            10: Chem.BondType.FOURANDAHALF,
            11: Chem.BondType.FIVEANDAHALF,
            12: Chem.BondType.AROMATIC,
            13: Chem.BondType.IONIC,
            14: Chem.BondType.HYDROGEN,
            15: Chem.BondType.THREECENTER,
            16: Chem.BondType.DATIVEONE,
            17: Chem.BondType.DATIVE,
            18: Chem.BondType.DATIVEL,
            19: Chem.BondType.DATIVER,
            20: Chem.BondType.OTHER,
            21: Chem.BondType.ZERO
        }
        if float(bond).is_integer():
            bondtype = bond_types[bond]
        else:
            half_orders = {
                0: 21, 0.5: 0, 1: 1, 1.5: 7, 2: 2, 2.5: 8,
                3: 3, 3.5: 9, 4: 4, 4.5: 10, 5: 5
            }
            nearest_half = round(bond * 2) / 2
            bondtype = bond_types[half_orders[nearest_half]]

        return bondtype

    @staticmethod
    def _process_molecule_rdkit(jmol):
        from rdkit import Chem

        # Handle errors
        if abs(jmol.molecular_charge) > 1.0e-6:
            raise InputError("RDKit does not currently support charged molecules.")

        if not jmol.connectivity:  # Check for empty list
            raise InputError("RDKit requires molecules to have a connectivity graph.")

        # Build out the base molecule
        base_mol = Chem.Mol()
        rw_mol = Chem.RWMol(base_mol)
        for sym in jmol.symbols:
            rw_mol.AddAtom(Chem.Atom(sym.title()))

        # Add in connectivity
        for atom1, atom2, bo in jmol.connectivity:
            bondtype = RDKitHarness.return_bondtype(bo)
            rw_mol.AddBond(atom1, atom2, bondtype)

        mol = rw_mol.GetMol()

        if jmol.name is not None:
            mol.SetProp("_Name", jmol.name)

        # Write out the conformer
        natom = len(jmol.symbols)
        conf = Chem.Conformer(natom)
        bohr2ang = ureg.conversion_factor("bohr", "angstrom")
        for line in range(natom):
            conf.SetAtomPosition(
                line,
                (
                    bohr2ang * jmol.geometry[line, 0],
                    bohr2ang * jmol.geometry[line, 1],
                    bohr2ang * jmol.geometry[line, 2],
                ),
            )

        mol.AddConformer(conf)

        # set sanitization flags to preserve aromaticity in heterocycles
        s_flags = Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP|Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY|Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES|Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS|Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
        # sanitize molecule
        Chem.rdmolops.SanitizeMol(mol, sanitizeOps=s_flags)

        return mol

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which_import(
            "rdkit",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install rdkit -c conda-forge`.",
        )

    def get_version(self) -> str:
        """Return the currently used version of RDKit."""
        self.found(raise_error=True)

        which_prog = which_import("rdkit")
        if which_prog not in self.version_cache:
            import rdkit

            self.version_cache[which_prog] = safe_version(rdkit.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs RDKit
        """

        self.found(raise_error=True)
        import rdkit
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Failure flag
        ret_data = {"success": False}

        # Build the Molecule
        jmol = input_data.molecule
        mol = self._process_molecule_rdkit(jmol)
        method = input_data.model.method.lower()
        driver = input_data.driver

        # CI tests fail without this property set
        ret_data["properties"] = {"calcinfo_natom": len(jmol.symbols)}

        if driver in ["energy", "gradient"]:
            if method == "uff":
                ff = AllChem.UFFGetMoleculeForceField(mol)
                all_params = AllChem.UFFHasAllMoleculeParams(mol)
                ret_data["provenance"] = Provenance(
                    creator="rdkit", version=rdkit.__version__, routine="rdkit.Chem.AllChem.UFFGetMoleculeForceField"
                )
            elif method in ["mmff94", "mmff94s"]:
                props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=input_data.model.method)
                ff = AllChem.MMFFGetMoleculeForceField(mol, props)
                all_params = AllChem.MMFFHasAllMoleculeParams(mol)
                ret_data["provenance"] = Provenance(
                    creator="rdkit", version=rdkit.__version__, routine="rdkit.Chem.AllChem.MMFFGetMoleculeForceField"
                )
            else:
                raise InputError(
                    "RDKit only supports the UFF, MMFF94, and MMFF94s methods for energy/gradient calculations currently."
                )
            if all_params is False:
                raise InputError("RDKit parameters not found for all atom types in molecule.")
            ff.Initialize()
            ret_data["properties"]["return_energy"] = ff.CalcEnergy() * ureg.conversion_factor("kJ / mol", "hartree")
            if driver == "energy":
                ret_data["return_result"] = ret_data["properties"]["return_energy"]
            elif driver == "gradient":
                coef = ureg.conversion_factor("kJ / mol", "hartree") * ureg.conversion_factor("angstrom", "bohr")
                ret_data["properties"]["return_gradient"] = [x * coef for x in ff.CalcGrad()]
                ret_data["return_result"] = ret_data["properties"]["return_gradient"]
            else:
                pass
        elif driver == "hessian":
            raise InputError("RDKit does not support hessian calculation yet.")
        elif driver == "properties":
            if method == "descriptors":
                get_descriptors = DescriptorConstructor(mol)
                descriptors = get_descriptors.construct_dict()
                ret_data["properties"]["descriptors"] = descriptors
                ret_data["return_result"] = descriptors
                ret_data["provenance"] = Provenance(
                    creator="rdkit", version=rdkit.__version__, routine="compute_descriptors"
                )
            else:
                raise InputError(
                    f'QCEngine has implemented the following methods for the "properties" driver in RDKit:\n\t"descriptors"'
                )
        else:
            raise InputError(f"RDKit does not support the {driver} driver.")

        ret_data["schema_name"] = "qcschema_output"
        ret_data["success"] = True

        # Form up a dict first, then sent to BaseModel to avoid repeat kwargs which don't override each other
        return AtomicResult(**{**input_data.dict(), **ret_data})

class DescriptorConstructor:

    def __init__(self, molecule):
        self.Chem = import_module('rdkit.Chem')
        self.Descriptors = import_module('rdkit.Chem.Descriptors')
        self.Descriptors3D = import_module('rdkit.Chem.Descriptors3D')
        self.molecule = self.Chem.Mol(molecule)
        self.process_mol()

    def process_mol(self):
        Chem = self.Chem
        Chem.AssignStereochemistryFrom3D(self.molecule)

        # Create RDKit(C++) parameters object
        parameters = Chem.AdjustQueryParameters()

        # Set parameters
        parameters.adjustConjugatedFiveRings = True
        parameters.aromatizeIfPossible = True
        parameters.useStereoCareForBonds = True
        parameters.makeAtomsGeneric = False

        # process the molecule according to query parameters
        mol = Chem.AdjustQueryProperties(self.molecule, params=parameters)

        self.molecule = Chem.Mol(mol)

    @property
    def stereochemistry(self):
        Chem = self.Chem
        mol = self.molecule

        def get_doublebond_stereo(mol):
            stereobonds = []
            for bond in mol.GetBonds():  # loop through the bonds
                is_stereo = False
                stereo = bond.GetStereo()  # get stereochemistry information from them
                # assign a label based on the BondStereo type
                if stereo == Chem.BondStereo.STEREOZ or stereo == Chem.BondStereo.STEREOCIS:
                    is_stereo = True
                    s_label = 'Z'
                elif stereo == Chem.BondStereo.STEREOE or stereo == Chem.BondStereo.STEREOTRANS:
                    is_stereo = True
                    s_label = 'E'
                else:  # skip over bonds without e/z stereochemistry
                    is_stereo = False
                if is_stereo:
                    # get the atom indices
                    idx1 = bond.GetBeginAtomIdx()
                    idx2 = bond.GetEndAtomIdx()
                    stereobonds += [(idx1, idx2, s_label)]  # add the bonds to the list
            if stereobonds != []:  # set a new property if not empty
                mol.SetProp("double-bond stereo", str(stereobonds))
            else:
                stereobonds = None

            return stereobonds

        db_stereo = get_doublebond_stereo(mol)
        if db_stereo is not None:
            stereochem = Chem.FindMolChiralCenters(mol) + get_doublebond_stereo(mol)
        else:
            stereochem = Chem.FindMolChiralCenters(mol)

        if stereochem == []:
            stereochem = None

        return stereochem

    @property
    def canonical_smiles(self):
        Chem = self.Chem
        mol = Chem.Mol(self.molecule)    # make a new copy so that RemoveHs won't edit the original in place.

        # Hydrogen should be implicit in SMILES string, but not expected to be in database entries.
        smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)

        return smiles

    @property
    def aromaticity(self):
        Chem = self.Chem
        mol = self.molecule
        arom_bonds = []
        arom_atoms = []

        # Get the bonds and find aromatic ones
        bonds = mol.GetBonds()
        for bond in bonds:
            if bond.GetIsAromatic():
                arom_bonds += [bond.GetIdx()]  # add bond index to list
                # Adjusting properties doesn't set atom aromaticity, just bonds, so set atoms too.
                atom1 = bond.GetBeginAtom()
                atom1.SetIsAromatic(True)
                atom2 = bond.GetEndAtom()
                atom2.SetIsAromatic(True)

        for atom in mol.GetAromaticAtoms():
            arom_atoms += [atom.GetIdx()]

        # If not empty, make a dictionary so aromatic atoms or bonds can be returned later
        if arom_bonds == []:
            indices = None    # for readability
        else:
            indices = {"bond_indices": arom_bonds, "atom_indices": arom_atoms}

        # Create property for aromaticity info
        mol.SetProp("aromaticity", str(indices))

        return indices

    @property
    def ring_info(self):
        Chem = self.Chem
        mol = self.molecule
        ring_info = {}

        def isRingAromatic(mol, ring):
            for idx in ring:
                if not mol.GetBondWithIdx(idx).GetIsAromatic():
                    return False

            return True

        def isRingHeterocyclic(mol, atom_ring):
            isHeterocycle = False

            for idx in atom_ring:
                atom = mol.GetAtomWithIdx(idx)
                symbol = atom.GetSymbol()
                if symbol == 'C' or isHeterocycle:
                    continue
                else:
                    isHeterocycle = True

            return isHeterocycle

        ri = mol.GetRingInfo()
        atom_rings = ri.AtomRings()
        bond_rings = ri.BondRings()

        nrings = len(atom_rings)

        if nrings == 0:
            ring_info = None
        else:
            ring_info['num_rings'] = nrings
            ring_info['rings'] = {}
            for index, ring in enumerate(bond_rings):
                ring_info['rings'][index] = {}
                ring_info['rings'][index]['size'] = len(atom_rings[index])
                ring_info['rings'][index]['aromatic'] = isRingAromatic(mol, ring)
                ring_info['rings'][index]['heterocyclic'] = isRingHeterocyclic(mol, atom_rings[index])
                ring_info['rings'][index]['atoms'] = list(atom_rings[index])

        return ring_info

    @property
    def gasteiger_charges(self):
        Chem = self.Chem
        mol = self.molecule

        mol.ComputeGasteigerCharges(nIter=48, throwOnParamFailure=True)
        g_charges = []

        # Charges are a property of the atoms in the molecule, so they must be iterated through to retrieve
        for atom in mol.GetAtoms():
            g_charges += [float(atom.GetProp("_GasteigerCharge"))]

        return g_charges

    def construct_dict(self):
        descriptors = {
            "canonical_smiles": self.canonical_smiles,
            "stereochemistry": self.stereochemistry,
            "aromaticity": self.aromaticity,
            "ring_info": self.ring_info,
            "gasteiger_charges": self.gasteiger_charges,
        }

        return descriptors

