"""
Calls the RDKit package.
"""

from typing import TYPE_CHECKING, Dict

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
        for atom1, atom2, bo in jmol.connectivity:
            rw_mol.AddBond(atom1, atom2, bond_types[bo])

        mol = rw_mol.GetMol()

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
        Chem.rdmolops.SanitizeMol(mol)

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

        def get_mol_descriptors(molecule):
            Chem.AssignStereochemistryFrom3D(molecule)
            descriptors = {
                # Hydrogen should be implicit in SMILES string, but not expected to be in database entries.
                "canonical_smiles": Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
            }

            return descriptors

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
            if driver == "energy":
                ret_data["return_result"] = ff.CalcEnergy() * ureg.conversion_factor("kJ / mol", "hartree")
                ret_data["properties"] = {"return_energy": ret_data["return_result"]}
            elif driver == "gradient":
                coef = ureg.conversion_factor("kJ / mol", "hartree") * ureg.conversion_factor("angstrom", "bohr")
                ret_data["return_result"] = [x * coef for x in ff.CalcGrad()]
                ret_data["properties"] = {"return_gradient": ret_data["return_result"]}
            else:
                pass
        elif driver == "hessian":
            raise InputError("RDKit does not support hessian calculation yet.")
        elif driver == "properties":
            if method == "descriptors":
                ret_data["properties"] = {}
                ret_data["return_result"] = get_mol_descriptors(mol)
                ret_data["provenance"] = Provenance(
                    creator="rdkit", version=rdkit.__version__, routine="get_molecular_descriptors"
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
