"""
Calls the OpenMM executable.

Requires RDKit
"""
import json
import os
from typing import Dict

from qcelemental.models import Result

from .model import ProgramHarness
from ..exceptions import InputError, RandomError, ResourceError, UnknownError
from ..util import execute, popen, temporary_directory
from ..units import ureg


class OpenMMHarness(ProgramHarness):

    # TODO: verify correctness of these default params
    # where are these used? By QCFractal?
    _defaults = {
        "name": "OpenMM",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        # TODO: this harness requires RDKit as well, so this needs checking too
        return which_import('simtk.openmm',
                     return_bool=True,
                     raise_error=raise_error,
                     raise_msg='Please install via `conda install openmm -c omnia`.')

    @staticmethod
    def _process_molecule_rdkit(jmol):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Handle errors
        if abs(jmol.molecular_charge) > 1.e-6:
            raise InputError("RDKit does not currently support charged molecules.")

        if not jmol.connectivity:  # Check for empty list
            raise InputError("RDKit requires molecules to have a connectivity graph.")

        # Build out the base molecule
        base_mol = Chem.Mol()
        rw_mol = Chem.RWMol(base_mol)
        for sym in jmol.symbols:
            rw_mol.AddAtom(Chem.Atom(sym.title()))

        # Add in connectivity
        bond_types = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
        for atom1, atom2, bo in jmol.connectivity:
            rw_mol.AddBond(atom1, atom2, bond_types[bo])

        mol = rw_mol.GetMol()

        # Write out the conformer
        natom = len(jmol.symbols)
        conf = Chem.Conformer(natom)
        bohr2ang = ureg.conversion_factor("bohr", "angstrom")
        for line in range(natom):
            conf.SetAtomPosition(line, (bohr2ang * jmol.geometry[line, 0],
                                        bohr2ang * jmol.geometry[line, 1],
                                        bohr2ang * jmol.geometry[line, 2]))  # yapf: disable

        mol.AddConformer(conf)
        Chem.rdmolops.SanitizeMol(mol)

        return mol

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs OpenMM on given structure, inputs, in vacuum.
        """
        self.found(raise_error=True)

        from simtk.openmm.app import Simulation
        from simtk.openmm import units

        from openforcefield.topology import Molecule
        from openforcefield.typing.engines.smirnoff import ForceField
        
        ret_data["success"] = False

        # Set workdir to scratch
        # Location resolution order config.scratch_dir, /tmp
        parent = config.scratch_directory
        with temporary_directory(parent=parent, suffix="_openmm_scratch") as tmpdir:

            # Grab molecule, forcefield
            jmol = input_data.molecule
            
            # TODO: If urls are supported by
            # `openforcefield.typing.engines.smirnoff.ForceField` already, we
            # can eliminate the `offxml` and `url` distinction
            # URL processing can happen there
            if input_data.model.offxml:
                offxml = input_data.offxml
            elif input_data.model.url:
                # TODO: add pull of offxml data from a url
                pass
            else:
                raise InputError("OpenMM requires either `model.offxml` or `model.url` to be set")

            # Process molecule with RDKit
            rdkit_mol = self._process_molecule_rdkit(jmol)

            # Create an OFF Molecule from the RDKit Molecule
            off_mol = Molecule(rdkit_mol)

            # Load an OFF Force Field
            off_forcefield = ForceField(offxml)

            # Create OpenMM system in vacuum from forcefield, molecule
            off_top = off_mol.to_topology()
            openmm_top = off_top.to_openmm()
            openmm_system = off_forcefield.create_openmm_system(off_top)
            
            # Need an integrator for simulation even if we don't end up using it really
            integrator = LangevinIntegrator(300*units.kelvin, 1/units.picosecond, 0.002*units.picoseconds)

            # Set platform to CPU explicitly
            platform = openmm.Platform.getPlatformByName('CPU')

            # Initialize simulation
            simulation = Simulation(openmm_top, openmm_system, integrator, platform)

            # TODO: need a way to generate positions from the topology
            # likely an existing tool in the modules we're using since this is
            # just a molecule floating in space with no solvent
            simulation.context.setPositions(pdb.positions)

            # Execute driver
            try:
                if input_data.driver == "energy":
                    simulation.context.setPositions(pdb.positions)
                    state = simulation.context.getState(getEnergy=True)
                    ret_data["return_result"] = state.getPotentialEnergy()
                elif input_data.driver == "gradient":
                    # TODO: calculate forces with OpenMM
                    pass
                else:
                    raise InputError(f"OpenMM can only compute energy and gradient driver methods. Found {input_data.driver}.")
            except:
                raise
            else:
                success = True

            # GET OUTPUT, PUT INTO FORM WE NEED
            # FINAL STRUCTURE?
            # PROBABLY IN AN XML FILE
            if success:
                # PROCESS OUTPUT
                output_data = dict()

        # Dispatch errors
        if output_data["success"] is False:
            error_message = output_data["error"]["error_message"]

            if ("SIGSEV" in error_message) or ("SIGSEGV" in error_message) or ("segmentation fault" in error_message):
                raise RandomError(error_message)
            elif "TypeError: set_global_option" in error_message:
                raise InputError(error_message)
            else:
                raise UnknownError(error_message)

        ## FINALIZE OUTPUT DATA

        # Move several pieces up a level
        output_data["provenance"]["nthreads"] = output_data.pop("nthreads")

        return Result(**output_data)
