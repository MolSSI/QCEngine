"""
Calls the OpenMM executable.

Requires RDKit
"""
import json
import os
from typing import Dict

from qcelemental.models import Provenance, Result

from .model import ProgramHarness
from ..exceptions import InputError, RandomError, ResourceError, UnknownError
from ..util import execute, popen, temporary_directory
from ..units import ureg

from .rdkit import RDKitHarness


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
        # this harness requires RDKit as well, so this needs checking too
        rdkit_found = RDKitHarness.found(raise_error=True)

        openmm_found = which_import('simtk.openmm',
                            return_bool=True,
                            raise_error=raise_error,
                            raise_msg='Please install via `conda install openmm -c omnia`.')


        return (rdkit_found and openmm_found)

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

        from simtk import openmm
        from simtk import unit

        from openforcefield.topology import Molecule
        from openforcefield.typing.engines.smirnoff import ForceField
        
        ret_data["success"] = False

        # get number of threads to use from `JobConfig; otherwise, try environment variable
        nthreads = config.nthreads
        if nthreads is None:
            nthreads = os.environ.get('OPENMM_CPU_THREADS')

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
            openmm_system = off_forcefield.create_openmm_system(off_top)
            
            # Need an integrator for simulation even if we don't end up using it really
            integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)

            # Set platform to CPU explicitly
            platform = openmm.Platform.getPlatformByName('CPU')

            # Initialize context
            context = openmm.Context(openmm_system, integrator, platform)

            # set number of threads
            if nthreads:
                platform.setPropertyValue(context, 'Threads', nthreads)

            # Set positions from our Open Force Field `Molecule`
            context.setPositions(off_mol.conformers[0])

            # Execute driver
            try:
                if input_data.driver == "energy":

                    # Compute the energy of the configuration
                    state = simulation.context.getState(getEnergy=True)

                    # Get the potential as a simtk.unit.Quantity, put into units of hartree
                    ret_data["return_result"] = state.getPotentialEnergy() / unit.hartree

                elif input_data.driver == "gradient":

                    # Get number of atoms
                    n_atoms = len(jmol.symbols)

                    # Compute the forces
                    state = simulation.context.getState(getForces=True)

                    # Get the gradient as a simtk.unit.Quantity with shape (n_atoms, 3)
                    gradient = state.getForces(asNumpy=True)

                    # Convert to hartree/bohr and reformat as 1D array
                    ret_data["return_result"] = (gradient / (unit.hartree/unit.bohr)).reshape([n_atoms * 3])

                else:
                    raise InputError(f"OpenMM can only compute energy and gradient driver methods. Found {input_data.driver}.")
            except:
                raise
            else:
                ret_data["success"] = True

        # Move several pieces up a level
        ret_data["provenance"] = Provenance(creator="openmm",
                                            version=openmm.__version__,
                                            nthreads=nthreads)

        return Result(**{**input_data.dict(), **ret_data})
