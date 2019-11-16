"""
Calls the OpenMM executable.

Requires RDKit
"""
import json
import os
import hashlib
from typing import Dict

from qcelemental.models import Provenance, Result
from qcelemental.util import which_import

from .model import ProgramHarness
from ..exceptions import InputError, RandomError, ResourceError, UnknownError
from ..util import execute, popen, temporary_directory
from ..units import ureg

from .rdkit import RDKitHarness


class OpenMMHarness(ProgramHarness):

    _CACHE = {}
    _CACHE_MAX_SIZE = 10

    _defaults = {
        "name": "OpenMM",
        "scratch": True,
        "thread_safe": True, # true if we use separate `openmm.Context` objects per thread
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }

    class Config(ProgramHarness.Config):
        pass

    def _get_off_forcefield(self, offxml):

        # perhaps want to avoid using a massive string as our key
        # if `offxml` is actual xml content?
        key = hashlib.sha256(offxml.encode()).hexdigest()

        off_forcefield = smirnoff.ForceField(offxml)

        return off_forcefield


    def _get_openmm_system(self, off_forcefield, off_top):

        #key = f"{mol_hash}-{input_model.method}-{hash(forcefield_schema)}"

        openmm_system = off_forcefield.create_openmm_system(off_top)

        return openmm_system

    def _cache_it(self, key, value):
        """Add to our LRU cache, possibly popping off least used key.

        """
        pass


    @staticmethod
    def found(raise_error: bool = False) -> bool:
        # this harness requires RDKit as well, so this needs checking too
        rdkit_found = RDKitHarness.found(raise_error=raise_error)

        openmm_found = which_import('.openmm',
                            return_bool=True,
                            raise_error=raise_error,
                            package='simtk',
                            raise_msg='Please install via `conda install openmm -c omnia`.')


        return (rdkit_found and openmm_found)

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs OpenMM on given structure, inputs, in vacuum.
        """
        self.found(raise_error=True)

        from simtk import openmm
        from simtk import unit

        import openforcefield.topology as offtop
        from openforcefield.typing.engines import smirnoff

        # temporary shim for making hartree available in `simtk.unit`
        if not hasattr(unit, 'hartree'):
            unit.hartree_base_unit = unit.ScaledUnit(4.3597447222071e-18, unit.joule, "hartree", "Ha")
            unit.hartree = unit.hartrees = unit.Unit({unit.hartree_base_unit: 1.0})
        
        # Failure flag
        ret_data = {"success": False}

        # get number of threads to use from `JobConfig.ncores`; otherwise, try environment variable
        nthreads = config.ncores
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
                offxml = input_data.model.offxml
            elif input_data.model.url:
                # TODO: add pull of offxml data from a url
                pass
            else:
                raise InputError("OpenMM requires either `model.offxml` or `model.url` to be set")

            # Process molecule with RDKit
            rdkit_mol = RDKitHarness._process_molecule_rdkit(jmol)

            # Create an Open Force Field `Molecule` from the RDKit Molecule
            off_mol = offtop.Molecule(rdkit_mol)

            # Load an Open Force Field `ForceField`
            off_forcefield = self._get_off_forcefield(offxml)

            # Create OpenMM system in vacuum from forcefield, molecule
            off_top = off_mol.to_topology()
            openmm_system = self._get_openmm_system(off_forcefield, off_top)
            
            # Need an integrator for simulation even if we don't end up using it really
            integrator = openmm.VerletIntegrator(1.0*unit.femtoseconds)

            # Set platform to CPU explicitly
            platform = openmm.Platform.getPlatformByName('CPU')

            # Initialize context
            context = openmm.Context(openmm_system, integrator, platform)

            # TODO: FIXME set number of threads;
            # this line currently throws an exception
            #if nthreads:
            #    platform.setPropertyValue(context, 'Threads', str(nthreads))

            # Set positions from our Open Force Field `Molecule`
            context.setPositions(off_mol.conformers[0])

            # Execute driver
            try:
                if input_data.driver == "energy":

                    # Compute the energy of the configuration
                    state = context.getState(getEnergy=True)

                    # Get the potential as a simtk.unit.Quantity, put into units of hartree
                    ret_data["return_result"] = state.getPotentialEnergy() / unit.hartree

                elif input_data.driver == "gradient":

                    # Get number of atoms
                    n_atoms = len(jmol.symbols)

                    # Compute the forces
                    state = context.getState(getForces=True)

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
