"""
Calls the OpenMM executable.

Requires RDKit
"""
import datetime
import hashlib
import os

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..config import TaskConfig
    from qcelemental.models import AtomicInput

from qcelemental.models import AtomicResult, Provenance
from qcelemental.util import which_import

from ..exceptions import InputError
from ..util import capture_stdout
from .model import ProgramHarness
from .rdkit import RDKitHarness


class OpenMMHarness(ProgramHarness):

    _CACHE = {}
    _CACHE_MAX_SIZE = 10

    _defaults = {
        "name": "OpenMM",
        "scratch": True,
        "thread_safe": True,  # true if we use separate `openmm.Context` objects per thread
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }

    class Config(ProgramHarness.Config):
        pass

    def _get_off_forcefield(self, hashstring, offxml):

        from openforcefield.typing.engines import smirnoff

        key = hashlib.sha256(hashstring.encode()).hexdigest()

        # get forcefield from cache, build new one if not present
        off_forcefield = self._get_cache(key) if key in self._CACHE else smirnoff.ForceField(offxml)

        # cache forcefield, no matter what
        # handles updating time touched, dropping items if cache too large
        self._cache_it(key, off_forcefield)

        return off_forcefield

    def _get_openmm_system(self, off_forcefield, off_top):

        # key = f"{mol_hash}-{input_model.method}-{hash(forcefield_schema)}"

        openmm_system = off_forcefield.create_openmm_system(off_top)

        return openmm_system

    def _get_cache(self, key):
        return self._CACHE[key]["value"]

    def _cache_it(self, key, value):
        """Add to our LRU cache, possibly popping off least used key.

        """
        self._CACHE[key] = {"value": value, "last_used": datetime.datetime.utcnow()}

        # if cache is beyond max size, whittle it down by dropping entry least
        # recently used
        while len(self._CACHE) > self._CACHE_MAX_SIZE:
            for i, (key, value) in enumerate(self._CACHE.items()):
                if i == 0:
                    oldest = value["last_used"]
                    oldest_key = key
                else:
                    oldest, oldest_key = (
                        (value["last_used"], key) if (value["last_used"] < oldest) else oldest,
                        oldest_key,
                    )

            self._CACHE.pop(oldest_key)

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        # this harness requires RDKit as well, so this needs checking too
        rdkit_found = RDKitHarness.found(raise_error=raise_error)

        openff_found = which_import(
            "openforcefield",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install openforcefield -c omnia`.",
        )

        openmm_found = which_import(
            ".openmm",
            return_bool=True,
            raise_error=raise_error,
            package="simtk",
            raise_msg="Please install via `conda install openmm -c omnia`.",
        )

        return rdkit_found and openff_found and openmm_found

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs OpenMM on given structure, inputs, in vacuum.
        """
        self.found(raise_error=True)

        from simtk import openmm
        from simtk import unit

        with capture_stdout():
            import openforcefield.topology as offtop

        # Failure flag
        ret_data = {"success": False}

        # generate basis, not given
        if not input_data.model.basis:
            raise InputError("Method must contain a basis set.")

        # Build a system based off input
        if input_data.model.basis.lower() == "smirnoff":

            ff_file = input_data.model.method + ".offxml"
            off_forcefield = self._get_off_forcefield(ff_file, ff_file)

            # Process molecule with RDKit
            rdkit_mol = RDKitHarness._process_molecule_rdkit(input_data.molecule)

            with capture_stdout():
                # Create an Open Force Field `Molecule` from the RDKit Molecule
                off_mol = offtop.Molecule(rdkit_mol)

                # Create OpenMM system in vacuum from forcefield, molecule
                off_top = off_mol.to_topology()
                openmm_system = self._get_openmm_system(off_forcefield, off_top)
        else:
            raise InputError("Accepted bases are: {'smirnoff', }")

        # Need an integrator for simulation even if we don't end up using it really
        integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)

        # Set platform to CPU explicitly
        platform = openmm.Platform.getPlatformByName("CPU")

        # Set number of threads to use
        # if `nthreads` is `None`, OpenMM default of all logical cores on
        # processor will be used
        nthreads = config.ncores
        if nthreads is None:
            nthreads = os.environ.get("OPENMM_CPU_THREADS")

        if nthreads:
            properties = {"Threads": str(nthreads)}
        else:
            properties = {}

        # Initialize context
        context = openmm.Context(openmm_system, integrator, platform, properties)

        # Set positions from our Open Force Field `Molecule`
        context.setPositions(off_mol.conformers[0])

        # Compute the energy of the configuration
        state = context.getState(getEnergy=True)

        # Get the potential as a simtk.unit.Quantity, put into units of hartree
        q = state.getPotentialEnergy() / unit.hartree / unit.AVOGADRO_CONSTANT_NA

        ret_data["properties"] = {"return_energy": q}

        # Execute driver
        if input_data.driver == "energy":
            ret_data["return_result"] = ret_data["properties"]["return_energy"]

        elif input_data.driver == "gradient":
            # Compute the forces
            state = context.getState(getForces=True)

            # Get the gradient as a simtk.unit.Quantity with shape (n_atoms, 3)
            gradient = state.getForces(asNumpy=True)

            # Convert to hartree/bohr and reformat as 1D array
            q = (gradient / (unit.hartree / unit.bohr)).reshape(-1) / unit.AVOGADRO_CONSTANT_NA

            # Force to gradient
            ret_data["return_result"] = -1 * q
        else:
            raise InputError(f"OpenMM can only compute energy and gradient driver methods. Found {input_data.driver}.")

        ret_data["success"] = True
        ret_data["extras"] = input_data.extras

        # Move several pieces up a level
        ret_data["provenance"] = Provenance(creator="openmm", version=openmm.__version__, nthreads=nthreads)

        return AtomicResult(**{**input_data.dict(), **ret_data})
