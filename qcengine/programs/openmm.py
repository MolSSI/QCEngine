"""
Calls the OpenMM executable.

Requires RDKit
"""
import datetime
import hashlib
import os
from typing import TYPE_CHECKING, Dict

import numpy as np
from qcelemental.models import AtomicResult, Provenance
from qcelemental.util import which_import

from ..exceptions import InputError
from ..util import capture_stdout
from .model import ProgramHarness
from .rdkit import RDKitHarness

if TYPE_CHECKING:
    from qcelemental.models import AtomicInput

    from ..config import TaskConfig


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

    # def _get_off_forcefield(self, hashstring, offxml):
    #
    #     from openforcefield.typing.engines import smirnoff
    #
    #     key = hashlib.sha256(hashstring.encode()).hexdigest()
    #
    #     # get forcefield from cache, build new one if not present
    #     off_forcefield = self._get_cache(key) if key in self._CACHE else smirnoff.ForceField(offxml)
    #
    #     # cache forcefield, no matter what
    #     # handles updating time touched, dropping items if cache too large
    #     self._cache_it(key, off_forcefield)
    #
    #     return off_forcefield

    # def _get_openmm_system(self, off_forcefield, off_top):
    #
    #     # key = f"{mol_hash}-{input_model.method}-{hash(forcefield_schema)}"
    #
    #     openmm_system = off_forcefield.create_openmm_system(off_top)
    #
    #     return openmm_system

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

        openmmff_found = which_import(
            ".generators",
            return_bool=True,
            raise_error=raise_error,
            package="openmmforcefields",
            raise_msg="Please install via `conda install openmmforcefields -c conda-forge`",
        )

        return rdkit_found and openff_found and openmm_found and openmmff_found

    def _generate_openmm_system(
        self, molecule: "offtop.Molecule", method: str, keywords: Dict = None
    ) -> "openmm.System":
        """
        Generate an OpenMM System object from the input molecule method and basis.
        """
        from openmmforcefields.generators import SystemGenerator
        from simtk import unit
        from simtk.openmm import app

        # create a hash based on the input options
        hashstring = molecule.to_smiles(isomeric=True, explicit_hydrogens=True, mapped=True) + method
        for value in keywords.values():
            hashstring += str(value)
        key = hashlib.sha256(hashstring.encode()).hexdigest()

        # now look for the system?
        if key in self._CACHE:
            system = self._get_cache(key)
        else:
            # make the system from the inputs
            # set up available options for openmm
            _constraint_types = {"hbonds": app.HBonds, "allbonds": app.AllBonds, "hangles": app.HAngles}
            _periodic_nonbond_types = {"ljpme": app.LJPME, "pme": app.PME, "ewald": app.Ewald}
            _non_periodic_nonbond_types = {"nocutoff": app.NoCutoff, "cutoffnonperiodic": app.CutoffNonPeriodic}

            if "constraints" in keywords:
                constraints = keywords["constraints"]
                try:
                    forcefield_kwargs = {"constraints": _constraint_types[constraints.lower()]}
                except (KeyError, AttributeError):
                    raise ValueError(
                        f"constraint '{constraints}' not supported, valid constraints are {_constraint_types.keys()}"
                    )
            else:
                forcefield_kwargs = None

            nonbondedmethod = keywords.get("nonbondedMethod", None)
            if nonbondedmethod is not None:
                if nonbondedmethod.lower() in _periodic_nonbond_types:
                    periodic_forcefield_kwargs = {"nonbondedMethod": _periodic_nonbond_types[nonbondedmethod.lower()]}
                    nonperiodic_forcefield_kwargs = None
                elif nonbondedmethod.lower() in _non_periodic_nonbond_types:
                    periodic_forcefield_kwargs = None
                    nonperiodic_forcefield_kwargs = {
                        "nonbondedMethod": _non_periodic_nonbond_types[nonbondedmethod.lower()]
                    }
                else:
                    raise ValueError(
                        f"nonbondedmethod '{nonbondedmethod}' not supported, valid nonbonded methods are periodic: {_periodic_nonbond_types.keys()}"
                        f" or non_periodic: {_non_periodic_nonbond_types.keys()}."
                    )
            else:
                periodic_forcefield_kwargs = None
                nonperiodic_forcefield_kwargs = None

            # now start the system generator
            system_generator = SystemGenerator(
                small_molecule_forcefield=method,
                forcefield_kwargs=forcefield_kwargs,
                nonperiodic_forcefield_kwargs=nonperiodic_forcefield_kwargs,
                periodic_forcefield_kwargs=periodic_forcefield_kwargs,
            )
            topology = molecule.to_topology()

            system = system_generator.create_system(topology=topology.to_openmm(), molecules=[molecule])
            self._cache_it(key, system)
        return system

    def compute(self, input_data: "AtomicInput", config: "TaskConfig") -> "AtomicResult":
        """
        Runs OpenMM on given structure, inputs, in vacuum.
        """
        self.found(raise_error=True)

        from simtk import openmm, unit

        with capture_stdout():
            from openforcefield import topology as offtop

        # Failure flag
        ret_data = {"success": False}

        # generate basis, not given
        if not input_data.model.basis:
            raise InputError("Method must contain a basis set.")

        # Make sure we are using smirnoff or antechamber
        basis = input_data.model.basis.lower()
        if basis in ["smirnoff", "antechamber"]:

            with capture_stdout():
                # try and make the molecule from the cmiles
                cmiles = None
                if input_data.molecule.extras:
                    cmiles = input_data.molecule.extras.get("canonical_isomeric_explicit_hydrogen_mapped_smiles", None)
                    if cmiles is None:
                        cmiles = input_data.molecule.extras.get("cmiles", {}).get(
                            "canonical_isomeric_explicit_hydrogen_mapped_smiles", None
                        )

                if cmiles is not None:
                    off_mol = offtop.Molecule.from_mapped_smiles(mapped_smiles=cmiles)
                    # add the conformer
                    conformer = unit.Quantity(value=np.array(input_data.molecule.geometry), unit=unit.bohr)
                    off_mol.add_conformer(conformer)
                else:
                    # Process molecule with RDKit
                    rdkit_mol = RDKitHarness._process_molecule_rdkit(input_data.molecule)

                    # Create an Open Force Field `Molecule` from the RDKit Molecule
                    off_mol = offtop.Molecule(rdkit_mol)

            # now we need to create the system
            openmm_system = self._generate_openmm_system(
                molecule=off_mol, method=input_data.model.method, keywords=input_data.keywords
            )
        else:
            raise InputError("Accepted bases are: {'smirnoff', 'antechamber', }")

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
