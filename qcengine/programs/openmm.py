"""
Calls the OpenMM executable.
"""
import json
import os
from typing import Dict

from qcelemental.models import Result

from .model import ProgramHarness
from ..exceptions import InputError, RandomError, ResourceError, UnknownError
from ..util import execute, popen, temporary_directory


class OpenMMHarness(ProgramHarness):

    # TODO: verify correctness of these default params
    # where are these used? By QCFractal?
    _defaults = {
        "name": "OpenMM",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False, # don't think openmm can run across multiple nodes
        "managed_memory": True,
    }

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which_import('simtk.openmm',
                     return_bool=True,
                     raise_error=raise_error,
                     raise_msg='Please install via `conda install openmm -c omnia`.')

    def compute(self, input_data: 'ResultInput', config: 'JobConfig') -> 'Result':
        """
        Runs OpenMM on given structure, inputs
        """
        self.found(raise_error=True)

        from simtk import openmm
        from simtk.openmm import app
        from simtk import unit
        

        ret_data["success"] = False

        # TODO: SET WORKDIR TO SCRATCH
        # Location resolution order config.scratch_dir, /tmp
        parent = config.scratch_directory
        with temporary_directory(parent=parent, suffix="_openmm_scratch") as tmpdir:

            # TODO: check what kind of processing we need to do of this input
            jmol = input_data.molecule

            ff_files = input_data.force_field + input_data.water_model
            forcefield = app.ForceField(*ff_files)

            # TODO: add nonbondedMethod, nonbondedCutoff as input params
            system = forcefield.createSystem(jmol.topology,
                                             nonbondedMethod=app.PME,
                                             nonbondedCutoff=1*unit.nanometer,
                                             constraints=input_data.constraints)

            # TODO: continue working here
            integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

            # initialize simulation
            simulation = Simulation(pdb.topology, system, integrator, platform)
            simulation.context.setPositions(pdb.positions)
            simulation.minimizeEnergy()

            # set up reporters
            simulation.reporters.append(PDBReporter('output.pdb', 1000))
            simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
                    potentialEnergy=True, temperature=True))

            # execute
            try:
                simulation.step(input_data.steps)
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
