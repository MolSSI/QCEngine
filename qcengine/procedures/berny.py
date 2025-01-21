import logging
import sys
import traceback
from io import StringIO
from typing import Any, ClassVar, Dict, Union

import numpy as np
from qcelemental.models.v2 import FailedOperation, Molecule, OptimizationInput, OptimizationResult
from qcelemental.util import which_import

import qcengine
from qcengine.exceptions import UnknownError

from ..config import TaskConfig
from .model import ProcedureHarness


class BernyProcedure(ProcedureHarness):
    _defaults: ClassVar[Dict[str, Any]] = {"name": "Berny", "procedure": "optimization"}

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "berny",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `pip install pyberny`.",
        )

    def build_input_model(
        self, data: Union[Dict[str, Any], "OptimizationInput"], *, return_input_schema_version: bool = False
    ) -> "OptimizationInput":
        return self._build_model(data, "OptimizationInput", return_input_schema_version=return_input_schema_version)

    def compute(
        self, input_model: "OptimizationInput", config: "TaskConfig"
    ) -> Union["OptimizationResult", "FailedOperation"]:
        try:
            import berny
        except ModuleNotFoundError:
            raise ModuleNotFoundError("Could not find Berny in the Python path.")

        # Get berny version from the installed package, use setuptools'
        # pkg_resources for python < 3.8
        if sys.version_info >= (3, 8):
            from importlib.metadata import distribution
        else:
            from pkg_resources import get_distribution as distribution
        berny_version = distribution("pyberny").version

        # Berny uses the stdlib logging module and by default uses per-module
        # loggers. For QCEngine, we create one logger per BernyProcedure
        # instance, by using the instance's id(), and send all logging messages
        # to a string stream
        log_stream = StringIO()
        log = logging.getLogger(f"{__name__}.{id(self)}")
        log.addHandler(logging.StreamHandler(log_stream))
        log.setLevel("INFO")
        input_data = input_model.model_dump()

        geom_qcng = input_data["initial_molecule"]
        comput = {"specification": input_data["specification"]["specification"], "molecule": geom_qcng}
        program = input_data["specification"]["specification"][
            "program"
        ]  # TODO don't need to collect when compute can work w/o 2nd arg
        task_config = config.dict()
        trajectory = []
        output_data = input_data.copy()
        try:
            # Pyberny uses angstroms for the Cartesian geometry, but atomic
            # units for everything else, including the gradients (hartree/bohr).
            geom_berny = berny.Geometry(geom_qcng["symbols"], geom_qcng["geometry"] / berny.angstrom)
            opt = berny.Berny(geom_berny, logger=log, **input_data["specification"]["keywords"])
            for geom_berny in opt:
                geom_qcng["geometry"] = np.stack(geom_berny.coords * berny.angstrom)
                ret = qcengine.compute(comput, program, task_config=task_config)
                if ret.success:
                    trajectory.append(ret.dict())
                    opt.send((ret.properties.return_energy, ret.return_result))
                else:
                    # qcengine.compute returned FailedOperation
                    raise UnknownError("Gradient computation failed")

        except UnknownError:
            error = ret.error.dict()  # ComputeError
        except Exception:
            error = {"error_type": "unknown", "error_message": f"Berny error:\n{traceback.format_exc()}"}
        else:
            final_molecule = trajectory[-1]["molecule"]
            output = {
                "input_data": input_model,
                "final_molecule": final_molecule,
                "properties": {
                    "nuclear_repulsion_energy": Molecule(**final_molecule).nuclear_repulsion_energy(),
                    "return_energy": trajectory[-1]["properties"]["return_energy"],
                    "return_gradient": trajectory[-1]["properties"]["return_gradient"],
                    "optimization_iterations": len(trajectory),
                },
                "trajectory_results": trajectory,
                "trajectory_properties": [r["properties"] for r in trajectory],
                "provenance": {"creator": "Berny", "routine": "berny.Berny", "version": berny_version},
                "stdout": log_stream.getvalue(),  # collect logged messages
                "success": True,
            }

            return OptimizationResult(**output)

        return FailedOperation(input_data=input_model, error=error)
