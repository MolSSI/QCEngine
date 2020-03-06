from typing import Any, Dict, Union

from qcelemental.models import OptimizationInput, OptimizationResult, AtomicInput
from qcelemental.util import parse_version
from ..exceptions import InputError, UnknownError

from .model import ProcedureHarness
from .. import programs


class QChemOptimizationProcedure(ProcedureHarness):
    # Is there a pythonic way to get these from programs.qchem.QChemHarness?
    _defaults: Dict[str, Any] = {
        "name": "QChem",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
        "procedure": "optimization",
    }

    class Config(ProcedureHarness.Config):
        pass

    def __init__(self, **kwargs):
        self.qchem_harness = programs.get_program("qchem", False)
        super().__init__(**{**self._defaults, **kwargs})

    def found(self, raise_error: bool = False) -> bool:
        return self.qchem_harness.found(raise_error)

    def build_input_model(self, data: Union[Dict[str, Any], "OptimizationInput"]) -> "OptimizationInput":
        return self._build_model(data, OptimizationInput)

    @staticmethod
    def build_atomic_input(input_model: "OptimizationInput") -> AtomicInput:
        keywords = input_model.keywords.copy()
        keywords.update(input_model.input_specification.keywords)
        extras = input_model.keywords.copy()
        extras.update(input_model.input_specification.extras)
        extras["_procedure"] = "optimization"
        return AtomicInput(
                molecule=input_model.initial_molecule,
                driver=input_model.input_specification.driver,
                model=input_model.input_specification.model,
                keywords=keywords,
                extras=extras,
            )

    def compute(self, input_data: "OptimizationInput", config: "TaskConfig") -> "OptimizationResult":
        self.found(raise_error=True)

        # build input file
        # run using qchem.execute
        # parse logfile repeatedly
        atomic_input = self.build_atomic_input(input_data)
        job_inputs = self.qchem_harness.build_input(atomic_input, config)

        # Check qchem version
        qceng_ver = "5.2"
        if parse_version(self.qchem_harness.get_version()) < parse_version(qceng_ver):
            raise TypeError(f"Q-Chem version <{qceng_ver} not supported (found version {self.get_version()})")
        # Run qchem
        exe_success, proc = self.qchem_harness.execute(job_inputs)

        # Determine whether the calculation succeeded
        if exe_success:
            # If execution succeeded, collect results
            result = self.parse_output(proc["outfiles"], input_data)
            return result
        else:
            outfile = proc["outfiles"]["dispatch.out"]
            if "fatal error occurred in module qparser" in outfile:
                raise InputError(proc["outfiles"]["dispatch.out"])
            else:
                # Return UnknownError for error propagation
                raise UnknownError(proc["outfiles"]["dispatch.out"])

    def parse_output(self, outfiles: Dict[str, str], input_model: "OptimizationInput"):
        atomic_input = self.build_atomic_input(input_model)
        atomic_input.extras.pop("_procedure")

