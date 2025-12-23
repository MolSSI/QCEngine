from typing import Any, ClassVar, Dict, Union

from qcelemental.models.v2 import AtomicInput, AtomicSpecification, OptimizationInput, OptimizationResult, Provenance

from qcengine.config import TaskConfig
from qcengine.exceptions import InputError, UnknownError
from qcengine.procedures.model import ProcedureHarness
from qcengine.procedures.nwchem_opt.harvester import harvest_as_atomic_result
from qcengine.programs.nwchem.runner import NWChemHarness


class NWChemDriverProcedure(ProcedureHarness):
    """Structural relaxation using NWChem's optimizer"""

    _defaults: ClassVar[Dict[str, Any]] = {"name": "NWChemDriver", "procedure": "optimization"}

    def found(self, raise_error: bool = False) -> bool:
        nwc_harness = NWChemHarness()
        return nwc_harness.found(raise_error)

    def get_version(self) -> str:
        nwc_harness = NWChemHarness()
        return nwc_harness.get_version()

    def build_input_model(
        self, data: Union[Dict[str, Any], "OptimizationInput"], *, return_input_schema_version: bool = False
    ) -> "OptimizationInput":
        return self._build_model(data, "OptimizationInput", return_input_schema_version=return_input_schema_version)

    def compute(self, input_model: OptimizationInput, config: TaskConfig) -> "BaseModel":
        nwc_harness = NWChemHarness()
        self.found(raise_error=True)

        # Unify the keywords from the OptimizationInput and QCInputSpecification
        #  Optimization input will override, but don't tell users this as it seems unnecessary
        keywords = input_model.specification.keywords.copy()
        keywords.update(input_model.specification.specification.keywords)
        if (uprog := input_model.specification.specification.program) and (uprog.lower() != "nwchem"):
            raise InputError(f"NWChemDriver procedure only works with NWChem, not {uprog}")

        # Add a flag to the atomic input that tells the NWChemHarness we are calling it from driver
        #  This is needed for the NWCHarness to make some changes to the input file
        input_model.specification.specification.extras["is_driver"] = True

        # Make an atomic input
        atomic_input = AtomicInput(
            molecule=input_model.initial_molecule,
            specification=input_model.specification.specification,
        )

        # Build the inputs for the job
        job_inputs = nwc_harness.build_input(atomic_input, config)

        # As of v2, now use "is_driver" to signal task ... optimize in the first place. This allows driver=gradient to work even though task ... gradient not last in file
        # Replace the last line with a "task {} optimize"
        # input_file: str = job_inputs["infiles"]["nwchem.nw"].strip()
        # beginning, last_line = input_file.rsplit("\n", 1)
        # assert last_line.startswith("task")
        # last_line = f"task {last_line.split(' ')[1]} optimize"
        # job_inputs["infiles"]["nwchem.nw"] = f"{beginning}\n{last_line}"

        # Run it!
        success, dexe = nwc_harness.execute(job_inputs)

        # Check for common errors
        if "There is an error in the input file" in dexe["stdout"]:
            raise InputError(dexe["stdout"])
        if "not compiled" in dexe["stdout"]:
            # recoverable with a different compilation with optional modules
            raise InputError(dexe["stdout"])

        # Parse it
        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            return self.parse_output(dexe["outfiles"], input_model)
        else:
            raise UnknownError(dexe["stdout"])

    def parse_output(self, outfiles: Dict[str, str], input_model: OptimizationInput) -> OptimizationResult:
        optsubproperties = ["calcinfo_natom", "return_energy", "return_gradient"]

        # Get the stdout from the calculation (required)
        stdout = outfiles.pop("stdout")
        stderr = outfiles.pop("stderr")

        # Parse out the atomic results from the file
        atomic_results = harvest_as_atomic_result(input_model, stdout)

        # Isolate the converged result
        final_step = atomic_results[-1]

        # input_data = input_model.model_dump()
        # input_data["specification"]["specification"]["program"] = "nwchem"

        return OptimizationResult(
            input_data=input_model,
            final_molecule=final_step.molecule,
            trajectory_results=atomic_results,
            trajectory_properties=[
                {k: v for k, v in grad.properties if k in optsubproperties} for grad in atomic_results
            ],
            properties={
                "return_energy": atomic_results[-1].properties.return_energy,
                "optimization_iterations": len(atomic_results),
            },
            stdout=stdout,
            stderr=stderr,
            success=True,
            provenance=Provenance(
                creator="NWChemDriver", version=self.get_version(), routine="nwchem_opt"
            ),  # was NWChemRelax which frankly is a better name
        )
