"""Compute quantum chemistry using ChronusQ executable."""

import io
from typing import Any, ClassVar, Dict, Optional, Tuple

import numpy as np
from qcelemental.models.v2 import AtomicInput, AtomicResult, BasisSet, Provenance
from qcelemental.util import safe_version, which

from ..exceptions import InputError
from ..util import execute
from .model import ProgramHarness


class ChronusQHarness(ProgramHarness):
    """Interface for ChronusQ project."""

    _defaults: ClassVar[Dict[str, Any]] = {
        "name": "ChronusQ",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "chronusq", return_bool=True, raise_error=raise_error, raise_msg="Please install via https://github.com/xsligroup/chronusq_public"
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("chronusq")
        if which_prog not in self.version_cache:
            success, output = execute([which_prog, "test.inp"], {"test.inp": """[Molecule]
charge = 0
mult = 1
geom:
 He               0  -0.07579184359               0
[QM]
reference = Real RHF
job = SCF
[BASIS]
basis = sto-3g  
[MISC]
memtype = os"""}, ["test.out"])

            if success:
                for line in output["outfiles"]["test.out"].splitlines():
                    if "Release Version:" in line:
                        branch = " ".join(line.strip().split()[2:])
                        break
                self.version_cache[which_prog] = safe_version(branch)

        return self.version_cache[which_prog]

    def compute(self, input_model: AtomicInput, config: "TaskConfig") -> AtomicResult:
        self.found(raise_error=True)

        job_inputs = self.build_input(input_model, config)
        success, dexe = self.execute(job_inputs)

        if success:
            dexe["outfiles"]["stdout"] = dexe["stdout"]
            dexe["outfiles"]["stderr"] = dexe["stderr"]
            dexe["outfiles"]["input"] = job_inputs["infiles"]["test.inp"]
            return self.parse_output(dexe["outfiles"], input_model)

    def build_input(
        self, input_model: AtomicInput, config: "TaskConfig", template: Optional[str] = None
    ) -> Dict[str, Any]:
        ChronusQrec = {
            "infiles": {},
            "scratch_directory": config.scratch_directory,
            "scratch_messy": config.scratch_messy,
        }

        # Handle molecule
        molcmd, moldata = input_model.molecule.to_string(dtype="ChronusQ", units="Angstrom", return_data=True)

        # Handle basis set
        if isinstance(input_model.specification.model.basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use string basis name.")
        if input_model.specification.model.basis is None:
            raise InputError("None for model.basis is not useable.")

        bascmd = ["[BASIS]", f"basis = {input_model.specification.model.basis}"]
        bascmd = "\n".join(bascmd)

        tempcmd = """
[QM]
reference = Real RHF
job = SCF
[MISC]
memtype = os"""

        # Handle conversion from schema (flat key/value) keywords into local format
        #optcmd = format_keywords(opts)

        ChronusQ = which("ChronusQ")
        ChronusQrec["infiles"]["test.inp"] = molcmd + bascmd + tempcmd
        ChronusQrec["command"] = [ChronusQ, "test.inp"]
        print(ChronusQrec["infiles"]["test.inp"])

        return ChronusQrec

    def execute(
        self, inputs: Dict[str, Any], *, extra_outfiles=None, extra_commands=None, scratch_name=None, timeout=None
    ) -> Tuple[bool, Dict]:

        # llel works b/c util.environ_context sets OMP_NUM_THREADS = config.ncores

        success, dexe = execute(
            inputs["command"],
            inputs["infiles"],
            ["test.out", "test.bin"],
            as_binary=["test.bin"],
            scratch_messy=True,
            scratch_directory=inputs["scratch_directory"],
        )
        return success, dexe

    def parse_output(
        self, outfiles: Dict[str, str], input_model: AtomicInput
    ) -> AtomicResult:  # lgtm: [py/similar-function]

        provenance = Provenance(creator="ChronusQ", version=self.get_version(), routine="ChronusQ").model_dump()

        output_data = {
            "schema_version": 2,
            "input_data": input_model,
            "molecule": input_model.molecule,  # overwrites with outfile Cartesians in case fix_*=F
            "extras": {},
            "native_files": {k: v for k, v in outfiles.items() if v is not None},
            "properties": {},
            "provenance": provenance,
            "return_result": 0,
            "success": True,
        }

        outbin = outfiles["test.bin"]

        # Parse HDF5 binary output file
        if outbin is not None:
            try:
                import h5py
            except ImportError:
                raise ImportError("h5py is required to parse ChronusQ binary output. Install with: pip install h5py")

            with h5py.File(io.BytesIO(outbin), "r") as hf:
                # -- Total energy --
                if "/SCF/TOTAL_ENERGY" in hf:
                    output_data["return_result"] = np.asarray(hf["/SCF/TOTAL_ENERGY"][()]).item()

        return AtomicResult(**output_data)
