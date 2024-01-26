import abc
import copy
import pprint
from ast import literal_eval
from typing import TYPE_CHECKING, Any, ClassVar, Dict, List, Tuple, Union, Literal, Optional

# printing and logging formatting niceties
from functools import partial
import numpy as np
pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
nppp = partial(np.array_str, max_line_width=120, precision=8, suppress_small=True)
nppp10 = partial(np.array_str, max_line_width=120, precision=10, suppress_small=True)

#import numpy as np
from pydantic import ConfigDict, field_validator, FieldValidationInfo, computed_field, BaseModel, Field
from qcelemental.models import FailedOperation, Molecule, DriverEnum, ProtoModel, AtomicResult, AtomicInput
from qcelemental.models.procedures_layered import BsseEnum, ManyBodyKeywords, ManyBodyInput, ManyBodyResult, ManyBodyResultProperties
from qcelemental.util import safe_version, which_import

from .model import ProcedureHarness
from ..extras import provenance_stamp
from ..exceptions import UnknownError

#if TYPE_CHECKING:
#    from ..config import TaskConfig


class ManyBodyProcedure(ProcedureHarness):

    _defaults: ClassVar[Dict[str, Any]] = {"name": "ManyBody", "procedure": "manybody"}

    version_cache: Dict[str, str] = {}

    def found(self, raise_error: bool = False) -> bool:
        return which_import(
            "psi4",  # since we're borrowing fns for now
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via `conda install psi4 -c conda-forge/label/libint_dev -c conda-forge`.",
        )

    def build_input_model(self, data: Union[Dict[str, Any], "ManyBodyInput"]) -> "ManyBodyInput":
        return self._build_model(data, ManyBodyInput)

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which_import("psi4")
        if which_prog not in self.version_cache:
            import psi4

            self.version_cache[which_prog] = safe_version(psi4.__version__)

        return self.version_cache[which_prog]

########
    def compute(self, input_model: "ManyBodyInput", config: "TaskConfig") -> "ManyBodyResult":

        self.found(raise_error=True)

        import psi4
        print("\n User Keywords (all)", input_model.specification.keywords.model_dump())
        print("\n User Keywords (exclude_unset)", input_model.specification.keywords.model_dump(exclude_unset=True))

#    # Build a packet
#    packet = {"molecule": molecule, "driver": driver, "method": method, "basis": basis, "keywords": keywords}
#
#    # First check for BSSE type
#    if kwargs.get("bsse_type", None) is not None:
#        levels = kwargs.pop('levels', None)
#
#        plan = ManyBodyComputer(**packet, **kwargs)
#        original_molecule = packet.pop("molecule")

        plan = ManyBodyComputerQCNG(
            molecule=input_model.molecule,
            driver=input_model.specification.driver,
            **input_model.specification.keywords.model_dump(exclude_unset=True),
            input_data=input_model,  # storage, to reconstitute ManyBodyResult
        )

        print("\n<<<  (2) QCEngine harness ManyBodyComputerQCNG init  >>>")
        print(type(plan), isinstance(plan, BaseModel))
        pprint.pprint(plan.model_dump(), width=200)

        for mc_level_idx, mtd in enumerate(plan.levels.values()):
            # TODO method, basis, cbsmeta = expand_cbs_methods(mtd, basis, driver, cbsmeta=cbsmeta, **kwargs)

            # analytic
            #logger.info(f"PLANNING MB:  {mc_level_idx=} {packet=}")
            plan.build_tasks(
                AtomicComputerQCNG,
                mc_level_idx=mc_level_idx,
                program=input_model.specification.specification.program,
                driver=input_model.specification.specification.driver,
                method=input_model.specification.specification.model.method,
                basis=input_model.specification.specification.model.basis,
                keywords=input_model.specification.specification.keywords,
                #**kwargs,
            )

        print("\n<<<  (3) QCEngine harness ManyBodyComputerQCNG build_levels  >>>")
        pp.pprint(plan.model_dump())

        plan.compute()

        try:
            results_list = {k: v.get_results() for k, v in plan.task_list.items()}
            all_good = all(v.success for v in results_list.values())

            if not all_good:
                for k, v in results_list.items():
                    if not v.success:
                        ret = v
                        raise UnknownError("ManyBody component computation failed")
        except UnknownError:
            error = ret.error.model_dump()  # ComputeError
        except Exception:
            error = {"error_type": "unknown", "error_message": f"Berny error:\n{traceback.format_exc()}"}
        else:
            output_model = plan.get_results()

            print("\n<<<  (4) QCEngine harness ManyBodyComputerQCNG get_results >>>")
            pp.pprint(output_model.model_dump())

            return output_model
        return FailedOperation(input_data=input_model, error=error)



#        output_data = input_model
#        input_data = input_model.model_dump()
#
#        # Set retries to two if zero while respecting local_config
#        local_config = config.model_dump()
#        local_config["retries"] = local_config.get("retries", 2) or 2
#        input_data["input_specification"]["extras"]["_qcengine_local_config"] = local_config
#
#        # Run the program
#        output_data = optking.optwrapper.optimize_qcengine(input_data)
#
#        output_data["schema_name"] = "qcschema_optimization_output"
#        output_data["input_specification"]["extras"].pop("_qcengine_local_config", None)
#        if output_data["success"]:
#            output_data = ManyBodyResult(**output_data)
#
#        return output_data

########

#    def _compute(self, input_model: "TorsionDriveInput", config: "TaskConfig"):
#
#        self.found(raise_error=True)
#
#        import torsiondrive.td_api
#
#        dihedrals = input_model.keywords.dihedrals
#        grid_spacing = input_model.keywords.grid_spacing
#
#        dihedral_ranges = input_model.keywords.dihedral_ranges
#
#        energy_decrease_thresh = input_model.keywords.energy_decrease_thresh
#        energy_upper_limit = input_model.keywords.energy_upper_limit
#
#        state = torsiondrive.td_api.create_initial_state(
#            dihedrals=dihedrals,
#            grid_spacing=grid_spacing,
#            elements=input_model.initial_molecule[0].symbols,
#            init_coords=[molecule.geometry.flatten().tolist() for molecule in input_model.initial_molecule],
#            dihedral_ranges=dihedral_ranges,
#            energy_upper_limit=energy_upper_limit,
#            energy_decrease_thresh=energy_decrease_thresh,
#        )
#
#        optimization_results = defaultdict(list)
#        error = None
#
#        # Spawn new optimizations at each grid points until convergence / an error.
#        while True:
#
#            next_jobs = torsiondrive.td_api.next_jobs_from_state(state, verbose=False)
#
#            if len(next_jobs) == 0:
#                break
#
#            grid_point_results = self._spawn_optimizations(next_jobs=next_jobs, input_model=input_model, config=config)
#
#            for grid_point, results in grid_point_results.items():
#
#                failed_results = [result for result in results if not result.success]
#
#                if len(failed_results) > 0:
#
#                    error_message = failed_results[0].error.error_message
#                    error = {
#                        "error_type": "unknown",
#                        "error_message": f"TorsionDrive error at {grid_point}:\n{error_message}",
#                    }
#                    break
#
#                optimization_results[grid_point].extend(results)
#
#            if error is not None:
#                break
#
#            task_results = {
#                grid_point: [
#                    (
#                        result.initial_molecule.geometry.flatten().tolist(),
#                        result.final_molecule.geometry.flatten().tolist(),
#                        result.energies[-1],
#                    )
#                    for result in results
#                ]
#                for grid_point, results in grid_point_results.items()
#            }
#
#            torsiondrive.td_api.update_state(state, {**task_results})
#
#        output_data = input_model.model_dump()
#        output_data["provenance"] = {
#            "creator": "TorsionDrive",
#            "routine": "torsiondrive.td_api.next_jobs_from_state",
#            "version": torsiondrive.__version__,
#        }
#        output_data["success"] = error is None
#
#        # even if we hit an error during the torsiondrive, we output what we can
#        output_data["final_energies"], output_data["final_molecules"] = {}, {}
#
#        for grid_point, results in optimization_results.items():
#
#            final_energy, final_molecule = self._find_final_results(results)
#
#            output_data["final_energies"][grid_point] = final_energy
#            output_data["final_molecules"][grid_point] = final_molecule
#
#        output_data["optimization_history"] = optimization_results
#
#        if error is not None:
#            output_data["error"] = error
#
#        return output_data
#
#    def compute(self, input_model: "TorsionDriveInput", config: "TaskConfig") -> "TorsionDriveResult":
#
#        # Capture the stdout and err here to avoid too much nesting in the _compute function
#        with io.StringIO() as stdout:
#            with io.StringIO() as stderr:
#
#                with redirect_stdout(stdout):
#                    with redirect_stderr(stderr):
#
#                        output_data = self._compute(input_model, config)
#
#                output_data["stdout"] = str(stdout.getvalue())
#                output_data["stderr"] = str(stderr.getvalue())
#
#        # these will get populated by the model below
#        output_data.pop("schema_name", None)
#        output_data.pop("schema_version", None)
#
#        output_data = TorsionDriveResult(**output_data)
#
#        return output_data
#
#    @staticmethod
#    def _spawn_optimization(
#        grid_point: str, job: List[float], input_model: "TorsionDriveInput", config: "TaskConfig"
#    ) -> Union[FailedOperation, OptimizationResult]:
#        """Spawns an optimization at a particular grid point and returns the result.
#
#        Parameters
#        ----------
#        grid_point
#            A string of the form 'dihedral_1_angle ... dihedral_n_angle' that encodes
#            the current dihedrals angles to optimize at.
#        job
#            The flattened conformer of the molecule to start the optimization at with
#            length=(n_atoms * 3)
#        input_model
#            The input model containing the relevant settings for how to optimize the
#            structure.
#        config
#            The configuration to launch the task using.
#
#        Returns
#        -------
#            The result of the optimization if successful, otherwise an error containing
#            object.
#        """
#
#        from qcengine import compute_procedure
#
#        input_molecule = input_model.initial_molecule[0].model_copy(deep=True).model_dump()
#        input_molecule["geometry"] = np.array(job).reshape(len(input_molecule["symbols"]), 3)
#        input_molecule = Molecule.from_data(input_molecule)
#
#        dihedrals = input_model.keywords.dihedrals
#        angles = grid_point.split()
#
#        keywords = {
#            **input_model.optimization_spec.keywords,
#            "constraints": {
#                "set": [
#                    {
#                        "type": "dihedral",
#                        "indices": dihedral,
#                        "value": int(angle),
#                    }
#                    for dihedral, angle in zip(dihedrals, angles)
#                ]
#            },
#        }
#
#        input_data = OptimizationInput(
#            keywords=keywords,
#            extras={},
#            protocols=input_model.optimization_spec.protocols,
#            input_specification=input_model.input_specification,
#            initial_molecule=input_molecule,
#        )
#
#        return compute_procedure(
#            input_data, procedure=input_model.optimization_spec.procedure, task_config=config.model_dump()
#        )
#
#    @staticmethod
#    def _find_final_results(
#        optimization_results: List[OptimizationResult],
#    ) -> Tuple[float, Molecule]:
#        """Returns the energy and final molecule of the lowest energy optimization
#        in a set."""
#
#        final_energies = np.array([result.energies[-1] for result in optimization_results])
#        lowest_energy_idx = final_energies.argmin()
#
#        return float(final_energies[lowest_energy_idx]), optimization_results[lowest_energy_idx].final_molecule
#
#    def _spawn_optimizations(
#        self, next_jobs: Dict[str, List[float]], input_model: "TorsionDriveInput", config: "TaskConfig"
#    ) -> Dict[str, List[Union[FailedOperation, OptimizationResult]]]:
#
#        grid_point_results = {
#            grid_point: [self._spawn_optimization(grid_point, job, input_model, config) for job in jobs]
#            for grid_point, jobs in next_jobs.items()
#        }
#        return grid_point_results


class BaseComputerQCNG(ProtoModel):
    """Base class for "computers" that plan, run, and process QC tasks."""

    @abc.abstractmethod
    def compute(self):
        pass

    @abc.abstractmethod
    def plan(self):
        pass

    # TODO can remove?
    model_config = ConfigDict(
        extra="allow",
        frozen=False,
    )


class AtomicComputerQCNG(BaseComputerQCNG):
    """Computer for analytic single-geometry computations."""

    molecule: Molecule = Field(..., description="The molecule to use in the computation.")
    basis: str = Field(..., description="The quantum chemistry basis set to evaluate (e.g., 6-31g, cc-pVDZ, ...).")
    method: str = Field(..., description="The quantum chemistry method to evaluate (e.g., B3LYP, MP2, ...).")
    driver: DriverEnum = Field(..., description="The resulting type of computation: energy, gradient, hessian, properties."
        "Note for finite difference that this should be the target driver, not the means driver.")
    keywords: Dict[str, Any] = Field(default_factory=dict, description="The keywords to use in the computation.")
    program: str = Field(..., description="Which program harness to run single-point with.")
    computed: bool = Field(False, description="Whether quantum chemistry has been run on this task.")
    result: Any = Field(default_factory=dict, description=":py:class:`~qcelemental.models.AtomicResult` return.")
    result_id: Optional[str] = Field(None, description="The optional ID for the computation.")

    @field_validator("basis")
    @classmethod
    def set_basis(cls, basis):
        return basis.lower()

    @field_validator("method")
    @classmethod
    def set_method(cls, method):
        return method.lower()

    @field_validator("keywords")
    @classmethod
    def set_keywords(cls, keywords):
        return copy.deepcopy(keywords)

    def plan(self) -> AtomicInput:
        """Form QCSchema input from member data."""

        atomic_model = AtomicInput(**{
            "molecule": self.molecule,
            "driver": self.driver,
            "model": {
                "method": self.method,
                "basis": self.basis
            },
            "keywords": self.keywords,
#            "protocols": {
#                "stdout": True,
#            },
#            "extras": {
#                "psiapi": True,
#                "wfn_qcvars_only": True,
#            },
        })

        return atomic_model

    def compute(self) -> None:
        """Run quantum chemistry single-point.

        NOTE: client removed (compared to psi4.driver.AtomicComputer)
        """
        from ..compute import compute as qcng_compute

        if self.computed:
            return

        #logger.info(f'<<< JSON launch ... {self.molecule.name} {self.molecule.nuclear_repulsion_energy()}')

        self.result = qcng_compute(
            self.plan(),
            self.program,
            raise_error=False,  #True,
            #task_config=task_config,
        )

        #pp.pprint(self.result.model_dump())
        #logger.debug(pp.pformat(self.result.model_dump()))
        self.computed = True

    def get_results(self) -> AtomicResult:
        """Return results as Atomic-flavored QCSchema.

        NOTE: client removed (compared to psi4.driver.AtomicComputer)
        """
        if self.result:
            return self.result


class ManyBodyComputerQCNG(BaseComputerQCNG):

    input_data: ManyBodyInput = Field(
        ...,
        description="Input schema containing the relevant settings for performing the many body "
            "expansion. This is entirely redundant with the piecemeal assembly of this Computer class "
            "and is only stored to be available for error handling and exact reconstruction of ManyBodyResult.",
    )
    bsse_type: List[BsseEnum] = Field(
        [BsseEnum.cp],
        description=ManyBodyKeywords.model_fields["bsse_type"].description,
    )
    molecule: Molecule = Field(
        ...,
        description="Target molecule for many body expansion (MBE) or interaction energy (IE) analysis. "
            "Fragmentation should already be defined in `fragments` and related fields.",
    )
    driver: DriverEnum = Field(
        ...,
        description="The computation driver; i.e., energy, gradient, hessian. In case of ambiguity (e.g., MBE gradient "
            "through finite difference energies or MBE energy through composite method), this field refers to the "
            "*target* derivative, not any *means* specification.",
    )
    embedding_charges: Dict[int, List[float]] = Field(
        {},
        description="Atom-centered point charges to be used on molecule fragments whose basis sets are not included in "
            "the computation. Keys: 1-based index of fragment. Values: list of atom charges for that fragment.",
        json_schema_extra={
            "shape": ["nfr", "<varies: nat in ifr>"],
        },
    )
    return_total_data: Optional[bool] = Field(  # after driver, embedding_charges
        None,
        validate_default=True,
        description=ManyBodyKeywords.model_fields["return_total_data"].description,
    )
    levels: Optional[Dict[Union[int, Literal["supersystem"]], str]] = Field(
        None,
        validate_default=True,
        description=ManyBodyKeywords.model_fields["levels"].description + \
            "Examples above are processed in the ManyBodyComputer, and once processed, only the values should be used. "
            "The keys turn into nbodies_per_mc_level, as notated below. "
            "* {1: 'ccsd(t)', 2: 'mp2', 'supersystem': 'scf'} -> nbodies_per_mc_level=[[1], [2], ['supersystem']] "
            "* {2: 'ccsd(t)/cc-pvdz', 3: 'mp2'} -> nbodies_per_mc_level=[[1, 2], [3]] ",
    )
    max_nbody: Optional[int] = Field(
        None,
        validate_default=True,
        description=ManyBodyKeywords.model_fields["max_nbody"].description,
    )
    supersystem_ie_only: Optional[bool] = Field(  # after max_nbody
        False,
        validate_default=True,
        description=ManyBodyKeywords.model_fields["supersystem_ie_only"].description,
    )
    task_list: Dict[str, Any] = {}  #MBETaskComputers] = {}

    @computed_field(description="Number of distinct fragments comprising full molecular supersystem.")
    @property
    def nfragments(self) -> int:
        return len(self.molecule.fragments)

    @field_validator("bsse_type", mode="before")
    @classmethod
    def set_bsse_type(cls, v: Any) -> List[BsseEnum]:
        if not isinstance(v, list):
            v = [v]
        # emulate ordered set
        return list(dict.fromkeys([bt.lower() for bt in v]))

    @field_validator("embedding_charges")
    @classmethod
    def set_embedding_charges(cls, v: Any, info: FieldValidationInfo) -> Dict[int, List[float]]:
        if len(v) != info.data["nfragments"]:
            raise ValueError("embedding_charges dict should have entries for each 1-indexed fragment.")

        return v

    @field_validator("return_total_data")
    @classmethod
    def set_return_total_data(cls, v: Optional[bool], info: FieldValidationInfo) -> bool:
        print(f"hit return_total_data validator with {v}", end="")
        if v is not None:
            rtd = v
        elif info.data["driver"] in ["gradient", "hessian"]:
            rtd = True
        else:
            rtd = False

        if info.data.get("embedding_charges", False) and rtd is False:
            raise ValueError("Cannot return interaction data when using embedding scheme.")

        print(f" ... setting rtd={rtd}")
        return rtd

    @field_validator("levels")
    @classmethod
    def set_levels(cls, v: Any, info: FieldValidationInfo) -> Dict[Union[int, Literal["supersystem"]], str]:
        print(f"hit levels validator with {v}", end="")

        if v is None:
            pass
            # TODO levels = {plan.max_nbody: method}
            #v = {info.data["nfragments"]: "???method???"}
            v = {len(info.data["molecule"].fragments): "???method???"}
        else:
            # rearrange bodies in order with supersystem last lest body count fail in organization loop below
            v = dict(sorted(v.items(), key=lambda item: 1000 if item[0] == "supersystem" else item[0]))

            # rm 1 for cp-only
            # We define cp as being a correction to only interaction energies
            # If only doing cp, we need to ignore any user-specified 1st (monomer) level
            #if 'cp' in kwargs.get("bsse_type", None) and 'nocp' not in kwargs.get("bsse_type", None):
            #    if 1 in levels.keys():
            #        removed_level = levels.pop(1)
            #        logger.info("NOTE: User specified exclusively 'cp' correction, but provided level 1 details")
            #        logger.info(f"NOTE: Removing level {removed_level}")
            #        logger.info("NOTE: For total energies, add 'nocp' to bsse_list")

        print(f" ... setting levels={v}")
        return v

    @computed_field(
        description="Distribution of active n-body levels among model chemistry levels. All bodies in range "
            "[1, self.max_nbody] must be present exactly once. Number of items in outer list is how many different "
            "modelchems. Each inner list specifies what n-bodies to be run at the corresponding modelchem (e.g., "
            "`[[1, 2]]` has max_nbody=2 and 1-body and 2-body contributions computed at the same level of theory; "
            "`[[1], [2]]` has max_nbody=2 and 1-body and 2-body contributions computed at different levels of theory. "
            "An entry 'supersystem' means all higher order n-body effects up to the number of fragments. The n-body "
            "levels are effectively sorted in the outer list, and any 'supersystem' element is at the end.")
        #json_schema_extra={
        #    "shape": ["nmc", "<varies>"],
        #},
    @property
    def nbodies_per_mc_level(self) -> List[List[Union[int, Literal["supersystem"]]]]:
        print(f"hit nbodies_per_mc_level", end="")

        # Organize nbody calculations into modelchem levels
        # * expand keys of `levels` into full lists of nbodies covered. save to plan, resetting max_nbody accordingly
        # * below, process values of `levels`, which are modelchem strings, into kwargs specs
        nbodies_per_mc_level = []
        prev_body = 0
        print("\nAAA", self.levels)
        for nb in self.levels:
            nbodies = []
            print("BBB bfore", nb, nbodies, prev_body)
            if nb == "supersystem":
                nbodies.append(nb)
            elif nb != (prev_body + 1):
                for m in range(prev_body + 1, nb + 1):
                    nbodies.append(m)
            else:
                nbodies.append(nb)
            print("BBB after", nb, nbodies)
            nbodies_per_mc_level.append(nbodies)
            prev_body = nb  # formerly buggy `+= 1`

        print(f" ... setting nbodies_per_mc_level={nbodies_per_mc_level}")
        return nbodies_per_mc_level

    @field_validator("max_nbody")
    @classmethod
    def set_max_nbody(cls, v: Any, info: FieldValidationInfo) -> int:
        print(f"hit max_nbody validator with {v}", end="")
        levels_max_nbody = max(nb for nb in info.data["levels"] if nb != "supersystem")
        nfr = len(info.data["molecule"].fragments)

        #ALT if v == -1:
        if v is None:
            v = levels_max_nbody
        elif v < 0 or v > nfr:
            raise ValueError(f"max_nbody={v} should be between 1 and {nfr}.")
        elif v != levels_max_nbody:
            #raise ValueError(f"levels={levels_max_nbody} contradicts user max_nbody={v}.")
            # TODO reconsider logic. move this from levels to here?
            info.data["levels"] = {v: "???method???"}
        else:
            pass
            # TODO once was           return min(v, nfragments)

        print(f" ... setting max_nbody={v}")
        return v

#       levels          max_nbody           F-levels        F-max_nbody     result
#
#       {stuff}         None                {stuff}         set from stuff  all consistent; max_nbody from levels
#       None            int                 {int: mtd}      int             all consistent; levels from max_nbody
#       None            None                {nfr: mtd}      nfr             all consistent; any order
#       {stuff}         int                 {stuff}         int             need to check consistency

    # TODO also, perhaps change nbodies_per_mc_level into dict of lists so that pos'n/label indexing coincides

    @field_validator("supersystem_ie_only")
    @classmethod
    def set_short_circuit_mbe(cls, v: Optional[bool], info: FieldValidationInfo) -> bool:
        sio = v
        nfr = len(info.data["molecule"].fragments)

        if (sio is True) and (info.data["max_nbody"] != nfr):
            raise ValueError(f"Cannot skip intermediate n-body jobs when {max_nbody=} != nfragments={nfr}.")

        return sio

    @classmethod
    def from_qcschema(cls, input_model: ManyBodyInput):

        computer_model = cls(
            molecule=input_model.molecule,
            driver=input_model.specification.driver,
            **input_model.specification.keywords.model_dump(exclude_unset=True),
            input_data=input_model,  # storage, to reconstitute ManyBodyResult
        )

        return computer_model

    def build_tasks(
        self,
        mb_computer: AtomicComputerQCNG, #MBETaskComputers,
        mc_level_idx: int,
        **kwargs: Dict[str, Any],
    ) -> int:
        """Adds to the task_list as many new unique tasks as necessary to treat a single model chemistry level at one
        or several n-body levels. New tasks are of type *mb_computer* with model chemistry level specified in *kwargs*
        and n-body levels accessed through *mc_level_idx*.

        Parameters
        ----------
#        mb_computer
#            Class of task computers to instantiate and add to self.task_list. Usually :class:`~psi4.driver.AtomicComputer` but may be other when wrappers are layered.
#        mc_level_idx
#            Position in field self.nbodies_per_mc_level used to obtain ``nbodies``, the list of n-body
#            levels (e.g., `[1]` or `[1, 2]` or `["supersystem"]`) to which the modelchem specified in **kwargs** applies.
#            That is, `nbodies = self.nbodies_per_mc_level[mc_level_idx]`.
#            Note the natural 1-indexing of ``nbodies`` _contents_, so `[1]` covers one-body contributions.
#            The corresponding user label is the 1-indexed counterpart, `mc_level_lbl = mc_level_idx + 1`
#            Formerly nlevel as in `nbody = self.nbody_list[nbody_level=nlevel]`.
#        kwargs
#            Other arguments for initializing **mb_computer**. In particular, specifies model chemistry.

        Returns
        -------
        count : int
            Number of new tasks planned by this call.
            Formerly, didn't include supersystem in count.

        """
        from psi4.driver.driver_nbody import build_nbody_compute_list

        # TODO method not coming from levels right

        # Get the n-body orders for this level. e.g., [1] or [2, 3] or ["supersystem"]
        nbodies = self.nbodies_per_mc_level[mc_level_idx]

#        info = "\n" + p4util.banner(f" ManyBody Setup: N-Body Levels {nbodies}", strNotOutfile=True) + "\n"
#        core.print_out(info)
#        logger.info(info)

#        for kwg in ['dft_functional']:
#            if kwg in kwargs:
#                kwargs['keywords']['function_kwargs'][kwg] = kwargs.pop(kwg)

        count = 0
        template = copy.deepcopy(kwargs)

        # Get compute list
        if nbodies == ["supersystem"]:
            # Add supersystem computation if requested -- always nocp
            data = template
            data["molecule"] = self.molecule
            key = f"supersystem_{self.nfragments}"
            self.task_list[key] = mb_computer(**data)
            count += 1

            compute_dict = build_nbody_compute_list(
                ["nocp"], list(range(1, self.max_nbody + 1)),
                self.nfragments, self.return_total_data, self.supersystem_ie_only)
        else:
            compute_dict = build_nbody_compute_list(
                self.bsse_type, nbodies,
                self.nfragments, self.return_total_data, self.supersystem_ie_only)

        def labeler(item) -> str:
#            mc_level_lbl = mc_level_idx + 1
#            return str(mc_level_lbl) + "_" + str(item)
            # note 0-index to 1-index shift for label
            return f"{mc_level_idx + 1}_{item}"

        print("HHHH", compute_dict)
        # Add current compute list to the master task list
        # * `pair` looks like `((1,), (1, 3))` where first is real (not ghost) fragment indices
        #    and second is basis set fragment indices, all 1-indexed
        for nb in compute_dict["all"]:
            for pair in compute_dict["all"][nb]:
                lbl = labeler(pair)
                if lbl in self.task_list:
                    continue

                data = template
                ghost = list(set(pair[1]) - set(pair[0]))
                # while psi4.core.Molecule.extract_subsets takes 1-indexed, qcel.models.Molecule.get_fragment takes 0-indexed.
                real0 = [(idx - 1) for idx in list(pair[0])]
                ghost0 = [(idx - 1) for idx in ghost]
                data["molecule"] = self.molecule.get_fragment(real=real0, ghost=ghost0, group_fragments=False)  # orient?
#                if self.embedding_charges:
#                    embedding_frags = list(set(range(1, self.nfragments + 1)) - set(pair[1]))
#                    charges = []
#                    for frag in embedding_frags:
#                        positions = self.molecule.extract_subsets(frag).geometry().np.tolist()
#                        charges.extend([[chg, i] for i, chg in zip(positions, self.embedding_charges[frag])])
#                    data['keywords']['function_kwargs'].update({'external_potentials': charges})

                self.task_list[lbl] = mb_computer(**data)
                count += 1

        return count

    def plan(self):
        # uncalled function
        return [t.plan() for t in self.task_list.values()]

    def compute(self):
        """Run quantum chemistry.

        NOTE: client removed (compared to psi4.driver.ManyBodyComputer)
        """
        from psi4.driver.p4util import banner

        info = "\n" + banner(f" ManyBody Computations ", strNotOutfile=True) + "\n"
        #logger.info(info)

        for t in self.task_list.values():
            t.compute()

    def prepare_results(
        self,
        results: Optional[Dict[str, "MBETaskComputers"]] = None,
    ) -> Dict[str, Any]:
        """Process the results from all n-body component molecular systems and model chemistry levels into final quantities.

        NOTE: client removed (compared to psi4.driver.ManyBodyComputer)

#        Parameters
#        ----------
#        results
#            A set of tasks to process instead of self.task_list. Used in multilevel processing to pass a subset of
#            self.task_list filtered to only one modelchem level.
#        client
#            QCFractal client if using QCArchive for distributed compute.
#
#        Returns
#        -------
#        nbody_results
#            When the ManyBodyComputer specifies a single model chemistry level (see self.nbodies_per_mc_level), the
#            return is a dictionary, nbody_results, described in the table below. Many of the items are actually filled
#            by successive calls to assemble_nbody_components(). When multiple model chemistry levels are specified, this
#            function diverts its return to driver_nbody_multilevel.prepare_results() wherein each mc level calls this
#            function again and collects separate nbody_results dictionaries and processes them into a final return that
#            is a small subset of the table below.
#
#
#                                       ptype_size = (1,)/(nat, 3)/(3 * nat, 3 * nat)
#                                        e/g/h := energy or gradient or Hessian
#                                        rtd := return_total_data
#
        """
        from psi4.driver.driver_nbody import assemble_nbody_components

        if results is None:
            results = {}

#        # formerly nlevels
        mc_level_labels = {i.split("_")[0] for i in self.task_list}
        if len(mc_level_labels) > 1 and not results:
            return psi4.driver.driver_nbody_multilevel.prepare_results(self, client)

        results_list = {k: v.get_results() for k, v in (results.items() or self.task_list.items())}
        trove = {  # AtomicResult.properties return None if missing
            "energy": {k: v.properties.return_energy for k, v in results_list.items()},
            "gradient": {k: v.properties.return_gradient for k, v in results_list.items()},
            "hessian": {k: v.properties.return_hessian for k, v in results_list.items()},
        }

#        # TODO: make assemble_nbody_components and driver_nbody_multilevel.prepare_results into class functions.
#        #   note that the former uses metadata as read-only (except for one solveable case) while the latter overwrites self (!).
        metadata = {
            "quiet": False, #self.quiet,
            "nbodies_per_mc_level": self.nbodies_per_mc_level,
            "bsse_type": self.bsse_type,
            "nfragments": self.nfragments,
            "return_total_data": self.return_total_data,
            "supersystem_ie_only": self.supersystem_ie_only,
            "molecule": self.molecule,
            "embedding_charges": bool(self.embedding_charges),
            "max_nbody": self.max_nbody,
        }
        if self.driver.name == "energy":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())

        elif self.driver.name == "gradient":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())
            nbody_results.update(assemble_nbody_components("gradient", trove["gradient"], metadata.copy()))

        elif self.driver.name == "hessian":
            nbody_results = assemble_nbody_components("energy", trove["energy"], metadata.copy())
            nbody_results.update(assemble_nbody_components("gradient", trove["gradient"], metadata.copy()))
            nbody_results.update(assemble_nbody_components("hessian", trove["hessian"], metadata.copy()))

        # save some mc_(frag, bas) component results
        # * formerly, intermediates_energy was intermediates2
        # * formerly, intermediates_gradient was intermediates_ptype
        # * formerly, intermediates_hessian was intermediates_ptype

        nbody_results["intermediates"] = {}
        for idx, task in results_list.items():
            mc, frag, bas = delabeler(idx)
            nbody_results["intermediates"][f"N-BODY ({frag})@({bas}) TOTAL ENERGY"] = task.properties.return_energy

        nbody_results["intermediates_energy"] = trove["energy"]

        if not all(x is None for x in trove["gradient"].values()):
            nbody_results["intermediates_gradient"] = trove["gradient"]

        if not all(x is None for x in trove["hessian"].values()):
            nbody_results["intermediates_hessian"] = trove["hessian"]

#        debug = False
#        if debug:
#            for k, v in nbody_results.items():
#                if isinstance(v, np.ndarray):
#                    print(f"CLS-prepared results >>> {k} {v.size}")
#                elif isinstance(v, dict):
#                    print(f"CLS-prepared results >>> {k} {len(v)}")
#                    for k2, v2 in v.items():
#                        if isinstance(v2, np.ndarray):
#                            print(f"CLS-prepared results      >>> {k2} {v2.size}")
#                        else:
#                            print(f"CLS-prepared results      >>> {k2} {v2}")
#                else:
#                    print(f"CLS-prepared results >>> {k} {v}")

        return nbody_results


    def get_results(self) -> AtomicResult:
        """Return results as ManyBody-flavored QCSchema.

        NOTE: client removed (compared to psi4.driver.ManyBodyComputer)
        """
        from psi4.driver.p4util import banner

        info = "\n" + banner(f" ManyBody Results ", strNotOutfile=True) + "\n"
        #logger.info(info)

        results = self.prepare_results()
        ret_energy = results.pop("ret_energy")
        ret_ptype = results.pop("ret_ptype")
        ret_gradient = results.pop("ret_gradient", None)

        # load QCVariables
        qcvars = {
            'NUCLEAR REPULSION ENERGY': self.molecule.nuclear_repulsion_energy(),
            'NBODY NUMBER': len(self.task_list),
        }

        properties = {
            "calcinfo_nmc": len(self.nbodies_per_mc_level),
            "calcinfo_nfr": self.nfragments,  # or len(self.molecule.fragments)
            "calcinfo_natom": len(self.molecule.symbols),
            "calcinfo_nmbe": len(self.task_list),
            "nuclear_repulsion_energy": self.molecule.nuclear_repulsion_energy(),
            "return_energy": ret_energy,
        }

        for k, val in results.items():
            qcvars[k] = val

        qcvars['CURRENT ENERGY'] = ret_energy
        if self.driver == 'gradient':
            qcvars['CURRENT GRADIENT'] = ret_ptype
            properties["return_gradient"] = ret_ptype
        elif self.driver == 'hessian':
            qcvars['CURRENT GRADIENT'] = ret_gradient
            qcvars['CURRENT HESSIAN'] = ret_ptype
            properties["return_gradient"] = ret_gradient
            properties["return_hessian"] = ret_ptype

        atprop = build_manybodyproperties(qcvars["nbody"])
#        print("ATPROP")
#        pp.pprint(atprop.model_dump())

        for qcv, val in qcvars.items():
            if isinstance(val, dict):
                if qcv in [
                ]:
                    for qcv2, val2 in val.items():
                            qcvars[str(qcv2)] = val2
            else:
                qcvars[qcv] = val

        component_results = self.model_dump()['task_list']  # TODO when/where include the indiv outputs
#        for k, val in component_results.items():
#            val['molecule'] = val['molecule'].to_schema(dtype=2)

        nbody_model = ManyBodyResult(
            **{
                'input_data': self.input_data,
                'properties': {**atprop.model_dump(), **properties},
                'provenance': provenance_stamp(__name__),
                'extras': {
                    'qcvars': qcvars,
#                    'component_results': component_results,
                },
                'return_result': ret_ptype,
                'success': True,
            })

#        logger.debug('\nNBODY QCSchema:\n' + pp.pformat(nbody_model.model_dump()))

        return nbody_model





def delabeler(item: str, return_obj: bool = False) -> Union[Tuple[str, str, str], Tuple[int, Tuple[int], Tuple[int]]]:
    """Transform labels like string "1_((2,), (1, 2))" into string tuple ("1", "2", "1, 2") or
    object tuple (1, (2,), (1, 2)).

    """
    mc, _, fragbas = item.partition("_")
    frag, bas = literal_eval(fragbas)

    if return_obj:
        return int(mc), frag, bas
    else:
        return mc, ", ".join(map(str, frag)), ", ".join(map(str, bas))


qcvars_to_manybodyproperties = {}
for skprop in ManyBodyResultProperties.model_fields.keys():
    qcvar = skprop.replace("_body", "-body").replace("_corr", "-corr").replace("_", " ").upper()
    qcvars_to_manybodyproperties[qcvar] = skprop
qcvars_to_manybodyproperties["CURRENT ENERGY"] = "return_energy"
qcvars_to_manybodyproperties["CURRENT GRADIENT"] = "return_gradient"
qcvars_to_manybodyproperties["CURRENT HESSIAN"] = "return_hessian"


def build_manybodyproperties(qcvars: Dict) -> ManyBodyResultProperties:
    """For results extracted from QC output in QCDB terminology, translate to QCSchema terminology.

    Parameters
    ----------
    qcvars : PreservingDict
        Dictionary of calculation information in QCDB QCVariable terminology.

    Returns
    -------
    atprop : ManyBodyResultProperties
        Object of calculation information in QCSchema ManyBodyResultProperties terminology.

    """
    atprop = {}
    for pv, dpv in qcvars.items():
        if pv in qcvars_to_manybodyproperties:
            atprop[qcvars_to_manybodyproperties[pv]] = dpv

    return ManyBodyResultProperties(**atprop)
