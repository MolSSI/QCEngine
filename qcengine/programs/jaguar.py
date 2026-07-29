"""QCEngine harness for Jaguar.

The Schrödinger Python modules are intentionally imported at run time so that
QCEngine itself does not depend on the Schrödinger suite.
"""

import os
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Dict, Optional

import numpy as np
from qcelemental import constants
from qcelemental.models.v2 import AtomicResult, BasisSet, Provenance
from qcelemental.util import safe_version, which_import

from ..exceptions import InputError, UnknownError
from ..util import environ_context, temporary_directory
from .model import ProgramHarness

if TYPE_CHECKING:
    from qcelemental.models.v2 import AtomicInput

    from ..config import TaskConfig


class JaguarHarness(ProgramHarness):
    """Run Jaguar calculations through the Schrödinger Python API.

    Jaguar is file based, so the harness uses a dedicated scratch directory.
    The current launch path temporarily changes process-wide state and is
    therefore not thread-safe. A single calculation can use multiple threads
    through Jaguar's ``-PARALLEL`` option, but the harness does not configure
    multi-node execution. Jaguar is treated as a memory-managed quantum
    chemistry program.
    """

    _defaults: ClassVar[Dict[str, Any]] = {
        "name": "Jaguar",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Check whether the Schrödinger Jaguar Python API is importable.

        Parameters
        ----------
        raise_error
            Raise an informative import error instead of returning ``False``
            when the API is unavailable.

        Returns
        -------
        bool
            Whether the Jaguar input module can be imported.
        """
        return which_import(
            "schrodinger.application.jaguar.input",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Run QCEngine from a Schrödinger Python environment with SCHRODINGER defined.",
        )

    def get_version(self) -> str:
        """Return the normalized Jaguar product version.

        Schrödinger stores this as an integer such as ``134039``, denoting
        Jaguar 13.4 release 39. The normalized value is cached by the path of
        the module that supplies the product version.

        Returns
        -------
        str
            A PEP 440-compatible Jaguar version such as ``13.4.39``.
        """
        self.found(raise_error=True)

        module_path = which_import("schrodinger.infra.mm")
        cache_key = f"jaguar:{module_path}"
        if cache_key not in self.version_cache:
            from schrodinger.infra import mm

            encoded_version = int(mm.mmfile_get_product_version("jaguar"))
            major_minor, release = divmod(encoded_version, 1000)
            major, minor = divmod(major_minor, 10)
            self.version_cache[cache_key] = safe_version(f"{major}.{minor}.{release:03d}")
        return self.version_cache[cache_key]

    def get_suite_version(self) -> str:
        """Return the normalized Schrödinger suite release.

        Returns
        -------
        str
            A PEP 440-compatible suite release such as ``2026.post4``.
        """
        self.found(raise_error=True)

        module_path = which_import("schrodinger.application.jaguar.constants")
        cache_key = f"suite:{module_path}"
        if cache_key not in self.version_cache:
            from schrodinger.application.jaguar.constants import VERSION

            self.version_cache[cache_key] = safe_version(VERSION)
        return self.version_cache[cache_key]

    @staticmethod
    def _validate_input(input_model: "AtomicInput") -> None:
        """Validate that an atomic input uses features supported by Jaguar.

        Parameters
        ----------
        input_model
            The QCSchema atomic input to validate.

        Raises
        ------
        InputError
            If the driver is unsupported, the basis is an explicit QCSchema
            ``BasisSet``, or a ghost dummy center is requested.
        """
        driver = input_model.specification.driver
        if driver not in {"energy", "gradient", "hessian"}:
            raise InputError(f"Driver {driver} not implemented for Jaguar.")

        basis = input_model.specification.model.basis
        if isinstance(basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use a string basis name.")

        for symbol, real in zip(input_model.molecule.symbols, input_model.molecule.real):
            if str(symbol).lower() == "x" and not real:
                raise InputError("A QCSchema dummy center (symbol X) cannot also be a ghost atom.")

    @staticmethod
    def _build_jaguar_input(input_model: "AtomicInput", job_name: str):
        """Translate a QCSchema atomic input into a Jaguar input object.

        Parameters
        ----------
        input_model
            The validated atomic input containing the molecule, model,
            driver, and native Jaguar keywords.
        job_name
            Basename or path to use for the Jaguar job and input file.

        Returns
        -------
        schrodinger.application.jaguar.input.JaguarInput
            A configured Jaguar input ready to be saved or run.

        Notes
        -----
        QCSchema coordinates are converted from bohr to angstrom. Model and
        driver settings override conflicting values in the native keywords.
        QCSchema dummy atoms (symbol ``X``) are represented as Schrödinger
        ``Du`` atoms, which Jaguar serializes with its standard ``X<n>`` label.
        Ghost atoms retain their element and are marked as counterpoise atoms,
        which Jaguar serializes by appending ``@`` to the atom label.
        """
        from schrodinger import structure
        from schrodinger.application.jaguar.input import JaguarInput

        molecule = input_model.molecule
        geometry = np.asarray(molecule.geometry) * constants.bohr2angstroms
        jaguar_structure = structure.create_new_structure()
        for symbol, xyz in zip(molecule.symbols, geometry):
            jaguar_symbol = "Du" if str(symbol).lower() == "x" else str(symbol)
            jaguar_structure.addAtom(jaguar_symbol, *map(float, xyz))

        keywords = {key.lower(): value for key, value in input_model.specification.keywords.items()}
        model = input_model.specification.model

        # Model and driver fields are authoritative over duplicate native keywords.
        keywords["dftname"] = "HF" if model.method.lower() in {"hf", "scf"} else model.method
        if model.basis:
            keywords["basis"] = model.basis

        driver = input_model.specification.driver
        keywords["igeopt"] = -1 if driver == "gradient" else 0
        keywords["ifreq"] = 1 if driver == "hessian" else 0
        keywords["isymm"] = 0

        jaguar_input = JaguarInput(name=job_name, structure=jaguar_structure, genkeys=keywords)
        for atom_index, real in enumerate(molecule.real, start=1):
            if not real:
                jaguar_input._setCounterpoise(atom_index, True)
        jaguar_input.setValue("molchg", int(molecule.molecular_charge))
        jaguar_input.setValue("multip", molecule.molecular_multiplicity)
        return jaguar_input

    @staticmethod
    def _optional_properties(jaguar_output: Any, input_model: "AtomicInput") -> Dict[str, Any]:
        """Translate available Jaguar results into QCSchema properties.

        Parameters
        ----------
        jaguar_output
            Parsed ``JaguarOutput`` returned by the Jaguar calculation.
        input_model
            The atomic input, used for molecule-derived metadata.

        Returns
        -------
        dict
            QCSchema ``AtomicResultProperties`` fields available in the
            Jaguar output. Missing optional quantities are omitted.
        """
        results = jaguar_output.last_results
        properties: Dict[str, Any] = {
            "calcinfo_natom": len(input_model.molecule.symbols),
            "return_energy": results.energy,
        }

        direct_mapping = {
            "nbasis": "calcinfo_nbasis",
            "nuclear_repulsion": "nuclear_repulsion_energy",
            "scf_energy": "scf_total_energy",
            "energy_one_electron": "scf_one_electron_energy",
            "energy_two_electron": "scf_two_electron_energy",
            "rimp2_ss_energy": "mp2_same_spin_correlation_energy",
            "rimp2_os_energy": "mp2_opposite_spin_correlation_energy",
            "rimp2_corr_energy": "mp2_correlation_energy",
            "rimp2_energy": "mp2_total_energy",
        }
        for source, target in direct_mapping.items():
            owner = jaguar_output if source == "nbasis" else results
            value = getattr(owner, source, None)
            if value is not None:
                properties[target] = value

        nalpha = getattr(jaguar_output, "num_occ_orbs_alpha", None)
        nbeta = getattr(jaguar_output, "num_occ_orbs_beta", None)
        if nalpha is None and nbeta is None:
            nocc = getattr(jaguar_output, "num_occ_orbs", None)
            nalpha = nbeta = nocc
        if nalpha is not None:
            properties["calcinfo_nalpha"] = nalpha
        if nbeta is not None:
            properties["calcinfo_nbeta"] = nbeta

        dipole = getattr(results, "dipole_qm", None)
        if dipole is not None and None not in (dipole.x, dipole.y, dipole.z):
            properties["scf_dipole_moment"] = (
                np.asarray([dipole.x, dipole.y, dipole.z]) / constants.dipmom_au2debye
            )

        return properties

    @staticmethod
    def _read_hessian(job_base: Path, natom: int):
        """Read a Cartesian Hessian from a generated Jaguar input file.

        Parameters
        ----------
        job_base
            Path to the Jaguar job without a filename extension.
        natom
            Number of atoms, used to validate the Hessian dimensions.

        Returns
        -------
        tuple[numpy.ndarray, pathlib.Path]
            The ``(3 * natom, 3 * natom)`` Hessian and the generated input
            file containing its ``&hess`` section.

        Raises
        ------
        UnknownError
            If no Hessian is found or its dimensions are inconsistent.
        """
        from schrodinger.application.jaguar.input import JaguarInput

        expected_shape = (3 * natom, 3 * natom)
        candidates = [job_base.with_suffix(".01.in"), *sorted(job_base.parent.glob(f"{job_base.name}*.01.in"))]
        for candidate in dict.fromkeys(candidates):
            if not candidate.is_file():
                continue
            hessian = JaguarInput(str(candidate), compute_connectivity=False).getHessian()
            if hessian is not None:
                hessian = np.asarray(hessian)
                if hessian.shape != expected_shape:
                    raise UnknownError(
                        f"Jaguar Hessian in {candidate.name} has shape {hessian.shape}; expected {expected_shape}."
                    )
                return hessian, candidate

        raise UnknownError("Jaguar completed the Hessian calculation but no &hess section was found.")

    @staticmethod
    def _collect_native_files(job_base: Path, restart_file: Optional[Path]) -> Dict[str, str]:
        """Read selected Jaguar files before the scratch directory is removed.

        Parameters
        ----------
        job_base
            Path to the Jaguar job without a filename extension.
        restart_file
            Generated input known to contain the requested Hessian, if any.

        Returns
        -------
        dict
            Text contents of the original input, generated numbered inputs,
            and Jaguar log, keyed by their native filenames. QCSchema output
            protocols determine which entries survive in the final result.
        """
        native_files = {}
        input_file = job_base.with_suffix(".in")
        if input_file.is_file():
            native_files["input"] = input_file.read_text(errors="replace")

        generated_inputs = [restart_file, *sorted(job_base.parent.glob(f"{job_base.name}.*.in"))]
        for generated_input in dict.fromkeys(generated_inputs):
            if generated_input is not None and generated_input.is_file():
                native_files[generated_input.name] = generated_input.read_text(errors="replace")

        log_file = job_base.with_suffix(".log")
        if log_file.is_file():
            native_files[log_file.name] = log_file.read_text(errors="replace")
        return native_files

    def compute(self, input_model: "AtomicInput", config: "TaskConfig") -> AtomicResult:
        """Execute Jaguar and construct a QCSchema atomic result.

        Parameters
        ----------
        input_model
            Atomic input specifying the molecule, model, driver, keywords,
            and output protocols.
        config
            QCEngine task resources and scratch-directory configuration.

        Returns
        -------
        AtomicResult
            Validated QCSchema result containing the requested energy,
            gradient, or Hessian and available auxiliary properties.

        Raises
        ------
        InputError
            If the requested input feature is unsupported.
        UnknownError
            If Jaguar fails, omits the requested result, or produces an
            inconsistent Hessian.
        """
        self.found(raise_error=True)
        self._validate_input(input_model)

        with temporary_directory(
            parent=config.scratch_directory, suffix="_jaguar_scratch", messy=config.scratch_messy
        ) as tmpdir:
            job_base = Path(tmpdir) / "dispatch"
            jaguar_input = self._build_jaguar_input(input_model, str(job_base))

            previous_directory = Path.cwd()
            try:
                # Job Control returns result files to the launch directory.
                os.chdir(tmpdir)
                with environ_context(config=config):
                    jaguar_output = jaguar_input.run(**{"-PARALLEL": str(config.ncores)})
            except Exception as exc:
                raise UnknownError(str(exc)) from exc
            finally:
                os.chdir(previous_directory)

            if jaguar_output.status != jaguar_output.OK:
                raise UnknownError(
                    f"Jaguar calculation failed with ERROR {jaguar_output.fatal_errorno}: "
                    f"{jaguar_output.fatal_error}"
                )

            properties = self._optional_properties(jaguar_output, input_model)
            results = jaguar_output.last_results
            driver = input_model.specification.driver
            restart_file = None

            if driver == "energy":
                return_result = results.energy
            elif driver == "gradient":
                if results.forces is None:
                    raise UnknownError("Jaguar completed the gradient calculation but returned no forces.")
                return_result = -np.asarray(results.forces)
                properties["return_gradient"] = return_result
            else:
                return_result, restart_file = self._read_hessian(job_base, len(input_model.molecule.symbols))
                properties["return_hessian"] = return_result

            output_file = job_base.with_suffix(".out")
            stdout = output_file.read_text(errors="replace") if output_file.is_file() else None
            native_files = self._collect_native_files(job_base, restart_file)

        provenance = Provenance(
            creator="Jaguar",
            version=self.get_version(),
            routine="schrodinger.application.jaguar.input.JaguarInput.run",
            nthreads=config.ncores,
        )
        return AtomicResult(
            input_data=input_model,
            molecule=input_model.molecule,
            properties=properties,
            return_result=return_result,
            stdout=stdout,
            native_files=native_files,
            success=True,
            provenance=provenance,
            extras={
                "jaguar": {
                    "jaguar_version": self.get_version(),
                    "suite_version": self.get_suite_version(),
                    "method": jaguar_output.method,
                    "functional": jaguar_output.functional,
                    "basis": jaguar_output.basis,
                    "point_group": jaguar_output.point_group,
                }
            },
        )
