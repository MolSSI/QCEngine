"""Compute quantum chemistry using the Valeev group's MPQC4 executable."""

import json
import os
import pprint
from typing import Any, ClassVar, Dict, Optional

from qcelemental.models.v2 import AtomicInput, AtomicResult, BasisSet
from qcelemental.util import safe_version, which

from ...config import TaskConfig
from ...exceptions import InputError, UnknownError
from ...util import create_mpi_invocation, execute
from ..model import ProgramHarness
from .germinate import EXCITATION_ENERGY, muster_modelchem
from .keywords import deep_merge, extract_reserved, format_keywords

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)

#: Default fundamental-constants set, pinned fallback so results are reproducible across
#: MPQC versions that change their own default.
_UNITS_SYSTEM = "2018CODATA"

#: MPQC's own ExcitationEnergy default.
_DEFAULT_N_ROOTS = 3


def _molecule_block(molecule) -> Dict[str, Any]:
    """Build MPQC's ``atoms`` block from a QCSchema Molecule.

    Geometry is emitted in Bohr, matching ``Molecule.geometry``, with
    ``units: bohr`` set explicitly because the ``atoms`` block defaults to
    angstrom. ``sort_input`` is pinned False so MPQC preserves atom order.

    Built directly rather than via ``qcel.molparse.to_string(dtype=...)``:
    that helper targets text formats, and MPQC's molecule spec is JSON.
    """
    geometry = molecule.geometry.reshape(-1, 3)
    return {
        "units": "bohr",
        "sort_input": False,
        "atoms": [
            {"element": symbol, "xyz": [float(x) for x in xyz]} for symbol, xyz in zip(molecule.symbols, geometry)
        ],
    }


class MPQCHarness(ProgramHarness):
    """Interface for the MPQC4 project.

    Notes
    -----
    * Input is a JSON KeyVal document. Results are read from the
      ``Output KeyVal (format=JSON)`` block MPQC writes to stdout.
    * Only ``Energy`` and ``ExcitationEnergy`` properties are supported;
      ``gradient`` and ``hessian`` drivers raise ``InputError``.
    * Basis sets resolve from each directory in ``MPQC_BASIS_PATH`` first
      (colon-separated, in listed order), and only then from libint2's
      bundled library, whose location ``LIBINT_DATA_PATH`` sets. So
      ``MPQC_BASIS_PATH`` shadows ``LIBINT_DATA_PATH`` rather than replacing
      it -- the libint2 library stays the last-resort fallback. Both are
      inherited from the parent environment; see ``_build_environment``.
    """

    _defaults: ClassVar[Dict[str, Any]] = {
        "name": "MPQC",
        "scratch": True,
        "thread_safe": False,
        "thread_parallel": True,
        "node_parallel": True,
        "managed_memory": True,
    }
    version_cache: Dict[str, str] = {}

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        return which(
            "mpqc",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install via https://github.com/ValeevGroup/mpqc4",
        )

    def get_version(self) -> str:
        self.found(raise_error=True)

        which_prog = which("mpqc")
        if which_prog not in self.version_cache:
            # `mpqc -v` prints e.g. "MPQC version 4.0.0-beta.1" and exits 0
            # without needing an input file.
            success, output = execute([which_prog, "-v"], {})
            if not success:
                raise UnknownError(output["stderr"])

            version = None
            for line in output["stdout"].splitlines():
                if line.startswith("MPQC version"):
                    version = line.split()[2]
                    break
            if version is None:
                raise UnknownError(f"Could not parse MPQC version from: {output['stdout']}")

            self.version_cache[which_prog] = safe_version(version)

        return self.version_cache[which_prog]

    def compute(self, input_model: AtomicInput, config: TaskConfig) -> AtomicResult:
        raise NotImplementedError("The MPQC harness cannot yet run calculations.")

    def build_input(
        self, input_model: AtomicInput, config: TaskConfig, template: Optional[str] = None
    ) -> Dict[str, Any]:
        spec = input_model.specification
        molecule = input_model.molecule

        # Guards, all before any subprocess launches.
        if isinstance(spec.model.basis, BasisSet):
            raise InputError("QCSchema BasisSet for model.basis not implemented. Use string basis name.")
        if spec.model.basis is None:
            raise InputError("None for model.basis is not useable. MPQC requires a basis set name.")
        if not all(molecule.real):
            raise InputError("Ghost atoms are not supported by the MPQC harness; MPQC has no ghost-atom facility.")

        opts, mpqc_input, mpqc_env, property_opts = extract_reserved(dict(spec.keywords))
        user_tree = format_keywords(opts)

        # Raises InputError for bad methods, derivative drivers, and a
        # property__type the wfn type cannot provide. muster_modelchem unwraps
        # DriverEnum itself. `type` is popped rather than read: the returned
        # property_type is authoritative and is written into the block below,
        # in MPQC's exact casing regardless of how the user spelled it.
        wfn_type, needs_ref, property_type, wfn_extras = muster_modelchem(
            spec.model.method, spec.driver, property_opts.pop("type", None)
        )

        if mpqc_input:
            tree = deep_merge(mpqc_input, user_tree)
        else:
            tree = deep_merge(self._skeleton(molecule, spec, wfn_type, needs_ref, wfn_extras), user_tree)

        # The molecule and the property block are always harness-owned, even
        # under the mpqc_input escape hatch.
        tree["atoms"] = _molecule_block(molecule)
        tree.setdefault("units", _UNITS_SYSTEM)

        property_block: Dict[str, Any] = {
            "type": property_type,
            "wfn": "$:wfn",
            "precision": 1.0e-10,
        }
        if property_type == EXCITATION_ENERGY:
            property_block["n_roots"] = _DEFAULT_N_ROOTS
        property_block.update(property_opts)
        # The property block is harness-owned, but the rest of the `mpqc`
        # subtree is not: `mpqc:file_prefix` and `mpqc:tasks` are real MPQC keys,
        # so a user's mpqc__* keyword must survive. Merging
        # rather than assigning keeps them.
        tree["mpqc"] = {**tree.get("mpqc", {}), "property": property_block}

        # A df method leaves a dangling $:dfbs unless the user supplied one.
        if tree.get("wfn", {}).get("method") == "df" or "DF" in wfn_type.upper():
            if "dfbs" not in tree or "name" not in tree.get("dfbs", {}):
                raise InputError(
                    f"MPQC method/type '{wfn_type}' uses density fitting but no density-fitting "
                    f"basis was given. Set the dfbs__name keyword, e.g. dfbs__name='cc-pVDZ-RI'."
                )
            # setdefault, not tree["wfn_world"][...]: under the mpqc_input
            # escape hatch there may be no wfn_world block to index into.
            tree.setdefault("wfn_world", {})["df_basis"] = "$:dfbs"

        if config.use_mpiexec:
            command = create_mpi_invocation(which("mpqc"), config)
        else:
            command = [which("mpqc")]
        command += ["-i", "mpqc.json"]

        return {
            "infiles": {"mpqc.json": json.dumps(tree, indent=2)},
            "command": command,
            "scratch_directory": config.scratch_directory,
            "scratch_messy": config.scratch_messy,
            "environment": os.environ.copy(),
        }

    @staticmethod
    def _skeleton(molecule, spec, wfn_type: str, needs_ref: bool, wfn_extras: Dict[str, Any]) -> Dict[str, Any]:
        """Generate the default MPQC block tree for a single-point calculation."""
        multiplicity = molecule.molecular_multiplicity
        scf_block = {
            "type": "SD",
            "wfn_world": "$:wfn_world",
            "atoms": "$:atoms",
            "fock": {"type": "DirectFockBuilder", "wfn_world": "$:wfn_world"},
            "charge": int(molecule.molecular_charge),
            "multiplicity": multiplicity,
            "spin_restricted": multiplicity == 1,
            "max_iter": 100,
        }

        tree: Dict[str, Any] = {
            "units": _UNITS_SYSTEM,
            "obs": {"name": spec.model.basis, "atoms": "$:atoms"},
            "wfn_world": {"atoms": "$:atoms", "basis": "$:obs"},
        }

        wfn_block: Dict[str, Any] = {
            "type": wfn_type,
            "wfn_world": "$:wfn_world",
            "atoms": "$:atoms",
        }
        if needs_ref:
            tree["scf"] = scf_block
            wfn_block["ref"] = "$:scf"
        else:
            # The wfn is the SCF; fold the reference settings into it.
            wfn_block.update({k: v for k, v in scf_block.items() if k not in ("type", "wfn_world", "atoms")})
        wfn_block.update(wfn_extras)
        tree["wfn"] = wfn_block
        return tree
