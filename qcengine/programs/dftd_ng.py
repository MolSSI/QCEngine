"""
Harness for the DFT-D dispersion correction.
This implementation interfaces with the dftd3 and dftd4 Python-API, which provides
native support for QCSchema.

Therefore, this harness only has to provide a thin wrapper to integrate the
respective dispersion correction.
"""

from typing import Dict

from qcelemental.models import AtomicInput, AtomicResult
from qcelemental.util import parse_version, safe_version, which_import

from ..config import TaskConfig
from ..exceptions import InputError
from .empirical_dispersion_resources import from_arrays, get_dispersion_aliases
from .model import ProgramHarness


class DFTD4Harness(ProgramHarness):
    """Calculation harness for the DFT-D4 dispersion correction."""

    _defaults = {
        "name": "dftd4",
        "scratch": False,
        "thread_safe": True,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": False,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Check for the availability of the Python API of dftd4"""

        return which_import(
            "dftd4",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install a dftd4 version with enabled Python API"
            + " (e.g. conda install dftd4-python -c conda-forge)",
        )

    def get_version(self) -> str:
        """Return the currently used version of dftd4"""
        self.found(raise_error=True)

        which_prog = which_import("dftd4")
        if which_prog not in self.version_cache:
            import dftd4

            self.version_cache[which_prog] = safe_version(dftd4.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_model: AtomicInput, config: TaskConfig) -> AtomicResult:
        """
        Actual interface to the dftd4 package. The compute function is just a thin
        wrapper around the native QCSchema interface of the dftd4 Python-API.
        """

        self.found(raise_error=True)

        import dftd4
        from dftd4.qcschema import run_qcschema

        # strip engine hint
        input_data = input_model.dict()
        method = input_model.model.method
        if method.startswith("d4-"):
            method = method[3:]
            input_data["model"]["method"] = method
        qcvkey = method.upper() if method is not None else None

        # send `from_arrays` the dftd4 behavior of functional specification overrides explicit parameters specification
        # * differs from dftd4 harness behavior where parameters extend or override functional
        # * stash the resolved plan in extras or, if errored, leave it for the proper dftd4 api to reject
        param_tweaks = None if method else input_model.keywords.get("params_tweaks", None)
        try:
            planinfo = from_arrays(
                verbose=1,
                name_hint=method,
                level_hint=input_model.keywords.get("level_hint", None),
                param_tweaks=param_tweaks,
                dashcoeff_supplement=input_model.keywords.get("dashcoeff_supplement", None),
            )
        except InputError:
            pass
        else:
            input_data["extras"]["info"] = planinfo

        # strip dispersion level from method
        for alias, d4 in get_dispersion_aliases().items():
            if d4 == "d4bjeeqatm" and method.lower().endswith(alias):
                method = method[: -(len(alias) + 1)]
                input_data["model"]["method"] = method

        # consolidate dispersion level aliases
        level_hint = input_model.keywords.get("level_hint", None)
        if level_hint and get_dispersion_aliases()[level_hint.lower()] == "d4bjeeqatm":
            level_hint = "d4"
            input_data["keywords"]["level_hint"] = level_hint

        input_model = AtomicInput(**input_data)

        # Run the Harness
        output = run_qcschema(input_model)

        if "info" in output.extras:
            qcvkey = output.extras["info"]["fctldash"].upper()

        calcinfo = {}
        energy = output.properties.return_energy
        calcinfo["CURRENT ENERGY"] = energy
        calcinfo["DISPERSION CORRECTION ENERGY"] = energy
        if qcvkey:
            calcinfo[f"{qcvkey} DISPERSION CORRECTION ENERGY"] = energy

        if output.driver == "gradient":
            gradient = output.return_result
            calcinfo["CURRENT GRADIENT"] = gradient
            calcinfo["DISPERSION CORRECTION GRADIENT"] = gradient
            if qcvkey:
                calcinfo[f"{qcvkey} DISPERSION CORRECTION GRADIENT"] = gradient

        if output.keywords.get("pair_resolved", False):
            pw2 = output.extras["dftd4"]["additive pairwise energy"]
            pw3 = output.extras["dftd4"]["non-additive pairwise energy"]
            assert abs(pw2.sum() + pw3.sum() - energy) < 1.0e-8, f"{pw2.sum()} + {pw3.sum()} != {energy}"
            calcinfo["2-BODY DISPERSION CORRECTION ENERGY"] = pw2.sum()
            calcinfo["3-BODY DISPERSION CORRECTION ENERGY"] = pw3.sum()
            calcinfo["2-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS"] = pw2
            calcinfo["3-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS"] = pw3

        output.extras["qcvars"] = calcinfo

        return output


class SDFTD3Harness(ProgramHarness):
    """
    Calculation harness for the DFT-D3 dispersion correction.

    This implementation of DFT-D3 supports the several damping functions, which
    are selected via the *level_hint* keyword. Damping parameter can be specified
    via the *param_tweaks* dictionary. If no *param_tweaks* are provided the
    functional parameters are obtained from the internal database of the library.

    The following damping function are available via *level_hint*:

    - ``d3bj``:
      Rational damping function for DFT-D3. The original scheme was proposed by
      Becke and Johnson and implemented in a slightly adjusted form using only
      the C8/C6 ratio in the critical radius for DFT-D3.
      Requires at least three parameters: *s8*, *a1*, and *a2*.
      The parameters *s6*, *s9*, and *alpha6* can be adjusted as well.
    - ``d3zero``:
      Original DFT-D3 damping function, based on a variant proposed by Chai and Head-Gordon.
      Requires at least two parameters: *s8* and *sr6*.
      The parameters *s6*, *s9*, *sr8*, and *alpha6* can be adjusted as well.
    - ``d3mbj``:
      Modified version of the rational damping parameters. The functional form of the
      damping function is *unmodified* with respect to the original rational damping scheme.
      However, for a number of functionals new parameters were introduced.
      Requires at least three parameters: *s8*, *a1*, and *a2*.
      The parameters *s6*, *s9*, and *alpha6* can be adjusted as well.
    - ``d3mzero``:
      Modified zero damping function for DFT-D3. This scheme adds an additional offset
      parameter to the zero damping scheme of the original DFT-D3.
      Requires at least three parameters: *s8*, *sr6*, and *beta*.
      The parameters *s6*, *s9*, *sr8*, and *alpha6* can be adjusted as well.
    - ``d3op``:
      Optimized power version of the rational damping function for DFT-D3.
      The functional form of the damping function is modified by adding an additional
      zero-damping like power function.
      Requires at least four parameters: *s8*, *a1*, *a2*, and *beta*.
      The parameters *s6*, *s9*, and *alpha6* can be adjusted as well.

    All damping functions by default *include* the ATM three-body contributions,
    it must be explicitly disabled by setting the *s9* value to zero.
    """

    _defaults = {
        "name": "s-dftd3",
        "scratch": False,
        "thread_safe": True,
        "thread_parallel": True,
        "node_parallel": False,
        "managed_memory": False,
    }
    version_cache: Dict[str, str] = {}

    class Config(ProgramHarness.Config):
        pass

    @staticmethod
    def found(raise_error: bool = False) -> bool:
        """Check for the availability of the Python API of dftd3"""

        return which_import(
            "dftd3",
            return_bool=True,
            raise_error=raise_error,
            raise_msg="Please install a dftd3 version with enabled Python API"
            + " (e.g. conda install dftd3-python -c conda-forge)",
        )

    def get_version(self) -> str:
        """Return the currently used version of dftd3"""
        self.found(raise_error=True)

        which_prog = which_import("dftd3")
        if which_prog not in self.version_cache:
            import dftd3

            self.version_cache[which_prog] = safe_version(dftd3.__version__)

        return self.version_cache[which_prog]

    def compute(self, input_model: AtomicInput, config: TaskConfig) -> AtomicResult:
        """
        Actual interface to the dftd3 package. The compute function is just a thin
        wrapper around the native QCSchema interface of the dftd3 Python-API.
        """
        self.found(raise_error=True)
        if parse_version(self.get_version()) < parse_version("0.5.1"):
            raise ResourceError("QCEngine's dftd3 wrapper requires version 0.5.1 or greater.")

        import dftd3
        from dftd3.qcschema import run_qcschema

        # strip engine hint
        input_data = input_model.dict()
        method = input_model.model.method
        if method.startswith("d3-"):
            method = method[3:]
            input_data["model"]["method"] = method
        qcvkey = method.upper() if method is not None else None

        # send `from_arrays` the s-dftd3 behavior of functional specification overrides explicit parameters specification
        # * differs from dftd3 harness behavior where parameters extend or override functional
        # * stash the resolved plan in extras or, if errored, leave it for the proper dftd3 api to reject
        param_tweaks = None if method else input_model.keywords.get("params_tweaks", None)
        try:
            planinfo = from_arrays(
                verbose=1,
                name_hint=method,
                level_hint=input_model.keywords.get("level_hint", None),
                param_tweaks=param_tweaks,
                dashcoeff_supplement=input_model.keywords.get("dashcoeff_supplement", None),
            )
        except InputError:
            pass
        else:
            input_data["extras"]["info"] = planinfo

        # strip dispersion level from method
        for alias, d3 in get_dispersion_aliases().items():
            if d3.startswith("d3") and method.lower().endswith(alias):
                method = method[: -(len(alias) + 1)]
                input_data["model"]["method"] = method

        # consolidate dispersion level aliases
        if input_model.keywords.pop("apply_qcengine_aliases", False):
            level_hint = input_model.keywords.get("level_hint", None)
            if level_hint:
                level_hint = get_dispersion_aliases()[level_hint.lower()]
                if level_hint.endswith("atm"):
                    level_hint = level_hint[:-3]
                if level_hint.endswith("2b"):
                    level_hint = level_hint[:-2]
                    input_data["keywords"]["params_tweaks"] = {**planinfo["dashparams"], "s9": 0.0}
                input_data["keywords"]["level_hint"] = level_hint

        input_model = AtomicInput(**input_data)

        # Run the Harness
        output = run_qcschema(input_model)

        if "info" in output.extras:
            qcvkey = output.extras["info"]["fctldash"].upper()

        calcinfo = {}
        energy = output.properties.return_energy
        calcinfo["CURRENT ENERGY"] = energy
        calcinfo["DISPERSION CORRECTION ENERGY"] = energy
        if qcvkey:
            calcinfo[f"{qcvkey} DISPERSION CORRECTION ENERGY"] = energy

        if output.driver == "gradient":
            gradient = output.return_result
            calcinfo["CURRENT GRADIENT"] = gradient
            calcinfo["DISPERSION CORRECTION GRADIENT"] = gradient
            if qcvkey:
                calcinfo[f"{qcvkey} DISPERSION CORRECTION GRADIENT"] = gradient

        if output.keywords.get("pair_resolved", False):
            pw2 = output.extras["dftd3"]["additive pairwise energy"]
            pw3 = output.extras["dftd3"]["non-additive pairwise energy"]
            assert abs(pw2.sum() + pw3.sum() - energy) < 1.0e-8, f"{pw2.sum()} + {pw3.sum()} != {energy}"
            calcinfo["2-BODY DISPERSION CORRECTION ENERGY"] = pw2.sum()
            calcinfo["3-BODY DISPERSION CORRECTION ENERGY"] = pw3.sum()
            calcinfo["2-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS"] = pw2
            calcinfo["3-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS"] = pw3

        output.extras["qcvars"] = calcinfo

        return output
