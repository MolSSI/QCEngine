"""Known errors for NWChem"""
from typing import Any, Dict

from qcelemental.models import AtomicInput

from qcengine.exceptions import SimpleKnownErrorException
from qcengine.programs.nwchem.germinate import xc_functionals


class GeomBinvrError(SimpleKnownErrorException):
    """GeomBinvr Error"""

    error_name = "geom_binvr"
    description = "Error when generating redundant atomic coordinates. Fixed by turning this feature off (`noautoz`)"

    @classmethod
    def _detect(cls, outputs: Dict[str, str]) -> bool:
        return "geom_binvr: #indep variables incorrect" in outputs["stdout"]

    def create_keyword_update(self, input_data: AtomicInput) -> Dict[str, Any]:
        return {"geometry__noautoz": True}


class ConvergenceFailedError(SimpleKnownErrorException):
    """Failed to converge on electronic steps. This will increase the limit by 2x"""

    error_name = "convergence_failed"
    description = "The computation failed to converge. We increased the number of steps by a factor of 2"

    @classmethod
    def _detect(cls, outputs: Dict[str, str]) -> bool:
        return (
            "This type of error is most commonly" in outputs["stdout"] and "convergence criteria" in outputs["stdout"]
        )

    def create_keyword_update(self, input_data: AtomicInput) -> Dict[str, Any]:
        # Fit the correct keyword we are looking to update is different for different methods
        method = input_data.model.method
        if method == "dft" or method.split()[0] in xc_functionals:
            if "dft__iterations" in input_data.keywords:
                kwd = "dft__iterations"
                cur_iter = input_data.keywords["dft__iterations"]
            elif "dft__maxiter" in input_data.keywords:
                kwd = "dft__maxiter"
                cur_iter = input_data.keywords["dft__maxiter"]
            else:
                kwd = "dft__maxiter"
                cur_iter = 20  # The NWChem default
        else:
            raise ValueError(f'Method "{method}" is not yet supported')

        return {kwd: cur_iter * 2}


# List of all the known errors
all_errors = [GeomBinvrError, ConvergenceFailedError]
