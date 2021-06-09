"""Known errors for NWChem"""
from typing import Any, Dict

from qcelemental.models import AtomicInput

from qcengine.exceptions import SimpleKnownErrorException


class GeomBinvrError(SimpleKnownErrorException):
    """GeomBinvr Error"""

    error_name = "geom_binvr"
    description = "Error when generating redundant atomic coordinates. Fixed by turning this feature off (`noautoz`)"

    @classmethod
    def _detect(cls, outputs: Dict[str, str]) -> bool:
        return "geom_binvr: #indep variables incorrect" in outputs["stdout"]

    def create_keyword_update(self, input_data: AtomicInput) -> Dict[str, Any]:
        return {"geometry__noautoz": True}


# List of all the known errors
all_errors = [GeomBinvrError]
