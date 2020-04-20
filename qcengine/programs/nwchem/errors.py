"""Known errors for NWChem"""
from typing import Dict

from qcengine.exceptions import SimpleKnownErrorException


class SymGeomProjectError(SimpleKnownErrorException):
    """sym_center_map is inconsistent with requested accuracy"""
    error_name = 'sym_geom_project'
    description = 'sym_center_map is inconsistent with requested accuracy. Occurs when geometry is loaded' \
                  ' and is resolved by lowering the tolerance when detecting symmetry.'

    @classmethod
    def _detect(cls, outputs):
        return "sym_geom_project: sym_center_map is inconsistent with requested accuracy" in outputs["stdout"]


class SymMapError(SimpleKnownErrorException):
    """sym_map: no match"""
    error_name = 'sym_map'
    description = 'Occurs when geometry is loaded and is resolved by lowering the tolerance when detecting symmetry.'

    @classmethod
    def _detect(cls, outputs):
        return "sym_map: no match" in outputs["stdout"]


class GeomBinvrError(SimpleKnownErrorException):
    """GeomBinvr Error"""

    error_name = 'geom_binvr'
    description = 'Error when generating redundant atomic coordinates. Fixed by turning this feature off (`noautoz`)'

    @classmethod
    def _detect(cls, outputs: Dict[str, str]) -> bool:
        return "geom_binvr: #indep variables incorrect" in outputs["stdout"]


# List of all the known errors
all_errors = [SymGeomProjectError, SymMapError, GeomBinvrError]
