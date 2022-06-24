"""Harness for TeraChem Frontend"""
import logging
from os import getenv
from typing import Any, Dict

from .terachem_pbs import TeraChemPBSHarness, _pbs_defaults

logger = logging.getLogger(__name__)

_fe_defaults = {
    "name": "terachem_fe",
}


class TeraChemFrontEndHarness(TeraChemPBSHarness):
    """QCEngine Harness for interfacing with the TeraChem Frontend (Protocol Buffer Server + file server)"""

    _defaults = {**_pbs_defaults, **_fe_defaults}
    _tcpb_min_version = "0.9.0"
    _tcpb_client = "TCFrontEndClient"
    _env_vars: Dict[str, Any] = {
        **TeraChemPBSHarness._env_vars,
        **{
            "frontend_host": getenv("TERACHEM_FE_HOST", "127.0.0.1"),
            "frontend_port": int(getenv("TERACHEM_FE_PORT", 80)),
            "uploads_prefix": getenv("TERACHEM_FE_UPLOADS_PREFIX", "uploads"),
        },
    }
    _env_vars_external: str = (
        TeraChemPBSHarness._env_vars_external + " and/or TERACHEM_FE_HOST, TERACHEM_FE_PORT, TERACHEM_FE_UPLOADS_PREFIX"
    )
