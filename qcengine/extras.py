
"""
Misc information and runtime information.
"""

from . import _version

__all__ = ["get_information"]

versions = _version.get_versions()

__info = {"version": versions['version'], "git_revision": versions['full-revisionid']}


def get_information(key):
    """
    Obtains a variety of runtime information about QCEngine.
    """
    key = key.lower()
    if key not in __info:
        raise KeyError("Information key '{}' not understood.".format(key))

    return __info[key]
