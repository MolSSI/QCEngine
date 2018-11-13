"""
A set of scripts to setup testing
"""

import pytest


def _plugin_import(plug):
    """
    Tests to see if a module is available
    """
    import sys
    if sys.version_info >= (3, 4):
        from importlib import util
        plug_spec = util.find_spec(plug)
    else:
        import pkgutil
        plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True

# Add flags
using_psi4 = pytest.mark.skipif(
    _plugin_import("psi4") is False, reason="Could not find psi4. Please install the package to enable tests")
using_rdkit = pytest.mark.skipif(
    _plugin_import("rdkit") is False, reason="Could not find rdkit. Please install the package to enable tests")
using_geometric = pytest.mark.skipif(
    _plugin_import("geometric") is False, reason="could not find geomeTRIC. Please install the package to enable tests")
using_torchani = pytest.mark.skipif(
    _plugin_import("torchani") is False, reason="Could not find TorchAni. Please install the package to enable tests")
