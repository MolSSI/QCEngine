import pytest


def pytest_addoption(parser):
    parser.addoption("--runnonroutine", action="store_true", default=False, help="run the nonroutine tests in stdsuite")


def pytest_collection_modifyitems(config, items):
    # psi4 on ECP systems needs libecpint compiled in and psi4's ECP fixes. checked only if such tests are collected
    psi4_ecp_items = [item for item in items if item.get_closest_marker("ecp") and item.get_closest_marker("psi4")]
    if psi4_ecp_items:
        from qcengine.testing import has_program

        if not has_program("psi4_ecp"):
            skip_psi4_ecp = pytest.mark.skip(
                reason="psi4 lacks libecpint (psi4.addons('ecpint')) or is older than v1.12"
            )
            for item in psi4_ecp_items:
                item.add_marker(skip_psi4_ecp)

    if config.getoption("--runnonroutine"):
        # --runnonroutine given in cli: do not skip nonroutine tests
        return
    skip_nonroutine = pytest.mark.skip(reason="need --runnonroutine option to run")
    for item in items:
        if "nonroutine" in item.keywords:
            item.add_marker(skip_nonroutine)
