import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runnonroutine", action="store_true", default=False, help="run the nonroutine tests in stdsuite"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runnonroutine"):
        # --runnonroutine given in cli: do not skip nonroutine tests
        return
    skip_nonroutine = pytest.mark.skip(reason="need --runnonroutine option to run")
    for item in items:
        if "nonroutine" in item.keywords:
            item.add_marker(skip_nonroutine)
