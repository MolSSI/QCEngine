"""
Imports the various compute backends
"""
# from .psi import run_psi4
# from .rdkit import run_rdkit
# from .torchani import run_torchani

import os
import re

programs_path = os.path.dirname(os.path.abspath(__file__))

# Patterns to ignore in the programs directory when looking for programs
ignore_files = r'^\.|^__init__.py$|^#'

DESCRIPTION = 'description'

programs = []
for file in os.listdir(programs_path):
    if file.endswith('.py') and not re.search(ignore_files, file):
        program = re.sub(r'.py$', '', file)
        programs.append(program)
programs.sort()


def _get_python_name(name):
    return name.replace('-', '_')


def _attr_setdefault(obj, name, value):
    """Like dict.setdefault, but for objects."""
    if not hasattr(obj, name):
        setattr(obj, name, value)
    return getattr(obj, name)


def get_module(name):
    """Imports the module for a particular program and returns it."""
    module_name = '%s.%s' % (__name__, name)
    module = __import__(module_name,
                        fromlist=[name, DESCRIPTION],
                        level=0)

    _attr_setdefault(module, DESCRIPTION, "")

    fn_name = _get_python_name(name)
    if not hasattr(module, fn_name):
        raise RuntimeError(
            "Program module %s (%s) must define function '%s'." % (module.__name__, module.__file__, fn_name))

    return module


def get_program(name):
    python_name = _get_python_name(name)
    return getattr(get_module(python_name), python_name)
