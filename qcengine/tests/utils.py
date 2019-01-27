import sys

from . import utils_compare


def _true_false_decorator(compare_fn, *args, **kwargs):
    """Turns `compare_fn` that returns `None` on success and raises
    `_TestComparisonError` on failure into a function that returns
    True/False, suitable for assertions in pytest.

    """

    def true_false_wrapper(*args, **kwargs):
        try:
            compare_fn(*args, **kwargs)
        except utils_compare._TestComparisonError as err:
            return False
        else:
            return True

    return true_false_wrapper


compare_dicts = _true_false_decorator(utils_compare.compare_dicts)
compare_arrays = _true_false_decorator(utils_compare.compare_arrays)
compare_values = _true_false_decorator(utils_compare.compare_values)
compare_molrecs = _true_false_decorator(utils_compare.compare_molrecs)
compare_strings = _true_false_decorator(utils_compare.compare_strings)
compare_integers = _true_false_decorator(utils_compare.compare_integers)


def tnm():
    """Returns the name of the calling function, usually name of test case."""

    return sys._getframe().f_back.f_code.co_name
