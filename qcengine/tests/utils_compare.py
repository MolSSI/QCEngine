import sys
import copy
import math
import pprint
import collections

import numpy as np


class _TestComparisonError(Exception):
    """Error when element or nuclide can't be identified."""

    def __init__(self, msg):
        self.message = '\nTestComparisonError: {}\n\n'.format(msg)


def _success(label):
    """Function to print a '*label*...PASSED' line to screen."""

    print('\t{0:.<66}PASSED'.format(label))
    sys.stdout.flush()


def compare_values(expected, computed, digits, label, passnone=False, verbose=1):
    """Function to compare two values to `digits` precision. Prints :py:func:`_success`.

    Raises
    ------
    _TestComparisonError
        If `computed` differs from `expected` by more than number of `digits`.
        (or by `digits` itself when <= 1 e.g. digits=0.04).

    """
    if passnone:
        if expected is None and computed is None:
            _success(label)
            return

    if digits > 1:
        thresh = 10**-digits
        message = """\t{}: computed value ({:.{digits1}f}) does not match ({:.{digits1}f}) to {digits} digits.""".format(
            label, computed, expected, digits1=int(digits) + 1, digits=digits)
    else:
        thresh = digits
        message = ("\t%s: computed value (%f) does not match (%f) to %f digits." % (label, computed, expected, digits))
    if abs(float(expected) - float(computed)) > thresh:
        raise _TestComparisonError(message)
    if math.isnan(computed):
        message += "\tprobably because the computed value is nan."
        raise _TestComparisonError(message)
    if verbose >= 1:
        _success(label)


def compare_integers(expected, computed, label, verbose=1):
    """Function to compare two integers or strings.

    Raises
    ------
    _TestComparisonError
        If `computed` doesn't exactly match `expected`.

    """
    if expected != computed:
        message = "\t{}: computed value ({}) does not match ({}).".format(label, computed, expected)
        raise _TestComparisonError(message)
    if verbose >= 1:
        _success(label)


compare_strings = compare_integers


def compare_dicts(expected, computed, tol, label, forgive=None, verbose=1):
    """Compares dictionaries `computed` to `expected` using DeepDiff Float
    comparisons made to `tol` significant decimal places. Note that a clean
    DeepDiff returns {}, which evaluates to False, hence the compare_integers.
    Keys in `forgive` may change between `expected` and `computed` without
    triggering failure.

    """
    try:
        import deepdiff
    except ImportError:
        raise ImportError("""Install deepdiff. `conda install deepdiff -c conda-forge` or `pip install deepdiff`""")

    if forgive is None:
        forgive = []
    forgiven = collections.defaultdict(dict)

    ans = deepdiff.DeepDiff(expected, computed, significant_digits=tol, verbose_level=2)

    for category in list(ans):
        for key in list(ans[category]):
            for fg in forgive:
                fgsig = "root['" + fg + "']"
                if key.startswith(fgsig):
                    forgiven[category][key] = ans[category].pop(key)
        if not ans[category]:
            del ans[category]

    clean = not bool(ans)
    if not clean:
        pprint.pprint(ans)
    if verbose >= 2:
        pprint.pprint(forgiven)
    return compare_integers(True, clean, label, verbose=verbose)


def compare_molrecs(expected, computed, tol, label, forgive=None, verbose=1, relative_geoms='exact'):
    """Function to compare Molecule dictionaries. Prints
#    :py:func:`util.success` when elements of `computed` match elements of
#    `expected` to `tol` number of digits (for float arrays).

    """
    thresh = 10**-tol if tol >= 1 else tol

    # Need to manipulate the dictionaries a bit, so hold values
    xptd = copy.deepcopy(expected)
    cptd = copy.deepcopy(computed)

    def massage_dicts(dicary):
        # deepdiff can't cope with np.int type
        #   https://github.com/seperman/deepdiff/issues/97
        if 'elez' in dicary:
            dicary['elez'] = [int(z) for z in dicary['elez']]
        if 'elea' in dicary:
            dicary['elea'] = [int(a) for a in dicary['elea']]
        # deepdiff w/py27 complains about unicode type and val errors
        if 'elem' in dicary:
            dicary['elem'] = [str(e) for e in dicary['elem']]
        if 'elbl' in dicary:
            dicary['elbl'] = [str(l) for l in dicary['elbl']]
        if 'fix_symmetry' in dicary:
            dicary['fix_symmetry'] = str(dicary['fix_symmetry'])
        if 'units' in dicary:
            dicary['units'] = str(dicary['units'])
        if 'fragment_files' in dicary:
            dicary['fragment_files'] = [str(f) for f in dicary['fragment_files']]
        # and about int vs long errors
        if 'molecular_multiplicity' in dicary:
            dicary['molecular_multiplicity'] = int(dicary['molecular_multiplicity'])
        if 'fragment_multiplicities' in dicary:
            dicary['fragment_multiplicities'] = [(m if m is None else int(m))
                                                 for m in dicary['fragment_multiplicities']]
        if 'fragment_separators' in dicary:
            dicary['fragment_separators'] = [(s if s is None else int(s)) for s in dicary['fragment_separators']]
        # forgive generator version changes
        if 'provenance' in dicary:
            dicary['provenance'].pop('version')
        # regularize connectivity ordering
        if 'connectivity' in dicary:
            conn = [(min(at1, at2), max(at1, at2), bo) for (at1, at2, bo) in dicary['connectivity']]
            conn.sort(key=lambda tup: tup[0])
            dicary['connectivity'] = conn

        return dicary

    xptd = massage_dicts(xptd)
    cptd = massage_dicts(cptd)

    if relative_geoms == 'exact':
        pass
    elif relative_geoms == 'align':
        raise FeatureNotImplemented(
            """compare_molrecs(..., relative_geoms='align') not available without B787 from qcdb.""")
        ## can't just expect geometries to match, so we'll align them, check that
        ##   they overlap and that the translation/rotation arrays jibe with
        ##   fix_com/orientation, then attach the oriented geom to computed before the
        ##   recursive dict comparison.
        #from .align import B787
        #cgeom = np.array(cptd['geom']).reshape((-1, 3))
        #rgeom = np.array(xptd['geom']).reshape((-1, 3))
        #rmsd, mill = B787(rgeom=rgeom,
        #                  cgeom=cgeom,
        #                  runiq=None,
        #                  cuniq=None,
        #                  atoms_map=True,
        #                  mols_align=True,
        #                  run_mirror=False,
        #                  verbose=0)
        #if cptd['fix_com']:
        #    compare_integers(1, np.allclose(np.zeros((3)), mill.shift, atol=thresh), 'null shift', verbose=verbose)
        #if cptd['fix_orientation']:
        #    compare_integers(1, np.allclose(np.identity(3), mill.rotation, atol=thresh), 'null rotation', verbose=verbose)
        #ageom = mill.align_coordinates(cgeom)
        #cptd['geom'] = ageom.reshape((-1))

    compare_dicts(xptd, cptd, tol, label, forgive=forgive, verbose=verbose)


def compare_arrays(expected, computed, digits, label, rtol=1.e-16, verbose=1):
    """Function to compare two numpy arrays to `digits` precision. Prints :py:func:`_success`.

    Sets rtol to zero to match expected Psi4 behaviour, otherwise measured as:
        absolute(computed - expected) <= (atol + rtol * absolute(expected))

    Raises
    ------
    _TestComparisonError
        If `computed` differs from `expected` by more than number of `digits`.
        (or by `digits` itself when <= 1 e.g. digits=0.04).

    """
    try:
        expected = np.asarray(expected)
        computed = np.asarray(computed)
        shape1 = expected.shape
        shape2 = computed.shape
    except:
        raise TypeError("Input objects do not have a shape attribute.")

    if shape1 != shape2:
        raise TypeError("Input shapes do not match.")

    tol = 10**(-digits)
    if not np.allclose(expected, computed, atol=tol, rtol=rtol):
        message = "\tArray difference norm is %12.6e." % np.linalg.norm(expected - computed)
        raise _TestComparisonError(message)

    if verbose >= 1:
        _success(label)
