import re

from ...extras import get_information

yes = re.compile(r'^(yes|true|on|1)', re.IGNORECASE)
no = re.compile(r'^(no|false|off|0)', re.IGNORECASE)
der0th = re.compile(r'^(0|none|energy)', re.IGNORECASE)
der1st = re.compile(r'^(1|first|gradient)', re.IGNORECASE)
der2nd = re.compile(r'^(2|second|hessian)', re.IGNORECASE)
der3rd = re.compile(r'^(3|third)', re.IGNORECASE)
der4th = re.compile(r'^(4|fourth)', re.IGNORECASE)
der5th = re.compile(r'^(5|fifth)', re.IGNORECASE)


def parse_dertype(dertype, max_derivative=2):
    """Apply generous regex to `dertype` to return regularized integer and driver values for derivative level.

    Parameters
    ----------
    dertype : int or str
        Interpretable as a derivative level, regardless of case or type.
    max_derivative : int, optional
        Derivative level above which should throw FeatureNotImplemented error.

    Returns
    -------
    (int, {'energy', 'gradient', 'hessian'})
        Returns dertype as an integer and a driver-valid string.

    """
    derdriver = dict(enumerate(['energy', 'gradient', 'hessian', 'third', 'fourth', 'fifth']))

    if der0th.match(str(dertype)):
        derint = 0
    elif der1st.match(str(dertype)):
        derint = 1
    elif der2nd.match(str(dertype)):
        derint = 2
    elif der3rd.match(str(dertype)):
        derint = 3
    elif der4th.match(str(dertype)):
        derint = 4
    elif der5th.match(str(dertype)):
        derint = 5
    else:
        raise ValidationError("""Requested derivative level ({}) not recognized.""".format(dertype))

    if derint > max_derivative:
        raise FeatureNotImplemented("""derivative level ({})""".format(derint))

    return (derint, derdriver[derint])


def provenance_stamp(routine):
    """Return dictionary satisfying QCSchema,
    https://github.com/MolSSI/QCSchema/blob/master/qcschema/dev/definitions.py#L23-L41
    with QCEngine's credentials for creator and version. The
    generating routine's name is passed in through `routine`.

    """
    return {'creator': 'QCEngine', 'version': get_information('version'), 'routine': routine}
