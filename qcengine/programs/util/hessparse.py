import numpy as np

from qcelemental.util import filter_comments

def load_hessian(shess, dtype):

    # list o'lines w/o comments or blanks
    shess = filter_comments(shess)
    lhess = list(filter(None, map(str.strip, shess.splitlines())))

    if dtype in ['fcmfinal', 'cfour']:
        nat = int(lhess[0].split()[0])
        ndof = 3 * nat
        datastr = '\n'.join(lhess[1:])
        nhess = np.fromstring(datastr, sep=' ')
        nhess = nhess.reshape(ndof, ndof)
    else:
        raise ValidationError('Unknown dtype: {}'.format(dtype))

    return nhess


def hess_to_string(hess, handle, dtype):
    nat = hess.shape[0] // 3
    assert hess.shape == (3 * nat, 3 * nat)

    header = '{:5}{:5}'.format(nat, 6 * nat)
    np.savetxt(handle, hess.reshape((-1, 3)), fmt='%20.10f', delimiter='', newline='\n', header=header, comments='')


#    fcm = fcm.splitlines()
#    Nat = int(fcm[0].split()[0])
#    Ndof = int(fcm[0].split()[1])
#
#    empty = True
#    hess = []
#    for df in range(Ndof):
#        for at in range(Nat):
#            lline = fcm[Ndof * at + at + 1].split()
#            if empty:
#                if (abs(float(lline[0])) > 1.0e-8) or \
#                   (abs(float(lline[1])) > 1.0e-8) or \
#                   (abs(float(lline[2])) > 1.0e-8):
#                    empty = False
#            fcm.append([float(lline[0]), float(lline[1]), float(lline[2])])
#
#    return None if empty else hess
