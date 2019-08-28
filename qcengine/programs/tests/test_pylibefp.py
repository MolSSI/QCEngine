from qcelemental.testing import compare, compare_recursive, compare_values

import qcengine as qcng
from qcengine.testing import using_pylibefp

b2a = 0.529177
a2b = 1.0 / b2a


@using_pylibefp
def test_total_1a():

    resi = {
        'molecule': {
            'geometry': [0, 0, 0],  # dummy
            'symbols': ['He'],  # dummy
            'extras': {
                'efp_molecule': {
                    'geometry': [],
                    'symbols': [],
                    'extras': {
                        'fragment_files': ['h2o', 'nh3'],
                        'hint_types': ['xyzabc', 'xyzabc'],
                        'geom_hints': [
                            [0.0 * a2b, 0.0 * a2b, 0.0 * a2b, 1.0, 2.0, 3.0],
                            [5.0 * a2b, 0.0 * a2b, 0.0 * a2b, 5.0, 2.0, 8.0],
                        ]
                    }
                }
            }
        },
        'driver': 'energy',
        'model': {
            'method': 'efpefp',
        },
        'keywords': {
            'elec': True,
            'elec_damp': 'screen',
            'xr': True,
            'pol': True,
            'disp': True,
            'disp_damp': 'tt',
        }
    }

    res = qcng.compute(resi, 'pylibefp', raise_error=True, return_dict=False)

    atol = 1.e-6
    assert compare_values(0.0001922903, res.return_result, atol=atol)


#    expected_ene = blank_ene()
#    expected_ene['elec'] = expected_ene['electrostatic'] = 0.0002900482
#    expected_ene['xr'] = expected_ene['exchange_repulsion'] = 0.0000134716
#    expected_ene['pol'] = expected_ene['polarization'] = 0.0002777238 - expected_ene['electrostatic']
#    expected_ene['disp'] = expected_ene['dispersion'] = -0.0000989033
#    expected_ene['total'] = 0.0001922903
#    assert compare(2, asdf.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag')
#    assert compare_values(0.0, asdf.get_frag_charge(1), sys._getframe().f_code.co_name + ': f_chg', atol=1.e-6)
#    assert compare(1, asdf.get_frag_multiplicity(1), sys._getframe().f_code.co_name + ': f_mult')
#    assert compare('NH3', asdf.get_frag_name(1), sys._getframe().f_code.co_name + ': f_name')
#    assert compare_recursive(expected_ene, ene, sys._getframe().f_code.co_name + ': ene', atol=1.e-6)
