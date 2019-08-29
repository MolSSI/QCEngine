import copy
import pytest

from qcelemental.testing import compare, compare_recursive, compare_values

import qcengine as qcng
from qcengine.testing import using_pylibefp

b2a = 0.529177
a2b = 1.0 / b2a

_blank_qm_mol = {
    'geometry': [0, 0, 0],  # dummy
    'symbols': ['He'],  # dummy
    'extras': {
        'efp_molecule': {
            'geometry': [],
            'symbols': [],
            'extras': None
        }
    }
}


def _blank_ene(label='libefp'):
    fields = [
        'charge_penetration', 'disp', 'dispersion', 'electrostatic', 'electrostatic_point_charges',
        'exchange_repulsion', 'polarization'
    ]
    if label == 'libefp':
        fields.extend(['elec', 'pol', 'xr'])
    elif label == 'psi':
        fields.extend(['elst', 'ind', 'exch'])

    return {f: 0.0 for f in fields}


@pytest.fixture
def system_1():
    kmol = copy.deepcopy(_blank_qm_mol)
    kmol['extras']['efp_molecule']['extras'] = {
        'fragment_files': ['h2o', 'nh3'],
        'hint_types': ['xyzabc', 'xyzabc'],
        'geom_hints': [
            [0.0 * a2b, 0.0 * a2b, 0.0 * a2b, 1.0, 2.0, 3.0],
            [5.0 * a2b, 0.0 * a2b, 0.0 * a2b, 5.0, 2.0, 8.0],
        ]
    }
    return kmol


@pytest.fixture
def system_2():
    kmol = copy.deepcopy(_blank_qm_mol)
    kmol['extras']['efp_molecule']['extras'] = {
        'fragment_files': ['h2o', 'nh3', 'h2o', 'h2o', 'nh3'],
        'hint_types': ['xyzabc'] * 5,
        'geom_hints': [
            [-1.0 * a2b,  3.7 * a2b,  0.4 * a2b, -1.3,  0.0,  7.0],
            [ 0.4 * a2b, -0.9 * a2b, -0.7 * a2b,  4.0,  1.6, -2.3],
            [ 1.7 * a2b,  2.0 * a2b,  3.3 * a2b, -1.2, -2.0,  6.2],
            [ 0.0 * a2b,  3.9 * a2b, -3.4 * a2b,  1.3,  5.2, -3.0],
            [-3.5 * a2b,  0.0 * a2b, -0.7 * a2b,  0.0, -2.7,  2.7],
        ]
    }  # yapf: disable

    return kmol


@pytest.fixture
def system_3():
    kmol = copy.deepcopy(_blank_qm_mol)
    kmol['extras']['efp_molecule']['extras'] = {
        'fragment_files': ['h2o', 'nh3', 'nh3', 'nh3', 'ch3oh', 'h2o', 'h2o', 'ch3oh', 'h2o'],
        'hint_types': ['xyzabc'] * 9,
        'geom_hints': [
            [i * a2b for i in [-3.394, -1.900, -3.700, -3.524, -1.089, -3.147, -2.544, -2.340, -3.445]],
            [i * a2b for i in [-5.515,  1.083,  0.968, -5.161,  0.130,  0.813, -4.833,  1.766,  0.609]],
            [i * a2b for i in [ 1.848,  0.114,  0.130,  1.966,  0.674, -0.726,  0.909,  0.273,  0.517]],
            [i * a2b for i in [-1.111, -0.084, -4.017, -1.941,  0.488, -3.813, -0.292,  0.525, -4.138]],
            [i * a2b for i in [-2.056,  0.767, -0.301, -2.999, -0.274, -0.551, -1.201,  0.360,  0.258]],
            [i * a2b for i in [-0.126, -2.228, -0.815,  0.310, -2.476,  0.037,  0.053, -1.277, -1.011]],
            [i * a2b for i in [-1.850,  1.697,  3.172, -1.050,  1.592,  2.599, -2.666,  1.643,  2.614]],
            [i * a2b for i in [ 1.275, -2.447, -4.673,  0.709, -3.191, -3.592,  2.213, -1.978, -4.343]],
            [i * a2b for i in [-5.773, -1.738, -0.926, -5.017, -1.960, -1.522, -5.469, -1.766,  0.014]],
        ]
    }  # yapf: disable

    return kmol


@pytest.fixture
def system_4():
    kmol = copy.deepcopy(_blank_qm_mol)
    kmol['extras']['efp_molecule']['extras'] = {
        'fragment_files':
        ['acetone', 'c2h5oh', 'c6h6', 'ccl4', 'ch3oh', 'ch4', 'cl2', 'dcm', 'dmso', 'h2', 'h2o', 'nh3'],
        'hint_types': ['xyzabc'] * 12,
        'geom_hints': [
            [ 0.0 * a2b,  0.0 * a2b, 0.0 * a2b, 0.0, 0.2, 0.3],
            [ 7.0 * a2b,  0.0 * a2b, 0.0 * a2b, 0.0, 2.0, 3.7],
            [14.0 * a2b,  0.0 * a2b, 0.0 * a2b, 3.1, 0.8, 2.0],
            [21.0 * a2b,  0.0 * a2b, 0.0 * a2b, 0.0, 8.0, 0.0],
            [ 0.0 * a2b,  6.0 * a2b, 0.0 * a2b, 0.7, 2.0, 1.0],
            [ 7.0 * a2b,  6.0 * a2b, 0.0 * a2b, 0.6, 0.0, 4.7],
            [14.0 * a2b,  6.0 * a2b, 0.0 * a2b, 0.0, 0.0, 0.3],
            [21.0 * a2b,  6.0 * a2b, 0.0 * a2b, 0.0, 0.4, 0.3],
            [ 0.0 * a2b, 12.0 * a2b, 0.0 * a2b, 0.8, 0.0, 0.0],
            [ 7.0 * a2b, 12.0 * a2b, 0.0 * a2b, 8.0, 0.7, 0.8],
            [14.0 * a2b, 12.0 * a2b, 0.0 * a2b, 0.0, 0.0, 0.0],
            [21.0 * a2b, 12.0 * a2b, 0.0 * a2b, 0.0, 2.0, 0.0],
        ]
    }  # yapf: disable

    return kmol


@pytest.fixture
def system_5():
    kmol = copy.deepcopy(_blank_qm_mol)
    kmol['extras']['efp_molecule']['extras'] = {
        'fragment_files': ['h2o', 'nh3'],
        'hint_types': ['xyzabc'] * 2,
        'geom_hints': [
            [ 0.0 * a2b,  0.0 * a2b,  0.0 * a2b, 3.0, 0.0, 7.0],
            [18.0 * a2b, 18.0 * a2b, 18.0 * a2b, 5.0, 4.0, 6.0],
        ]
    }  # yapf: disable

    return kmol


@pytest.fixture
def system_6():
    kmol = copy.deepcopy(_blank_qm_mol)
    kmol['extras']['efp_molecule']['extras'] = {
        'fragment_files': ['h2o', 'ch3oh', 'h2o', 'ch3oh', 'nh3'],
        'hint_types': ['xyzabc'] * 5,
        'geom_hints': [
            [ 0.0 * a2b,  0.0 * a2b,  0.0 * a2b, 0.0, 0.0, 0.0],
            [19.0 * a2b,  0.0 * a2b,  0.0 * a2b, 0.0, 0.0, 0.0],
            [ 0.0 * a2b, 19.0 * a2b,  0.0 * a2b, 0.0, 0.0, 0.0],
            [ 0.0 * a2b,  0.0 * a2b, 19.0 * a2b, 0.0, 0.0, 0.0],
            [18.0 * a2b, 18.0 * a2b, 18.0 * a2b, 0.0, 0.0, 0.0],
        ]
    }  # yapf: disable

    sys.prepare()
    return sys


@pytest.fixture
def system_qm1():
    kmol = copy.deepcopy(_blank_qm_mol)
    kmol['extras']['efp_molecule']['extras'] = {
        'fragment_files': ['h2o', 'c6h6', 'nh3'],
        'hint_types': ['xyzabc'] * 3,
        'geom_hints': [
            [-1.6 * a2b,  4.7 * a2b,  1.4 * a2b, -1.3,  0.1,  7.0],
            [ 0.4 * a2b, -0.9 * a2b, -0.7 * a2b,  2.3,  1.6, -2.3],
            [-3.5 * a2b, -2.0 * a2b, -0.7 * a2b,  0.0,  2.2,  2.7],
        ]
    }  # yapf: disable

    return kmol


@pytest.fixture
def system_qm2():
    kmol = copy.deepcopy(_blank_qm_mol)
    kmol['extras']['efp_molecule']['extras'] = {
        'fragment_files': ['ch3oh', 'dmso', 'dmso', 'acetone', 'dcm', 'acetone', 'acetone'],
        'hint_types': ['xyzabc'] * 7,
        'geom_hints': [
            [ 0.0 * a2b, -1.0 * a2b,  0.0 * a2b, 0.0, 1.1, 2.0],
            [-5.0 * a2b, 12.0 * a2b,  0.0 * a2b, 3.0, 0.2, 5.0],
            [ 0.0 * a2b, -3.0 * a2b,  5.0 * a2b, 6.0, 2.3, 8.0],
            [-5.0 * a2b, -4.0 * a2b, -5.0 * a2b, 9.0, 0.4, 1.0],
            [-9.0 * a2b, -5.0 * a2b,  1.0 * a2b, 2.0, 1.5, 4.0],
            [ 7.0 * a2b, -2.0 * a2b, 11.0 * a2b, 5.0, 0.6, 7.0],
            [-9.0 * a2b, -7.0 * a2b, -9.0 * a2b, 8.0, 2.7, 0.0],
        ]
    }  # yapf: disable

    return kmol


@using_pylibefp
@pytest.mark.parametrize("keywords_label,results_label", [
    ("libefp", "libefp"),
    ("libefp", "psi"),
    ("psi", "libefp"),
    ("psi", "psi"),
])
def test_total_1a(system_1, keywords_label, results_label):

    resi = {
        'molecule': system_1,
        'driver': 'energy',
        'model': {
            'method': 'efpefp',
        },
        'keywords': {
            'keywords_label': keywords_label,
            'results_label': results_label,
        }
    }

    # keywords_label changes both keywords names and defaults they modify
    if keywords_label == 'libefp':
        resi['keywords'].update({
            'elec': True,
            'elec_damp': 'screen',
            'xr': True,
            'pol': True,
            'disp': True,
            'disp_damp': 'tt',
        })
    elif keywords_label == 'psi':
        resi['keywords'].update({
            'disp_damping': 'tt',
        })

    res = qcng.compute(resi, 'pylibefp', raise_error=True, return_dict=False)

    atol = 1.e-6
    assert compare_values(0.0001922903, res.return_result, atol=atol)

    elst = 0.0002900482
    exch = 0.0000134716
    ind = 0.0002777238
    disp = -0.0000989033
    tot = 0.0001922903

    if results_label == 'libefp':
        expected_ene = _blank_ene()
        expected_ene['elec'] = expected_ene['electrostatic'] = elst
        expected_ene['xr'] = expected_ene['exchange_repulsion'] = exch
        expected_ene['pol'] = expected_ene['polarization'] = ind - expected_ene['electrostatic']
    elif results_label == 'psi':
        expected_ene = _blank_ene(label='psi')
        expected_ene['elst'] = expected_ene['electrostatic'] = elst
        expected_ene['exch'] = expected_ene['exchange_repulsion'] = exch
        expected_ene['ind'] = expected_ene['polarization'] = ind - expected_ene['electrostatic']
    expected_ene['disp'] = expected_ene['dispersion'] = disp
    expected_ene['total'] = tot

    # assert compare(2, asdf.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag')
    # assert compare_values(0.0, asdf.get_frag_charge(1), sys._getframe().f_code.co_name + ': f_chg', atol=1.e-6)
    # assert compare(1, asdf.get_frag_multiplicity(1), sys._getframe().f_code.co_name + ': f_mult')
    # assert compare('NH3', asdf.get_frag_name(1), sys._getframe().f_code.co_name + ': f_name')
    assert compare_recursive(expected_ene, res.extras['local_properties'], atol=atol)


@pytest.mark.parametrize("keywords,errmsg", [
    ({'nonsense_key': 'harmless', 'elec_damp': 'nonsense'}, "invalid value for [screen/overlap/off] elec_damp: nonsense"),
    ({'elec': 'yEs'}, "invalid value for [T/F] elec: yEs"),
])  # yapf: disable
def test_opts_fail_1(keywords, errmsg, system_1):
    resi = {
        'molecule': system_1,
        'driver': 'energy',
        'model': {
            'method': 'efpefp',
        },
        'keywords': keywords,
    }

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcng.compute(resi, 'pylibefp', raise_error=True, return_dict=False)

    assert errmsg in str(e.value)


"""
start by generating various halogenated benzenes
opt by b3lyp
run makefp to generate those files
run them in dimers in libefp
run them in dimers in sapt for reference
upload to qca
"""
