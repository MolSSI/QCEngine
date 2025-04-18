import numpy as np
import pytest
import qcelemental as qcel

import qcengine as qcng
from qcengine.testing import using


@using("dftd4")
def test_dftd4():
    #! Exercises the various DFT-D corrections, both through python directly and through c++
    import psi4
    import numpy as np

    ref_d2         = [-0.00390110, -0.00165271, -0.00058118]
    ref_d3zero     = [-0.00285088, -0.00084340, -0.00031923]
    ref_d3bj       = [-0.00784595, -0.00394347, -0.00226683]
    ref_d4bj       = [-0.00625445, -0.00306407, -0.00176150]
    ref_d4bj_2body = [-0.00625366, -0.00306413, -0.00176146]

    ref_pbe_d2     = [-0.00278650, -0.00118051, -0.00041513]
    ref_pbe_d3zero = [-0.00175474, -0.00045421, -0.00016839]
    ref_pbe_d3bj   = [-0.00475937, -0.00235265, -0.00131239]
    ref_pbe_d4bj   = [-0.00399176, -0.00190682, -0.00108739]

    gref_d4bj = [
                np.array(
                    [
                        [-0.0, -0.000121388195, -0.000285720303],
                        [0.0, 0.000121388195, -0.000285720303],
                        [0.000013500589, -0.000035224119, -0.000051913946],
                        [-0.000013500589, -0.000035224119, -0.000051913946],
                        [-0.000013500589, 0.000035224119, -0.000051913946],
                        [0.000013500589, 0.000035224119, -0.000051913946],
                        [0.0, -0.0, 0.000162779428],
                        [0.0, -0.0, 0.00042515118],
                        [0.0, -0.0, 0.000180544011],
                        [0.0, -0.0, 0.000010621775],
                    ]
                ),
                np.array(
                    [
                        [-0.0, -0.000077703214, -0.000000117582],
                        [0.0, 0.000077703214, -0.000000117582],
                        [-0.000003720109, -0.000014373039, 0.000000058791],
                        [0.000003720109, -0.000014373039, 0.000000058791],
                        [0.000003720109, 0.000014373039, 0.000000058791],
                        [-0.000003720109, 0.000014373039, 0.000000058791],
                    ]
                ),
                np.array(
                    [
                        [0.0, 0.0, 0.000044445504],
                        [0.0, 0.0, -0.000044316404],
                        [0.0, 0.0, -0.000011452466],
                        [0.0, 0.0, 0.000011323366],
                    ]
                ),
    ]


    eneyne = psi4.geometry("""
    C   0.000000  -0.667578  -2.124659
    C   0.000000   0.667578  -2.124659
    H   0.923621  -1.232253  -2.126185
    H  -0.923621  -1.232253  -2.126185
    H  -0.923621   1.232253  -2.126185
    H   0.923621   1.232253  -2.126185
    --
    C   0.000000   0.000000   2.900503
    C   0.000000   0.000000   1.693240
    H   0.000000   0.000000   0.627352
    H   0.000000   0.000000   3.963929
    """)

    print('  -D correction from Py-side')
    eneyne.update_geometry()
    E, G = eneyne.run_dftd4('b3lyp', 'd4bj')
    assert psi4.compare_values(ref_d4bj[0], E, 7, 'Ethene-Ethyne -D4 (bj)')
    mA = eneyne.extract_subsets(1)
    E, G = mA.run_dftd4('b3lyp', 'd4bj')
    assert psi4.compare_values(ref_d4bj[1], E, 7, 'Ethene -D4 (bj)')
    assert psi4.compare_values(gref_d4bj[1], G, 7, 'Ethene -D4 (bj) grad')
    mB = eneyne.extract_subsets(2)
    E, G = mB.run_dftd4('b3lyp', 'd4bj')
    assert psi4.compare_values(ref_d4bj[2], E, 7, 'Ethyne -D4 (bj)')

    psi4.set_options({
        "basis": "sto-3g",
        "scf_type": "df",
        "dft_radial_points": 50,  # use really bad grid for speed since all we want is the -D value
        "dft_spherical_points": 110,
    })

    print('  -D correction from C-side')
    psi4.activate(mA)
    psi4.gradient('b3lyp-d4bj')
    assert psi4.compare_values(ref_d4bj[1], psi4.variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D4')
    assert psi4.compare_values(gref_d4bj[1], psi4.variable('DISPERSION CORRECTION GRADIENT'), 7, 'Ethene -D4 grad')

    mAgB = eneyne.extract_subsets(1, 2)
    psi4.activate(mAgB)
    psi4.energy('b3lyp-d4')
    assert psi4.compare_values(ref_d4bj[1], psi4.variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D4 (alias)')

    print('  non-default -D correction from C-side')
    psi4.activate(mB)
    psi4.set_options({"dft_dispersion_parameters": [0.40868035, 4.53807137, 16.0, 1.0, 2.02929367, 1.0]})  # b3lyp-d4
    psi4.energy('b3lyp-d4(bj)')
    assert psi4.compare_values(ref_d4bj[2], psi4.variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D4')

    psi4.activate(eneyne)
    psi4.set_options({"dft_dispersion_parameters": [0.40868035, 4.53807137, 16.0, 1.0, 2.02929367, 0.0]})  # b3lyp-d4-2body
    psi4.energy('b3lyp-d4(bj)')
    assert psi4.compare_values(ref_d4bj_2body[0], psi4.variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D4-2body')

    gAmB = eneyne.extract_subsets(2, 1)
    psi4.activate(gAmB)
    psi4.set_options({"dft_dispersion_parameters": [0.75]})
    psi4.energy('b3lyp-d2')
    assert psi4.compare_values(ref_pbe_d2[2], psi4.variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D2 (alias)')

    psi4.set_options({"dft_dispersion_parameters": [1.0,  0.722, 1.217, 14.0]})
    psi4.energy('b3lyp-d3')
    assert psi4.compare_values(ref_pbe_d3zero[2], psi4.variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D4 (alias)')

    psi4.set_options({"dft_dispersion_parameters": [0.38574991, 4.80688534, 16.0, 1.0, 0.95948085, 1.0]})  # pbe-d4
    psi4.energy('b3lyp-d4')
    assert psi4.compare_values(ref_pbe_d4bj[2], psi4.variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene -D4 (alias)')
    psi4.activate(mA)
    psi4.set_options({"dft_dispersion_parameters": [1.0]})
    psi4.energy('wb97x-d')
    assert psi4.compare_values(-0.000834247063, psi4.variable('DISPERSION CORRECTION ENERGY'), 7, 'Ethene wb97x-d (chg)')

