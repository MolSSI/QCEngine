import numpy as np

# This file contains reference geometries, energies, and gradients
#  for comparing QC program results. It is placed in QCEngine but is also
#  used directly by QCDB and Psi4. New or filled-out systems+modelchems
#  welcome, but it's good to start a discussion on GitHub to check that
#  its downstream roles not impinged.

std_refs = ["rhf", "uhf", "rohf"]

std_molecules = {
    "hf": """
                H
                F 1 0.917
symmetry c1  # TODO for nwchem, manage another way
              """,
    "bh3p": """
                1 2
                B     0.10369114     0.00000000     0.00000000
                H    -1.13269886     0.00000000     0.00000000
                H     3.00000000     0.37149000     0.00000000
                H     3.00000000    -0.37149000     0.00000000
              """,
    "h2o": """
                O
                H 1 R
                H 1 R 2 A
            
                R=0.958
                A=104.5
              """,
    "nh2": """
                0 2
                N
                H 1 R
                H 1 R 2 A

                R=1.008
                A=105.0
              """,
    "h2o-xyz": """
 # R=0.958 A=104.5
 H                  0.000000000000     1.431430901356     0.984293362719
 O                  0.000000000000     0.000000000000    -0.124038860300
 H                  0.000000000000    -1.431430901356     0.984293362719
 units au
              """,
    "nh2-xyz": """
 # R=1.008 #A=105.0
 0 2
 N   0.000000000000000   0.000000000000000  -0.145912918634892
 H   0.000000000000000  -1.511214298139000   1.013682596946108
 H   0.000000000000000   1.511214298139000   1.013682596946108
 units au
              """,
    "hf-xyz": """
        H    0.          0.         -1.64558411
        F    0.          0.          0.08729475
        symmetry c1
        units bohr
               """,
}
std_molecules["bh3p-xyz"] = std_molecules["bh3p"]

_std_generics = {
    "hf_cc-pvdz_ae": (19, 19, 5, 5),
    "hf_cc-pvdz_fc": (19, 19, 5, 5),
    "bh3p_cc-pvdz_ae": (29, 29, 4, 3),
    "bh3p_cc-pvdz_fc": (29, 29, 4, 3),
    "h2o_aug-cc-pvdz_ae": (41, 41, 5, 5),
    "h2o_aug-cc-pvdz_fc": (41, 41, 5, 5),
    "nh2_aug-cc-pvdz_ae": (41, 41, 5, 4),
    "nh2_aug-cc-pvdz_fc": (41, 41, 5, 4),
    "h2o_cfour-qz2p_ae": (48, 48, 5, 5),
    "h2o_cfour-qz2p_fc": (48, 48, 5, 5),
    "nh2_cfour-qz2p_ae": (48, 48, 5, 4),
    "nh2_cfour-qz2p_fc": (48, 48, 5, 4),
}
_std_generics = {
    k: dict(zip(["N BASIS FUNCTIONS", "N MOLECULAR ORBITALS", "N ALPHA ELECTRONS", "N BETA ELECTRONS"], v))
    for k, v in _std_generics.items()
}

_scf_hf_dz_pk_rhf = -100.01941126902270
_scf_bh3p_dz_pk_uhf = -25.94513842869638
_scf_bh3p_dz_pk_rohf = -25.943614318546

_scf_hf_dz_df_rhf = -100.019400605629
_scf_bh3p_dz_df_uhf = -25.945130559147
_scf_bh3p_dz_df_rohf = -25.943606522029

_scf_hf_dz_cd_rhf = -100.01939270219628
_scf_bh3p_dz_cd_uhf = -25.94511891510799
_scf_bh3p_dz_cd_rohf = -25.943595251664313


_scf_h2o_qz2p_pk_rhf = -76.0627484601
_scf_nh2_qz2p_pk_uhf = -55.5893469688
_scf_nh2_qz2p_pk_rohf = -55.5847372601

_scf_h2o_qz2p_df_rhf = -76.06274142753659
_scf_nh2_qz2p_df_uhf = -55.58934323208402
_scf_nh2_qz2p_df_rohf = -55.58473319013903

_scf_h2o_qz2p_cd_rhf = -76.06277445978574
_scf_nh2_qz2p_cd_uhf = -55.58934916135869
_scf_nh2_qz2p_cd_rohf = -55.58473942870229


_scf_h2o_adz_pk_rhf = -76.0413815332
_scf_nh2_adz_pk_uhf = -55.57513805247548
_scf_nh2_adz_pk_rohf = -55.570724348574

_scf_h2o_adz_df_rhf = -76.04136132628614
_scf_nh2_adz_df_uhf = -55.57512538464817
_scf_nh2_adz_df_rohf = -55.57071142443952

_scf_h2o_adz_cd_rhf = -76.04132169763341
_scf_nh2_adz_cd_uhf = -55.57506886675886
_scf_nh2_adz_cd_rohf = -55.57065536578708


_std_suite = [
    # <<<  CONV-AE-CONV  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.203781911950,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05427697023782003,
            "MP2 TOTAL GRADIENT": np.array(
                [  # fnocc findif-5 ae pk+conv
                    0.0000000000,
                    0.0000000000,
                    0.0028193375,
                    0.0000000000,
                    0.0000000000,
                    -0.0028193375,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2099060277,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.048339903547,  # fnocc
            "LCCSD CORRELATION ENERGY": -0.2107436391,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "LCCSD SAME-SPIN CORRELATION ENERGY": -0.048460183760,  # fnocc
            "CCSD CORRELATION ENERGY": -0.208743643,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04857419039,
            "CCSD TOTAL GRADIENT": np.array([0.0, 0.0, 0.001989217717, 0.0, 0.0, -0.001989217717,]).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.0019363896542312043,
            "OLCCD REFERENCE CORRECTION ENERGY": 0.0005522939,  # p4n
            "OLCCD CORRELATION ENERGY": -0.2104417743,  # p4n
            "OLCCD SAME-SPIN CORRELATION ENERGY": -0.0484443079,  # occ
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.2218977246,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05669988343022163,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.009624481085, 0.0, 0.005505796371, -0.004812240542, 0.0, -0.005505796371, -0.004812240542,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2318870702,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.049937236558,  # fnocc
            "LCCSD CORRELATION ENERGY": -0.2341051403,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "LCCSD SAME-SPIN CORRELATION ENERGY": -0.050442387759,  # fnocc
            "CCSD CORRELATION ENERGY": -0.2294105794,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.050177977945205,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.007512595487, 0.0, 0.004613769715, -0.003756297743, 0.0, -0.004613769715, -0.003756297743,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00523856,
            "OLCCD REFERENCE CORRECTION ENERGY": 0.0011895155,  # p4n
            "OLCCD CORRELATION ENERGY": -0.2330452995,  # p4n
            "OLCCD SAME-SPIN CORRELATION ENERGY": -0.0503175223,  # occ
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.2701916672,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06530131,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, -0.000531535533, 0.0, -0.000960201925, 0.000265767766, 0.0, 0.000960201925, 0.000265767766,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2786913134,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.057792990490,  # fnocc
            "LCCSD CORRELATION ENERGY": -0.2808517417,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "LCCSD SAME-SPIN CORRELATION ENERGY": -0.058297242512,  # fnocc
            "CCSD CORRELATION ENERGY": -0.275705491773,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.058006927914493,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, -0.003374258422, 0.0, -0.002334452569, 0.001687129211, 0.0, 0.002334452569, 0.001687129211,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.007263596331,
            "OLCCD REFERENCE CORRECTION ENERGY": 0.0013521561,  # p4n
            "OLCCD CORRELATION ENERGY": -0.2800053174,  # p4n
            "OLCCD SAME-SPIN CORRELATION ENERGY": -0.0582676514,  # occ
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.05948928003552,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001918693775,
            "MP2 TOTAL GRADIENT": np.array(
                [  # occ findif-5 ae pk+conv
                    0.00000000000000,
                    0.00000000000000,
                    0.01250561195911,
                    0.00000000000000,
                    0.00000000000000,
                    -0.01206536529299,
                    0.00000000000000,
                    0.01033165380573,
                    -0.00022012333306,
                    0.00000000000000,
                    -0.01033165380573,
                    -0.00022012333306,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.0834347185,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0024003297,  # occ
            "LCCSD CORRELATION ENERGY": -0.0848110820,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.08217287869,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.002377557359,
            "CCSD TOTAL GRADIENT": np.array(
                [
                    0.0,
                    0.0,
                    0.005209606766,
                    0.0,
                    0.0,
                    -0.005071403517,
                    0.0,
                    0.014880198292,
                    -0.000069101625,
                    0.0,
                    -0.014880198292,
                    -0.000069101625,
                ]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00062614,
            "OLCCD REFERENCE CORRECTION ENERGY": 0.0014842084,  # p4n
            "OLCCD CORRELATION ENERGY": -0.0847413506,  # p4n
            "OLCCD SAME-SPIN CORRELATION ENERGY": -0.0024486744,  # occ
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.15485993330517828,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03520162545964887,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.025490652204, 0.0, 0.013491755791, -0.012745326102, 0.0, -0.013491755791, -0.012745326102,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1770086091,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0341268118,  # occ
            "LCCSD CORRELATION ENERGY": -0.1786081472,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.17387203707017695,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.033935818857082,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.029278727285, 0.0, 0.015813927533, -0.014639363642, 0.0, -0.015813927533, -0.014639363642,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00384378,
            "OLCCD REFERENCE CORRECTION ENERGY": 0.0011118724,  # p4n
            "OLCCD CORRELATION ENERGY": -0.1781057943,  # p4n
            "OLCCD SAME-SPIN CORRELATION ENERGY": -0.0344689234,  # occ
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.195530391293,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04161633,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.013731673196, 0.0, 0.005352105826, -0.006865836598, 0.0, -0.005352105826, -0.006865836598,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2167878305,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0401306050,  # occ
            "LCCSD CORRELATION ENERGY": -0.2185061347,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.213298055172,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.039907245914335,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.016842165003, 0.0, 0.007150136873, -0.008421082502, 0.0, -0.007150136873, -0.008421082502,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00516659,
            "OLCCD REFERENCE CORRECTION ENERGY": 0.0012856903,  # p4n
            "OLCCD CORRELATION ENERGY": -0.2180560836,  # p4n
            "OLCCD SAME-SPIN CORRELATION ENERGY": -0.0405122800,  # occ
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.060478115157,
            "MP2 SINGLES ENERGY": -0.000694049865,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.002004909679,
            "MP2 TOTAL GRADIENT": np.array(
                [
                    # switches sign from unkn ref
                    0.000000000000000,
                    0.000000000000000,
                    0.013594741747853,
                    0.000000000000000,
                    0.000000000000000,
                    -0.013127629532095,
                    0.000000000000000,
                    0.010308255599051,
                    -0.000233556107879,
                    0.000000000000000,
                    -0.010308255599051,
                    -0.000233556107879,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.0834094914,  # p4n
            "LCCSD CORRELATION ENERGY": -0.0861427228,  # p4n
            "CCSD CORRELATION ENERGY": -0.08357160616,
            "CCSD SINGLES ENERGY": -0.0011743271,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.00244892164,
            "CCSD TOTAL GRADIENT": np.array(
                [
                    0.0,
                    0.0,
                    0.005568141758,
                    0.0,
                    0.0,
                    -0.005430974166,
                    0.0,
                    0.014884143028,
                    -0.000068583796,
                    0.0,
                    -0.014884143028,
                    -0.000068583796,
                ]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.000713766189,
            "OLCCD REFERENCE CORRECTION ENERGY": -0.0000399018,  # p4n
            "OLCCD CORRELATION ENERGY": -0.0862654609,  # p4n
            "OLCCD SAME-SPIN CORRELATION ENERGY": -0.0024486744,  # occ
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.15949744108346664,
            "MP2 SINGLES ENERGY": -0.0028296307982793997,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03541709278508698,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.025609525826, 0.0, 0.013506941035, -0.012804762913, 0.0, -0.013506941035, -0.012804762913,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1791714105,  # p4n
            "LCCSD CORRELATION ENERGY": -0.1830545845,  # p4n
            "CCSD CORRELATION ENERGY": -0.178236032911,
            "CCSD SINGLES ENERGY": -0.00327524740575,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.033982707798170,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.029273628227, 0.0, 0.015808308241, -0.014636814114, 0.0, -0.015808308241, -0.014636814114,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.003901085777,
            "OLCCD REFERENCE CORRECTION ENERGY": -0.0033018315,  # p4n
            "OLCCD CORRELATION ENERGY": -0.1825194982,  # p4n
            "OLCCD SAME-SPIN CORRELATION ENERGY": -0.0344689234,  # occ
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.2005395272,
            "MP2 SINGLES ENERGY": -0.00298375,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04178535,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.0138883429, 0.0, 0.005389090661, -0.00694417145, 0.0, -0.005389090661, -0.00694417145,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2191039411,  # p4n
            "LCCSD CORRELATION ENERGY": -0.2231241199,  # p4n
            "CCSD CORRELATION ENERGY": -0.217849506326,
            "CCSD SINGLES ENERGY": -0.00338286103325,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.039891470497466,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.016833254665, 0.0, 0.007144029475, -0.008416627332, 0.0, -0.007144029475, -0.008416627332,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.005233938447,
            "OLCCD REFERENCE CORRECTION ENERGY": -0.0033240178,  # p4n
            "OLCCD CORRELATION ENERGY": -0.2226657917,  # p4n
            "OLCCD SAME-SPIN CORRELATION ENERGY": -0.0405122800,  # occ
        },
    },
    # <<<  CONV-FC-CONV  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.201627516796,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0534895025483176,
            "MP2 TOTAL GRADIENT": np.array(
                [  # fnocc findif-5 fc pk+conv
                    0.00000000000000,
                    0.00000000000000,
                    0.00317450456474,
                    0.00000000000000,
                    0.00000000000000,
                    -0.00317450456474,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2079585027,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.047635656759,  # fnocc
            "LCCSD CORRELATION ENERGY": -0.2087915976,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "LCCSD SAME-SPIN CORRELATION ENERGY": -0.047754723454,  # fnocc
            "CCSD CORRELATION ENERGY": -0.2068152041,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.0478712079,
            "CCSD TOTAL GRADIENT": np.array([0.0, 0.0, 0.002335204281, 0.0, 0.0, -0.002335204281,]).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.0019205007159748158,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.2194081478,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.055833980855745646,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.010245839621, 0.0, 0.005893268945, -0.00512291981, 0.0, -0.005893268945, -0.00512291981,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2296135965,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.049154543318,  # fnocc
            "LCCSD CORRELATION ENERGY": -0.2318316308,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "LCCSD SAME-SPIN CORRELATION ENERGY": -0.049659952324,  # fnocc
            "CCSD CORRELATION ENERGY": -0.2271733460,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.049398348010672,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.008118157882, 0.0, 0.004988381189, -0.004059078941, 0.0, -0.004988381189, -0.004059078941,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00521238,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.24515185206,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06126410,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.00033347691, 0.0, -0.00056224437, -0.000166738455, 0.0, 0.00056224437, -0.000166738455,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2531942099,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.053842594884,  # fnocc
            "LCCSD CORRELATION ENERGY": -0.2553008820,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "LCCSD SAME-SPIN CORRELATION ENERGY": -0.054321637599,  # fnocc
            "CCSD CORRELATION ENERGY": -0.250330548844,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.054051928864870,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, -0.002486174824, 0.0, -0.001923330621, 0.001243087412, 0.0, 0.001923330621, 0.001243087412,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.007096579721,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.058423513790,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001767468898,
            "MP2 TOTAL GRADIENT": np.array(
                # switched sign from unkn origin
                [
                    0.000000000000000,
                    0.000000000000000,
                    0.012305278627642,
                    0.000000000000000,
                    0.000000000000000,
                    -0.011851332672482,
                    0.000000000000000,
                    0.010327045553422,
                    -0.000226972977580,
                    0.000000000000000,
                    -0.010327045553422,
                    -0.000226972977580,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.0824313452,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0022532590,  # occ
            "LCCSD CORRELATION ENERGY": -0.0837903430,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.08117105566,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.002231965267,
            "CCSD TOTAL GRADIENT": np.array(
                [
                    0.0,
                    0.0,
                    0.00496423512,
                    0.0,
                    0.0,
                    -0.004814203262,
                    0.0,
                    0.014877060204,
                    -0.000075015929,
                    0.0,
                    -0.014877060204,
                    -0.000075015929,
                ]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00060401,
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.15242755400188052,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03445360441348938,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.026279427993, 0.0, 0.013998590506, -0.013139713997, 0.0, -0.013998590506, -0.013139713997,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1747537294,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0334378864,  # occ
            "LCCSD CORRELATION ENERGY": -0.1763496376,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.1716495276680232,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.033248190929062,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.030055915902, 0.0, 0.016307167756, -0.015027957951, 0.0, -0.016307167756, -0.015027957951,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00381116,
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.171184123093,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03822454,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.014740098324, 0.0, 0.005852228009, -0.007370049162, 0.0, -0.005852228009, -0.007370049162,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1917024115,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0367596656,  # occ
            "LCCSD CORRELATION ENERGY": -0.1933416962,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.188317222733,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.036526852874970,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.017883390799, 0.0, 0.00765987541, -0.0089416954, 0.0, -0.00765987541, -0.0089416954,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00498265,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.059407254257,
            "MP2 SINGLES ENERGY": -0.000688368657,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001851937488,
            "MP2 TOTAL GRADIENT": np.array(
                [
                    # switched sign from unkn ref
                    0.000000000000000,
                    0.000000000000000,
                    0.013388410166131,
                    0.000000000000000,
                    0.000000000000000,
                    -0.012907368096590,
                    0.000000000000000,
                    0.010303507439169,
                    -0.000240521034770,
                    0.000000000000000,
                    -0.010303507439169,
                    -0.000240521034770,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.0824056198,  # p4n
            "LCCSD CORRELATION ENERGY": -0.0851177481,  # p4n
            "CCSD CORRELATION ENERGY": -0.08256719,
            "CCSD SINGLES ENERGY": -0.00117001688,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.00230304,
            "CCSD TOTAL GRADIENT": np.array(
                [
                    0.0,
                    0.0,
                    0.005323074361,
                    0.0,
                    0.0,
                    -0.005174249172,
                    0.0,
                    0.014881203442,
                    -0.000074412594,
                    0.0,
                    -0.014881203442,
                    -0.000074412594,
                ]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00068823,  # cfour only
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.15702660833165538,
            "MP2 SINGLES ENERGY": -0.0028059971624814647,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03466304269235235,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.026398091851, 0.0, 0.014012163884, -0.013199045925, 0.0, -0.014012163884, -0.013199045925,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1769020687,  # p4n
            "LCCSD CORRELATION ENERGY": -0.1807707740,  # p4n
            "CCSD CORRELATION ENERGY": -0.175988485854028,
            "CCSD SINGLES ENERGY": -0.003256808469230,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.033291143258924,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.030051791297, 0.0, 0.016301545337, -0.015025895649, 0.0, -0.016301545337, -0.015025895649,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.003863167899,  # cfour only
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.1761163066,
            "MP2 SINGLES ENERGY": -0.00294339,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03837483,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.014894057335, 0.0, 0.005886660707, -0.007447028667, 0.0, -0.005886660707, -0.007447028667,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1939920915,  # p4n
            "LCCSD CORRELATION ENERGY": -0.1979175937,  # p4n
            "CCSD CORRELATION ENERGY": -0.19282621471297376,
            "CCSD SINGLES ENERGY": -0.003354603508621,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.036502859698546,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.017873897449, 0.0, 0.007653541045, -0.008936948724, 0.0, -0.007653541045, -0.008936948724,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00504351,  # cfour only
        },
    },
    # <<<  CONV-AE-CD  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.203778449,
            "MP2 SINGLES ENERGY": 0.0,
            "LCCD CORRELATION ENERGY": -0.20990784,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.20874537,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00193646,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.22188844,
            "MP2 SINGLES ENERGY": 0.0,
            "LCCD CORRELATION ENERGY": -0.23188996,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.22941330,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00523874,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.27018509,
            "MP2 SINGLES ENERGY": 0.0,
            "LCCD CORRELATION ENERGY": -0.27869144,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.27570541,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00726403,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.059477703268,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001918940186,
            "LCCD CORRELATION ENERGY": -0.08343267,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.00240067,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.15485159,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03520588,
            "LCCD CORRELATION ENERGY": -0.17701281,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.03413088,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.19552518,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04162160,
            "LCCD CORRELATION ENERGY": -0.21678793,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.04013550,  # dfocc
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.0604664810,
            "MP2 SINGLES ENERGY": -0.000694049858,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.002005152902,
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.15948893,
            "MP2 SINGLES ENERGY": -0.00282963,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03542136,
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.20053428,
            "MP2 SINGLES ENERGY": -0.00298375,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04179065,
        },
    },
    # <<<  CONV-FC-CD  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.20162439774,
            "MP2 SINGLES ENERGY": 0.0,
            "LCCD CORRELATION ENERGY": -0.20796060,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.20681721,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00192057,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.21939933,
            "MP2 SINGLES ENERGY": 0.0,
            "LCCD CORRELATION ENERGY": -0.22961687,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.22717646,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00521255,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.24514540,
            "MP2 SINGLES ENERGY": 0.0,
            "LCCD CORRELATION ENERGY": -0.25319438,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "CCSD CORRELATION ENERGY": -0.25033052,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00709694,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.05841222894,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0017676971,
            "LCCD CORRELATION ENERGY": -0.08242955,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.00225358,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.15241971,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03445776,
            "LCCD CORRELATION ENERGY": -0.17475833,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.03344184,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.17117906,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03822980,
            "LCCD CORRELATION ENERGY": -0.19170259,  # dfocc
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.03676455,  # dfocc
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.059395907176,
            "MP2 SINGLES ENERGY": -0.00068836865,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001852162877,
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.15701860,
            "MP2 SINGLES ENERGY": -0.00280600,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03466721,
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.17611121,
            "MP2 SINGLES ENERGY": -0.00294339,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03838011,
        },
    },
    # <<<  CONV-AE-DF  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.20377997248921056,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05431321036920538,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.0027998, 0.0, 0.0, -0.0027998]  # dfmp2 findif-5 ae pk+df
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2100497124,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.20888438,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00193859,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.22188894,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05674808,
            "LCCD CORRELATION ENERGY": -0.2320261414,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.22954333,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00524393,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.27018057,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06530212,
            "LCCD CORRELATION ENERGY": -0.2786878429,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.27570207,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00726375,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.05945820694747983,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0019203155958724552,
            "LCCD CORRELATION ENERGY": -0.0835080983,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0024018298,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.15484736,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03523329,
            "LCCD CORRELATION ENERGY": -0.1771107929,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0340809591,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.19551918,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04161696,
            "LCCD CORRELATION ENERGY": -0.2167841215,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0401306929,  # dfocc
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.0604460449537298,
            "MP2 SINGLES ENERGY": -0.0006940498589629459,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0020066877639503184,
            "LCCD CORRELATION ENERGY": -0.0834825821,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.15948289,
            "MP2 SINGLES ENERGY": -0.00282963,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03544835,
            "LCCD CORRELATION ENERGY": -0.1792713801,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.20052829,
            "MP2 SINGLES ENERGY": -0.00298375,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04178599,
            "LCCD CORRELATION ENERGY": -0.2191002183,  # p4n
        },
    },
    # <<<  CONV-FC-DF  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.20162566806258586,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05352569481658172,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.00315485, 0.0, 0.0, -0.00315485]  # dfmp2 findif-5 fc pk+df
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2081020566,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.20695586,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00192267,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.21939942,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05588210,
            "LCCD CORRELATION ENERGY": -0.2297524911,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.22730597,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00521769,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_pk_rhf,
            "MP2 CORRELATION ENERGY": -0.24514425,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06126481,
            "LCCD CORRELATION ENERGY": -0.2531939249,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.25033030,
            "CCSD SINGLES ENERGY": 0.0,
            "(T) CORRECTION ENERGY": -0.00709666,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.058392397606538686,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0017690135626491292,
            "LCCD CORRELATION ENERGY": -0.0825046579,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0022547041,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.15241501,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03448519,
            "LCCD CORRELATION ENERGY": -0.1748557523,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0333918420,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_uhf,
            "MP2 CORRELATION ENERGY": -0.17117615,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03822512,
            "LCCD CORRELATION ENERGY": -0.1917015960,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0367596684,  # dfocc
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.05937514348825628,
            "MP2 SINGLES ENERGY": -0.0006883686516107368,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0018536363586657242,
            "LCCD CORRELATION ENERGY": -0.0824786458,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.15701209,
            "MP2 SINGLES ENERGY": -0.00280600,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03469422,
            "LCCD CORRELATION ENERGY": -0.1770018748,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "pk",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_pk_rohf,
            "MP2 CORRELATION ENERGY": -0.17610830,
            "MP2 SINGLES ENERGY": -0.00294339,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03837541,
            "LCCD CORRELATION ENERGY": -0.1939912613,  # p4n
        },
    },
    # <<<  CD-AE-CD  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.20377328786815951,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05427570,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 ae cd+cd
                [0.0, 0.0, 0.00281146, 0.0, 0.0, -0.00281146]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.20990226,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.20873986012771106,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04857381,
            "(T) CORRECTION ENERGY": -0.0019363109218456449,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.22188817,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05670210,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 ae cd+cd
                [0.0, 0.0, 0.009643414073, 0.0, 0.005501440694, -0.004821707036, 0.0, -0.005501440694, -0.004821707036,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.23188949,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.22941290,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.05017955,
            "(T) CORRECTION ENERGY": -0.00523867,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.27018399,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06530655,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 ae cd+cd
                [0.0, 0.0, -0.000546229785, 0.0, -0.000967320028, 0.000273114892, 0.0, 0.000967320028, 0.000273114892,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.27869015,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.27570421,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.05801141,
            "(T) CORRECTION ENERGY": -0.00726395,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.059476326350818454,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0019188791023,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 ae cd+cd
                [
                    0.0,
                    0.0,
                    0.0125029,
                    0.0,
                    0.0,
                    -0.01205882,
                    0.0,
                    0.01033888,
                    -0.00022204,
                    0.0,
                    -0.01033888,
                    -0.00022204,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.08343038,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.00240059,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.15485101,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03520580,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfocc findif-5 ae cd+cd
                    0.0,
                    0.0,
                    0.025470063809,
                    0.0,
                    0.013535107677,
                    -0.012735031905,
                    0.0,
                    -0.013535107677,
                    -0.012735031905,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.17701192,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.03413070,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.19552441,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04162127,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfocc findif-5 ae cd+cd
                    0.0,
                    0.0,
                    0.013727424376,
                    0.0,
                    0.005348487843,
                    -0.006863712188,
                    0.0,
                    -0.005348487843,
                    -0.006863712188,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.21678706,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.04013515,  # dfocc
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.06046475293245379,
            "MP2 SINGLES ENERGY": -0.00069387098844,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.002005071400,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 ae cd+cd
                [0.0, 0.0, 0.01359215, 0.0, 0.0, -0.01312116, 0.0, 0.01031541, -0.0002355, 0.0, -0.01031541, -0.0002355]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.15948823,
            "MP2 SINGLES ENERGY": -0.00282948,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03542128,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfocc findif-5 ae cd+cd
                    0.0,
                    0.0,
                    0.025588961002,
                    0.0,
                    0.013550360249,
                    -0.012794480501,
                    0.0,
                    -0.013550360249,
                    -0.012794480501,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.20053352,
            "MP2 SINGLES ENERGY": -0.00298373,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04179032,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfocc findif-5 ae cd+cd
                    0.0,
                    0.0,
                    0.013884053665,
                    0.0,
                    0.005385412795,
                    -0.006942026833,
                    0.0,
                    -0.005385412795,
                    -0.006942026833,
                ]
            ).reshape((-1, 3)),
        },
    },
    # <<<  CD-FC-CD  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.201619244596,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05348825,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 fc cd+cd
                [0.0, 0.0, 0.00316665, 0.0, 0.0, -0.00316665]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.20795503,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.2068117080298787,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04787083,
            "(T) CORRECTION ENERGY": -0.0019204203743072874,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.21939907,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05583617,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 fc cd+cd
                [0.0, 0.0, 0.010264703011, 0.0, 0.00588885358, -0.005132351506, 0.0, -0.00588885358, -0.005132351506,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.22961642,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.22717607,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04939986,
            "(T) CORRECTION ENERGY": -0.00521248,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.24514436,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06126931,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 fc cd+cd
                [0.0, 0.0, 0.000318778691, 0.0, -0.000569356625, -0.000159389346, 0.0, 0.000569356625, -0.000159389346,]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.25319315,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.25032939,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.05405638,
            "(T) CORRECTION ENERGY": -0.00709686,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.058410863785614,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001767642489,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 fc cd+cd
                [
                    0.0,
                    0.0,
                    0.01230315,
                    0.0,
                    0.0,
                    -0.01184537,
                    0.0,
                    0.01033427,
                    -0.00022889,
                    0.0,
                    -0.01033427,
                    -0.00022889,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.08242726,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.00225350,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.15241915,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03445770,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfocc findif-5 fc cd+cd
                    0.0,
                    0.0,
                    0.026258239074,
                    0.0,
                    0.01404196652,
                    -0.013129119537,
                    0.0,
                    -0.01404196652,
                    -0.013129119537,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.17475747,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0334416820,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.17117831,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03822948,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfocc findif-5 fc cd+cd
                    0.0,
                    0.0,
                    0.014735846129,
                    0.0,
                    0.005848618964,
                    -0.007367923065,
                    0.0,
                    -0.005848618964,
                    -0.007367923065,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.19170174,  # dfocc
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.03676422,  # dfocc
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.05939419492939635,
            "MP2 SINGLES ENERGY": -0.0006881934,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0018520881544,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 fc cd+cd
                [
                    0.0,
                    0.0,
                    0.01338641,
                    0.0,
                    0.0,
                    -0.01290149,
                    0.0,
                    0.01031066,
                    -0.00024246,
                    0.0,
                    -0.01031066,
                    -0.00024246,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.15701792,
            "MP2 SINGLES ENERGY": -0.00280584,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03466715,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfocc findif-5 fc cd+cd
                    0.0,
                    0.0,
                    0.026376923581,
                    0.0,
                    0.014055606253,
                    -0.01318846179,
                    0.0,
                    -0.014055606253,
                    -0.01318846179,
                ]
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "cd",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.17611046,
            "MP2 SINGLES ENERGY": -0.00294336,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03837979,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfocc findif-5 fc cd+cd
                    0.0,
                    0.0,
                    0.014889762324,
                    0.0,
                    0.00588299146,
                    -0.007444881162,
                    0.0,
                    -0.00588299146,
                    -0.007444881162,
                ]
            ).reshape((-1, 3)),
        },
    },
    # <<<  CD-AE-DF  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.2037748110768,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.054311870576,
            "MP2 TOTAL GRADIENT": np.array(
                # dfmp2 findif-5 ae cd+df
                [0.0, 0.0, 0.00279182, 0.0, 0.0, -0.00279182]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2100441271,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.20887885,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04845784,
            "(T) CORRECTION ENERGY": -0.00193844,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.22188866,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05674801,
            "LCCD CORRELATION ENERGY": -0.2320256729,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.22954292,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.05010092,
            "(T) CORRECTION ENERGY": -0.00524386,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.27017947,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06530177,
            "LCCD CORRELATION ENERGY": -0.2786865554,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.27570087,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.05800702,
            "(T) CORRECTION ENERGY": -0.00726367,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.059456828193,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.00192025457659,
            "LCCD CORRELATION ENERGY": -0.0835057932,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0024017496,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.15484678,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03523320,
            "LCCD CORRELATION ENERGY": -0.1771099018,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0340807883,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.19551841,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04161663,
            "LCCD CORRELATION ENERGY": -0.2167832515,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0401303480,  # dfocc
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.06044431529,
            "MP2 SINGLES ENERGY": -0.00069387098844,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0020066063,
            "LCCD CORRELATION ENERGY": -0.0834800819,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.15948219,
            "MP2 SINGLES ENERGY": -0.00282948,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03544828,
            "LCCD CORRELATION ENERGY": -0.1792705171,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.20052752,
            "MP2 SINGLES ENERGY": -0.00298373,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04178566,
            "LCCD CORRELATION ENERGY": -0.2190993784,  # p4n
        },
    },
    # <<<  CD-FC-DF  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.2016205147678,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0535243575,
            "MP2 TOTAL GRADIENT": np.array(
                # dfmp2 findif-5 fc cd+df
                [0.0, 0.0, 0.00314686, 0.0, 0.0, -0.00314686]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2080964757,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.20695033,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04775464,
            "(T) CORRECTION ENERGY": -0.00192252,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.21939916,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05588204,
            "LCCD CORRELATION ENERGY": -0.2297520405,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.22730558,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04932106,
            "(T) CORRECTION ENERGY": -0.00521762,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.24514320,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06126448,
            "LCCD CORRELATION ENERGY": -0.2531926943,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.25032917,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.05405189,
            "(T) CORRECTION ENERGY": -0.00709658,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.05839103061,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.00176895897,
            "LCCD CORRELATION ENERGY": -0.0825023638,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0022546311,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.15241445,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03448512,
            "LCCD CORRELATION ENERGY": -0.1748548876,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0333916888,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_cd_uhf,
            "MP2 CORRELATION ENERGY": -0.17117540,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03822480,
            "LCCD CORRELATION ENERGY": -0.1917007514,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0367593319,  # dfocc
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.05937342969795,
            "MP2 SINGLES ENERGY": -0.0006881934,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001853561678,
            "LCCD CORRELATION ENERGY": -0.0824761581,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.15701141,
            "MP2 SINGLES ENERGY": -0.00280584,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03469416,
            "LCCD CORRELATION ENERGY": -0.1770010376,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "cd",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_cd_rohf,
            "MP2 CORRELATION ENERGY": -0.17610756,
            "MP2 SINGLES ENERGY": -0.00294336,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03837509,
            "LCCD CORRELATION ENERGY": -0.1939904460,  # p4n
        },
    },
    # <<<  DF-AE-DF  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.2037649370559149,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05430875283333263,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.00279211492833,
                    0.00000000000000,
                    0.00000000000000,
                    -0.00279211492833,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2100337333,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.20886884012911314,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04845491,
            "CCSD TOTAL GRADIENT": np.array([0.0, 0.0, 0.001970675302, 0.0, 0.0, -0.001970675302,]).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.0019380186429220421,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.22187976,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05674571,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.0,
                    0.0,
                    0.00962182765,
                    0.0,
                    0.005498317937,
                    -0.004810913825,
                    0.0,
                    -0.005498317937,
                    -0.004810913825,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2320149229,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.22953289,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.05009877,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.007518759967, 0.0, 0.004613106602, -0.003759379983, 0.0, -0.004613106602, -0.003759379983,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00524345,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_df_rhf,
            "MP2 CORRELATION ENERGY": -0.27016105,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06529808,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.0,
                    0.0,
                    -0.000566657943,
                    0.0,
                    -0.000968877215,
                    0.000283328971,
                    0.0,
                    0.000968877215,
                    0.000283328971,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2786671617,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.27568236,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.05800380,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, -0.003408844165, 0.0, -0.002343169064, 0.001704422083, 0.0, 0.002343169064, 0.001704422083]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00726213,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
            "MP2 CORRELATION ENERGY": -0.0594557966607590,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001920220330437888,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.01252024755551,
                    0.00000000000000,
                    0.00000000000000,
                    -0.01207773525598,
                    0.00000000000000,
                    0.01032204616770,
                    -0.00022125614977,
                    0.00000000000000,
                    -0.01032204616770,
                    -0.00022125614977,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.0835030877,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0024016379,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_df_uhf,
            "MP2 CORRELATION ENERGY": -0.15483909,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03523134,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.0,
                    0.0,
                    0.025476049585,
                    0.0,
                    0.013480567736,
                    -0.012738024793,
                    0.0,
                    -0.013480567736,
                    -0.012738024793,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1770997033,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0340788149,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_df_uhf,
            "MP2 CORRELATION ENERGY": -0.19550726,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04161470,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.0,
                    0.0,
                    0.013708831104,
                    0.0,
                    0.005340400162,
                    -0.006854415552,
                    0.0,
                    -0.005340400162,
                    -0.006854415552,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2167706529,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0401283617,  # dfocc
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
            "MP2 CORRELATION ENERGY": -0.0604436327328384,
            "MP2 SINGLES ENERGY": -0.0006940750313001934,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0020065676314032863,
            "MP2 TOTAL GRADIENT": np.array(
                [  # occ findif-5 ae df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.01361287313486,
                    0.00000000000000,
                    0.00000000000000,
                    -0.01314329502424,
                    0.00000000000000,
                    0.01029838165151,
                    -0.00023478905531,
                    0.00000000000000,
                    -0.01029838165151,
                    -0.00023478905531,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.0834776542,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_df_rohf,
            "MP2 CORRELATION ENERGY": -0.15947485,
            "MP2 SINGLES ENERGY": -0.00282982,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03544639,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.0,
                    0.0,
                    0.025593521597,
                    0.0,
                    0.013495283342,
                    -0.012796760798,
                    0.0,
                    -0.013495283342,
                    -0.012796760798,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1792603912,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "rohf",
            "fcae": "ae",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_df_rohf,
            "MP2 CORRELATION ENERGY": -0.20051655,
            "MP2 SINGLES ENERGY": -0.00298400,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.04178365,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 ae df+df
                    0.0,
                    0.0,
                    0.013865245912,
                    0.0,
                    0.005377216253,
                    -0.006932622956,
                    0.0,
                    -0.005377216253,
                    -0.006932622956,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2190866990,  # p4n
        },
    },
    # <<<  DF-FC-DF  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.201610660387,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0535212487451535,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.00314716362539,
                    0.00000000000000,
                    0.00000000000000,
                    -0.00314716362539,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2080860831,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.20694032546082639,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04775171,
            "CCSD TOTAL GRADIENT": np.array([0.0, 0.0, 0.002316563628, 0.0, 0.0, -0.002316563628,]).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.001922093564526723,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_adz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.21939028,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05587974,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.0,
                    0.0,
                    0.010243193827,
                    0.0,
                    0.005885789424,
                    -0.005121596913,
                    0.0,
                    -0.005885789424,
                    -0.005121596913,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2297412879,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.22729554,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04931891,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.008124347934, 0.0, 0.004987676555, -0.004062173967, 0.0, -0.004987676555, -0.004062173967,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00521721,
        },
    },
    {
        "meta": {
            "system": "h2o",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_h2o_qz2p_df_rhf,
            "MP2 CORRELATION ENERGY": -0.24512893,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.06126089,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.0,
                    0.0,
                    0.000298272081,
                    0.0,
                    -0.000570968013,
                    -0.00014913604,
                    0.0,
                    0.000570968013,
                    -0.00014913604,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2531777549,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "CCSD CORRELATION ENERGY": -0.25031508,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.05404876,
            "CCSD TOTAL GRADIENT": np.array(
                [0.0, 0.0, -0.002520920562, 0.0, -0.001932133533, 0.001260460281, 0.0, 0.001932133533, 0.001260460281,]
            ).reshape((-1, 3)),
            "(T) CORRECTION ENERGY": -0.00709505,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
            "MP2 CORRELATION ENERGY": -0.058390006825,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001768919072594215,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.00000000000000,
                    0.00000000000000,
                    0.01231996225662,
                    0.00000000000000,
                    0.00000000000000,
                    -0.01186374280678,
                    0.00000000000000,
                    0.01031743020277,
                    -0.00022810972492,
                    0.00000000000000,
                    -0.01031743020277,
                    -0.00022810972492,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.0824996438,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0022545103,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_df_uhf,
            "MP2 CORRELATION ENERGY": -0.15240678,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03448325,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.0,
                    0.0,
                    0.026264866471,
                    0.0,
                    0.013987430104,
                    -0.013132433236,
                    0.0,
                    -0.013987430104,
                    -0.013132433236,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1748446809,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0333897039,  # dfocc
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "uhf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_df_uhf,
            "MP2 CORRELATION ENERGY": -0.17116675,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03822296,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.0,
                    0.0,
                    0.01471721142,
                    0.0,
                    0.005840479593,
                    -0.00735860571,
                    0.0,
                    -0.005840479593,
                    -0.00735860571,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1916908596,  # p4n
            "LCCD SINGLES ENERGY": 0.0000000000,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.0367574293,
        },
    },
    {
        "meta": {
            "system": "bh3p",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
            "MP2 CORRELATION ENERGY": -0.059372748391,
            "MP2 SINGLES ENERGY": -0.000688391888527046,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0018535174789756292,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.0,
                    0.0,
                    0.01340658,
                    0.0,
                    0.0,
                    -0.01292306,
                    0.0,
                    0.01029363,
                    -0.00024176,
                    0.0,
                    -0.01029363,
                    -0.00024176,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.0824737155,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "aug-cc-pvdz",
            "scf_type": "df",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_adz_df_rohf,
            "MP2 CORRELATION ENERGY": -0.15700408,
            "MP2 SINGLES ENERGY": -0.00280619,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03469227,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.0,
                    0.0,
                    0.026382129796,
                    0.0,
                    0.014000533629,
                    -0.013191064898,
                    0.0,
                    -0.014000533629,
                    -0.013191064898,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1769909051,  # p4n
        },
    },
    {
        "meta": {
            "system": "nh2",
            "basis": "cfour-qz2p",
            "scf_type": "df",
            "reference": "rohf",
            "fcae": "fc",
            "corl_type": "df",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_nh2_qz2p_df_rohf,
            "MP2 CORRELATION ENERGY": -0.17609909,
            "MP2 SINGLES ENERGY": -0.00294363,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.03837317,
            "MP2 TOTAL GRADIENT": np.array(
                [  # dfmp2 findif-5 fc df+df
                    0.0,
                    0.0,
                    0.014870916178,
                    0.0,
                    0.00587474124,
                    -0.007435458089,
                    0.0,
                    -0.00587474124,
                    -0.007435458089,
                ]
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.1939804718,  # p4n
        },
    },
    # <<<  lopsided SCF/CORL algorithms  >>>
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.201612517228,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05348507322421174,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.00316682, 0.0, 0.0, -0.00316682]  # occ findif-5 fc df+conv
            ).reshape((-1, 3)),
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "fc",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.20162236483,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.053488165399,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.00316652, 0.0, 0.0, -0.00316652]  # occ findif-5 fc cd+conv
            ).reshape((-1, 3)),
        },
    },
    # {
    #    "meta": {
    #        "system": "hf",
    #        "basis": "cc-pvdz",
    #        "scf_type": "df",
    #        "reference": "rhf",
    #        "fcae": "fc",
    #        "corl_type": "cd",
    #    },
    #    "data": {
    #        "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
    #        "MP2 CORRELATION ENERGY": -0.201609396752,
    #        #            "MP2 TOTAL ENERGY": -100.221010002381,
    #        "MP2 SINGLES ENERGY": 0.0,
    #        "MP2 SAME-SPIN CORRELATION ENERGY": -0.4,
    #    },
    # },
    #    {
    #        "meta": {
    #            "system": "bh3p",
    #            "basis": "cc-pvdz",
    #            "scf_type": "df",
    #            "reference": "uhf",
    #            "fcae": "fc",
    #            "corl_type": "conv",
    #        },
    #        "data": {
    #            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
    #            "MP2 CORRELATION ENERGY": -0.058421122206,
    #            #            "MP2 TOTAL ENERGY": -26.003551681354,
    #            "MP2 SINGLES ENERGY": 0.0,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -0.5,
    #        },
    #    },
    #    {
    #        "meta": {
    #            "system": "bh3p",
    #            "basis": "cc-pvdz",
    #            "scf_type": "df",
    #            "reference": "uhf",
    #            "fcae": "fc",
    #            "corl_type": "cd",
    #        },
    #        "data": {
    #            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
    #            "MP2 CORRELATION ENERGY": -0.058409837177,
    #            #            "MP2 TOTAL ENERGY": -26.003540396324,
    #            "MP2 SINGLES ENERGY": 0.0,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -0.7,
    #        },
    #    },
    #    {
    #        "meta": {
    #            "system": "bh3p",
    #            "basis": "cc-pvdz",
    #            "scf_type": "df",
    #            "reference": "rohf",
    #            "fcae": "fc",
    #            "corl_type": "conv",
    #        },
    #        "data": {
    #            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
    #            "MP2 CORRELATION ENERGY": -0.060939211739,
    #            #            "MP2 TOTAL ENERGY": -26.004545733768,
    #            "MP2 SINGLES ENERGY": 1.1,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -1.1,
    #        },
    #    },
    #    {
    #        "meta": {
    #            "system": "bh3p",
    #            "basis": "cc-pvdz",
    #            "scf_type": "df",
    #            "reference": "rohf",
    #            "fcae": "fc",
    #            "corl_type": "cd",
    #        },
    #        "data": {
    #            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
    #            "MP2 CORRELATION ENERGY": -0.059393510962,
    #            #            "MP2 TOTAL ENERGY": -26.003000032991,
    #            "MP2 SINGLES ENERGY": 1.3,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -1.3,
    #        },
    #    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "df",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
            "MP2 CORRELATION ENERGY": -0.2037668844651997,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05427252944164894,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.00281165, 0.0, 0.0, -0.00281165]  # occ findif-5 ae df+conv
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2098900858,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.048336089041,  # fnocc
            "LCCSD CORRELATION ENERGY": -0.2107275173,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "LCCSD SAME-SPIN CORRELATION ENERGY": -0.048456320034,  # fnocc
            "CCSD CORRELATION ENERGY": -0.20872812,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04857038,
        },
    },
    {
        "meta": {
            "system": "hf",
            "basis": "cc-pvdz",
            "scf_type": "cd",
            "reference": "rhf",
            "fcae": "ae",
            "corl_type": "conv",
        },
        "data": {
            "HF TOTAL ENERGY": _scf_hf_dz_cd_rhf,
            "MP2 CORRELATION ENERGY": -0.2037767503537,
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05427563053,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.00281136, 0.0, 0.0, -0.00281136]  # occ findif-5 ae cd+conv
            ).reshape((-1, 3)),
            "LCCD CORRELATION ENERGY": -0.2099004485,  # p4n
            "LCCD SINGLES ENERGY": 0.0,
            "LCCD SAME-SPIN CORRELATION ENERGY": -0.048339111990,  # fnocc
            "LCCSD CORRELATION ENERGY": -0.2107380019,  # p4n
            "LCCSD SINGLES ENERGY": 0.0,
            "LCCSD SAME-SPIN CORRELATION ENERGY": -0.048459381537,  # fnocc
            "CCSD CORRELATION ENERGY": -0.20873814,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04857333,
        },
    },
    #    {
    #        "meta": {
    #            "system": "hf",
    #            "basis": "cc-pvdz",
    #            "scf_type": "df",
    #            "reference": "rhf",
    #            "fcae": "ae",
    #            "corl_type": "cd",
    #        },
    #        "data": {
    #            "HF TOTAL ENERGY": _scf_hf_dz_df_rhf,
    #            "MP2 CORRELATION ENERGY": -2.3,
    #            "MP2 SINGLES ENERGY": 0.0,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.3,
    #        },
    #    },
    #    {
    #        "meta": {
    #            "system": "bh3p",
    #            "basis": "cc-pvdz",
    #            "scf_type": "df",
    #            "reference": "uhf",
    #            "fcae": "ae",
    #            "corl_type": "conv",
    #        },
    #        "data": {
    #            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
    #            "MP2 CORRELATION ENERGY": -2.4,
    #            "MP2 SINGLES ENERGY": 0.0,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.4,
    #        },
    #    },
    #    {
    #        "meta": {
    #            "system": "bh3p",
    #            "basis": "cc-pvdz",
    #            "scf_type": "df",
    #            "reference": "uhf",
    #            "fcae": "ae",
    #            "corl_type": "cd",
    #        },
    #        "data": {
    #            "HF TOTAL ENERGY": _scf_bh3p_dz_df_uhf,
    #            "MP2 CORRELATION ENERGY": -2.5,
    #            "MP2 SINGLES ENERGY": 0.0,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.5,
    #        },
    #    },
    #    {
    #        "meta": {
    #            "system": "bh3p",
    #            "basis": "cc-pvdz",
    #            "scf_type": "df",
    #            "reference": "rohf",
    #            "fcae": "ae",
    #            "corl_type": "conv",
    #        },
    #        "data": {
    #            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
    #            "MP2 CORRELATION ENERGY": -2.7,
    #            "MP2 SINGLES ENERGY": -2.7,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.7,
    #        },
    #    },
    #    {
    #        "meta": {
    #            "system": "bh3p",
    #            "basis": "cc-pvdz",
    #            "scf_type": "df",
    #            "reference": "rohf",
    #            "fcae": "ae",
    #            "corl_type": "cd",
    #        },
    #        "data": {
    #            "HF TOTAL ENERGY": _scf_bh3p_dz_df_rohf,
    #            "MP2 CORRELATION ENERGY": -2.8,
    #            "MP2 SINGLES ENERGY": -2.8,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.8,
    #        },
    #    },
]


for calc in _std_suite:
    if calc["data"]:
        if "MP2 CORRELATION ENERGY" in calc["data"]:
            calc["data"]["MP2 TOTAL ENERGY"] = calc["data"]["MP2 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            if "MP2 SINGLES ENERGY" in calc["data"]:
                calc["data"]["MP2 DOUBLES ENERGY"] = (
                    calc["data"]["MP2 CORRELATION ENERGY"] - calc["data"]["MP2 SINGLES ENERGY"]
                )
                if "MP2 SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["MP2 CORRELATION ENERGY"]
                        - calc["data"]["MP2 SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["MP2 SINGLES ENERGY"]
                    )
                    calc["data"]["SCS-MP2 CORRELATION ENERGY"] = (
                        (1 / 3) * calc["data"]["MP2 SAME-SPIN CORRELATION ENERGY"]
                        + (6 / 5) * calc["data"]["MP2 OPPOSITE-SPIN CORRELATION ENERGY"]
                        + calc["data"]["MP2 SINGLES ENERGY"]
                    )
                    calc["data"]["SCS-MP2 TOTAL ENERGY"] = (
                        calc["data"]["SCS-MP2 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
                    )

        if "LCCD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["LCCD TOTAL ENERGY"] = (
                calc["data"]["LCCD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            if "LCCD SINGLES ENERGY" in calc["data"]:
                calc["data"]["LCCD DOUBLES ENERGY"] = (
                    calc["data"]["LCCD CORRELATION ENERGY"] - calc["data"]["LCCD SINGLES ENERGY"]
                )
                if "LCCD SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["LCCD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["LCCD CORRELATION ENERGY"]
                        - calc["data"]["LCCD SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["LCCD SINGLES ENERGY"]
                    )

        if "LCCSD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["LCCSD TOTAL ENERGY"] = (
                calc["data"]["LCCSD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            if "LCCSD SINGLES ENERGY" in calc["data"]:
                calc["data"]["LCCSD DOUBLES ENERGY"] = (
                    calc["data"]["LCCSD CORRELATION ENERGY"] - calc["data"]["LCCSD SINGLES ENERGY"]
                )
                if "LCCSD SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["LCCSD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["LCCSD CORRELATION ENERGY"]
                        - calc["data"]["LCCSD SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["LCCSD SINGLES ENERGY"]
                    )

        if "CCSD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCSD TOTAL ENERGY"] = (
                calc["data"]["CCSD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            if "CCSD SINGLES ENERGY" in calc["data"]:
                calc["data"]["CCSD DOUBLES ENERGY"] = (
                    calc["data"]["CCSD CORRELATION ENERGY"] - calc["data"]["CCSD SINGLES ENERGY"]
                )
                if "CCSD SAME-SPIN CORRELATION ENERGY" in calc["data"]:
                    calc["data"]["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                        calc["data"]["CCSD CORRELATION ENERGY"]
                        - calc["data"]["CCSD SAME-SPIN CORRELATION ENERGY"]
                        - calc["data"]["CCSD SINGLES ENERGY"]
                    )

        if "(T) CORRECTION ENERGY" in calc["data"]:
            calc["data"]["CCSD(T) CORRELATION ENERGY"] = (
                calc["data"]["CCSD CORRELATION ENERGY"] + calc["data"]["(T) CORRECTION ENERGY"]
            )
            calc["data"]["CCSD(T) TOTAL ENERGY"] = (
                calc["data"]["CCSD(T) CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )

        if "OLCCD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["OLCCD TOTAL ENERGY"] = (
                calc["data"]["OLCCD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            calc["data"]["OLCCD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                calc["data"]["OLCCD CORRELATION ENERGY"]
                - calc["data"]["OLCCD REFERENCE CORRECTION ENERGY"]
                - calc["data"]["OLCCD SAME-SPIN CORRELATION ENERGY"]
            )

    calc["data"].update(_std_generics[f"{calc['meta']['system']}_{calc['meta']['basis']}_{calc['meta']['fcae']}"])


def answer_hash(**kwargs):
    system = kwargs.pop("system")
    basis = kwargs.pop("basis")
    scf_type = kwargs.pop("scf_type")
    reference = kwargs.pop("reference")
    fcae = kwargs.pop("fcae")
    corl_type = kwargs.pop("corl_type")

    return "_".join([system, basis, scf_type, reference, fcae, corl_type])


std_suite = {answer_hash(**calc["meta"]): calc["data"] for calc in _std_suite}
