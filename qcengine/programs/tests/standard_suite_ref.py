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

_scf_h2o_adz_pk_rhf = -76.0413815332
_scf_nh2_adz_pk_uhf = -55.57513805247548
_scf_nh2_adz_pk_rohf = -55.570724348574

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
            "CCSD CORRELATION ENERGY": -0.208743643,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.04857419039,
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
            "CCSD CORRELATION ENERGY": -0.08217287869,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.002377557359,
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
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.002004909679,
            "MP2 SINGLES ENERGY": -0.000694049865,
            "MP2 CORRELATION ENERGY": -0.060478115157,
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
            "CCSD CORRELATION ENERGY": -0.08357160616,
            "CCSD SINGLES ENERGY": -0.0011743271,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.00244892164,
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
            "CCSD CORRELATION ENERGY": -0.275705491773,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.058006927914493,
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
            "CCSD CORRELATION ENERGY": -0.2294105794,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.050177977945205,
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
            "CCSD CORRELATION ENERGY": -0.213298055172,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.039907245914335,
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
            "CCSD CORRELATION ENERGY": -0.17387203707017695,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.033935818857082,
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
            "CCSD CORRELATION ENERGY": -0.217849506326,
            "CCSD SINGLES ENERGY": -0.00338286103325,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.039891470497466,
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
            "CCSD CORRELATION ENERGY": -0.178236032911,
            "CCSD SINGLES ENERGY": -0.00327524740575,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.033982707798170,
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
            "CCSD CORRELATION ENERGY": -0.2068152041,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.0478712079,
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
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001767468898,
            "MP2 CORRELATION ENERGY": -0.058423513790,
            "MP2 SINGLES ENERGY": 0.0,
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
            "CCSD CORRELATION ENERGY": -0.08117105566,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.002231965267,
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
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001851937488,
            "MP2 SINGLES ENERGY": -0.000688368657,
            "MP2 CORRELATION ENERGY": -0.059407254257,
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
            "CCSD CORRELATION ENERGY": -0.08256719,
            "CCSD SINGLES ENERGY": -0.00117001688,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.00230304,
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
            "CCSD CORRELATION ENERGY": -0.250330548844,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.054051928864870,
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
            "CCSD CORRELATION ENERGY": -0.2271733460,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.049398348010672,
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
            "CCSD CORRELATION ENERGY": -0.188317222733,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.036526852874970,
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
            "CCSD CORRELATION ENERGY": -0.1716495276680232,
            "CCSD SINGLES ENERGY": 0.0,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.033248190929062,
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
            "CCSD CORRELATION ENERGY": -0.19282621471297376,
            "CCSD SINGLES ENERGY": -0.003354603508621,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.036502859698546,
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
            "CCSD CORRELATION ENERGY": -0.175988485854028,
            "CCSD SINGLES ENERGY": -0.003256808469230,
            "CCSD SAME-SPIN CORRELATION ENERGY": -0.033291143258924,
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
            "MP2 SAME-SPIN CORRELATION ENERGY": 0.13,
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
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -2.9,
            "MP2 CORRELATION ENERGY": -0.20162439774,
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
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0017676971,
            "MP2 CORRELATION ENERGY": -0.05841222894,
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
            "MP2 SINGLES ENERGY": -0.00068836865,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.001852162877,
            "MP2 CORRELATION ENERGY": -0.059395907176,
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
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 CORRELATION ENERGY": -0.20377997248921056,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05431321036920538,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.0027998, 0.0, 0.0, -0.0027998]  # dfmp2 findif-5 ae pk+df
            ).reshape((-1, 3)),
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
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 CORRELATION ENERGY": -0.05945820694747983,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0019203155958724552,
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
            "MP2 SINGLES ENERGY": -0.0006940498589629459,
            "MP2 CORRELATION ENERGY": -0.0604460449537298,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0020066877639503184,
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
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.05352569481658172,
            "MP2 CORRELATION ENERGY": -0.20162566806258586,
            "MP2 TOTAL GRADIENT": np.array(
                [0.0, 0.0, 0.00315485, 0.0, 0.0, -0.00315485]  # dfmp2 findif-5 fc pk+df
            ).reshape((-1, 3)),
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
            "MP2 SINGLES ENERGY": 0.0,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0017690135626491292,
            "MP2 CORRELATION ENERGY": -0.058392397606538686,
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
            "MP2 SINGLES ENERGY": -0.0006883686516107368,
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.0018536363586657242,
            "MP2 CORRELATION ENERGY": -0.05937514348825628,
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
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.87,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 ae cd+cd
                [0.0, 0.0, 0.00281146, 0.0, 0.0, -0.00281146]
            ).reshape((-1, 3)),
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
            "MP2 SAME-SPIN CORRELATION ENERGY": -0.46,
            "MP2 TOTAL GRADIENT": np.array(
                # dfocc findif-5 fc cd+cd
                [0.0, 0.0, 0.00316665, 0.0, 0.0, -0.00316665]
            ).reshape((-1, 3)),
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
    #            "MP2 SINGLES ENERGY": 0.0,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.3,
    #            "MP2 CORRELATION ENERGY": -2.3,
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
    #            "MP2 SINGLES ENERGY": 0.0,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.4,
    #            "MP2 CORRELATION ENERGY": -2.4,
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
    #            "MP2 SINGLES ENERGY": 0.0,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.5,
    #            "MP2 CORRELATION ENERGY": -2.5,
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
    #            "MP2 SINGLES ENERGY": -2.7,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.7,
    #            "MP2 CORRELATION ENERGY": -2.7,
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
    #            "MP2 SINGLES ENERGY": -2.8,
    #            "MP2 SAME-SPIN CORRELATION ENERGY": -2.8,
    #            "MP2 CORRELATION ENERGY": -2.8,
    #        },
    #    },
]


for calc in _std_suite:
    if calc["data"]:
        calc["data"]["MP2 TOTAL ENERGY"] = calc["data"]["MP2 CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
        calc["data"]["MP2 DOUBLES ENERGY"] = calc["data"]["MP2 CORRELATION ENERGY"] - calc["data"]["MP2 SINGLES ENERGY"]
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

        if "CCSD CORRELATION ENERGY" in calc["data"]:
            calc["data"]["CCSD TOTAL ENERGY"] = (
                calc["data"]["CCSD CORRELATION ENERGY"] + calc["data"]["HF TOTAL ENERGY"]
            )
            calc["data"]["CCSD DOUBLES ENERGY"] = (
                calc["data"]["CCSD CORRELATION ENERGY"] - calc["data"]["CCSD SINGLES ENERGY"]
            )
            calc["data"]["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = (
                calc["data"]["CCSD CORRELATION ENERGY"]
                - calc["data"]["CCSD SAME-SPIN CORRELATION ENERGY"]
                - calc["data"]["CCSD SINGLES ENERGY"]
            )


def answer_hash(**kwargs):
    system = kwargs.pop("system")
    basis = kwargs.pop("basis")
    scf_type = kwargs.pop("scf_type")
    reference = kwargs.pop("reference")
    fcae = kwargs.pop("fcae")
    corl_type = kwargs.pop("corl_type")

    return "_".join([system, basis, scf_type, reference, fcae, corl_type])


std_suite = {answer_hash(**calc["meta"]): calc["data"] for calc in _std_suite}
