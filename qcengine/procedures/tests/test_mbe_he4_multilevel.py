import copy
import pprint

import pytest

from qcelemental import constants
from qcelemental.models import Molecule
from qcelemental.models.procedures_manybody import AtomicSpecification, ManyBodyKeywords, ManyBodyInput
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using

def skprop(qcvar):
    return qcng.procedures.manybody.qcvars_to_manybodyproperties[qcvar]


@pytest.fixture(scope="function")
def mbe_data_multilevel_631g():
    # note that spherical/cartesian irrelevant for He & 6-31G, and fc/ae irrelevant for He
    c4_kwds = {}
    gms_kwds = {"basis__ngauss": 6, "ccinp__ncore": 0, "ccinp__iconv": 9, "scf__conv": 9}
    nwc_kwds = {"scf__thresh": 1.0e-8, "ccsd__thresh": 1.e-8}
    p4_kwds = {"scf_type": "pk", "mp2_type": "conv"}

    protocols = {"stdout": False}
    return {
        "specification": {
            "specification": {
                "c4-hf": {
                    "model": {
                        "method": "hf",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "cfour",
                    "keywords": c4_kwds,
                    "protocols": protocols,
                },
                "c4-mp2": {
                    "model": {
                        "method": "mp2",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "cfour",
                    "keywords": c4_kwds,
                    "protocols": protocols,
                },
                "c4-ccsd": {
                    "model": {
                        "method": "ccsd",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "cfour",
                    "keywords": c4_kwds,
                    "protocols": protocols,
                },
                "gms-hf": {
                    "model": {
                        "method": "hf",
                        "basis": "n31",
                    },
                    "driver": "energy",
                    "program": "gamess",
                    "keywords": gms_kwds,
                    "protocols": protocols,
                },
                "gms-mp2": {
                    "model": {
                        "method": "mp2",
                        "basis": "n31",
                    },
                    "driver": "energy",
                    "program": "gamess",
                    "keywords": gms_kwds,
                    "protocols": protocols,
                },
                "gms-ccsd": {
                    "model": {
                        "method": "ccsd",
                        "basis": "n31",
                    },
                    "driver": "energy",
                    "program": "gamess",
                    "keywords": gms_kwds,
                    "protocols": protocols,
                },
                "nwc-hf": {
                    "model": {
                        "method": "hf",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "nwchem",
                    "keywords": nwc_kwds,
                    "protocols": protocols,
                },
                "nwc-mp2": {
                    "model": {
                        "method": "mp2",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "nwchem",
                    "keywords": nwc_kwds,
                    "protocols": protocols,
                },
                "nwc-ccsd": {
                    "model": {
                        "method": "ccsd",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "nwchem",
                    "keywords": nwc_kwds,
                    "protocols": protocols,
                },
                "p4-hf": {
                    "model": {
                        "method": "hf",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "psi4",
                    "keywords": p4_kwds,
                    "protocols": protocols,
                },
                "p4-mp2": {
                    "model": {
                        "method": "mp2",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "psi4",
                    "keywords": p4_kwds,
                    "protocols": protocols,
                },
                "p4-ccsd": {
                    "model": {
                        "method": "ccsd",
                        "basis": "6-31g",
                    },
                    "driver": "energy",
                    "program": "psi4",
                    "keywords": p4_kwds,
                    "protocols": protocols,
                },
            },
            "keywords": None,
            "driver": "energy",
        },
        "molecule": None,
    }


he4_refs_conv_multilevel_631g = {
    # 1: ccsd; 2,3: mp2; 4: hf, all 6-31G
    "121": {
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":       -11.480648555603,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":       -11.472000052247,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":       -11.472089645469,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":       -11.472068853166,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":   0.008648503357,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":   0.008558910134,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":   0.008579702437,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":       0.008648503357,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":      -0.000089593222,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":       0.000020792303,

        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.480648555603,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.471058574581,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         -11.471324608815,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.471272244751,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":   0.009589981022,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":   0.009323946788,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":   0.009376310852,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":       0.009589981022,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":       -0.000266034234,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":       0.000052364064,
    },
    # 1,2: ccsd; 3,4: mp2, all 6-31G
    "22": {
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":    -11.480648555603,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":    -11.471764016410,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":    -11.471853609632,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":    -11.471834096023,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":    0.008884539193,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":    0.008794945971,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":    0.008814459580,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":    0.008884539193,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":    -0.000089593222,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":    0.000019513609,

        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":      -11.480648555603,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":      -11.470705938773,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":      -11.470971973006,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":      -11.470913449084,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":      0.009942616831,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":      0.009676582597,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":      0.009735106519,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":      0.009942616831,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":      -0.000266034234,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":      0.000058523922,
    },
}

# only here for keys
he4_refs_conv = {
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":           -11.530668717083888,
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":           -11.522467757090013,
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":           -11.522702864080149,
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":           -11.522639870651439,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":       0.008200959993875045,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":       0.007965853003739198,
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":       0.008028846432448944,
        "CP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":           0.008200959993875045,
        "CP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":          -0.00023510699013584713,
        "CP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":           6.299342870974556e-05,

        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.530668717083888,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.522851206178828,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         -11.523095269671348,
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.523038093664368,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     0.007817510905059777,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":     0.0075734474125397355,
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     0.007630623419519367,
        "NOCP-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         0.007817510905059777,
        "NOCP-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        -0.00024406349252004134,
        "NOCP-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":         5.717600697963121e-05,

        "VMFC-CORRECTED TOTAL ENERGY THROUGH 1-BODY":         -11.530668717083888,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 2-BODY":         -11.52244892169719,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 3-BODY":         -11.52268452228489,
        "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY":         -11.522621528856181,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 2-BODY":     0.00821979538669737,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 3-BODY":     0.007984194798996924,
        "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY":     0.00804718822770667,
        "VMFC-CORRECTED 2-BODY CONTRIBUTION TO ENERGY":         0.00821979538669737,
        "VMFC-CORRECTED 3-BODY CONTRIBUTION TO ENERGY":        -0.00023560058770044634,
        "VMFC-CORRECTED 4-BODY CONTRIBUTION TO ENERGY":         6.299342870974556e-05,
}

sumdict = {
    "4b_all": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED TOTAL ENERGY": "VMFC-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "VMFC-CORRECTED INTERACTION ENERGY": "VMFC-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocpcp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocp_rtd_sio": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocp_sio": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_cp_rtd_sio": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_cp_sio": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocp_rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_cp_rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "4b_cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
    },
    "3b_nocp_rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b_nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b_cp_rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "3b_cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
    },
    "2b_nocp_rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b_nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b_cp_rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "2b_cp": {
        "CP-CORRECTED INTERACTION ENERGY": "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
    },
    "1b_nocp_rtd": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "zero",  #"NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b_nocp": {
        "NOCP-CORRECTED TOTAL ENERGY": "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "NOCP-CORRECTED INTERACTION ENERGY": "zero",  #"NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b_cp_rtd": {
        "CP-CORRECTED TOTAL ENERGY": "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        "CP-CORRECTED INTERACTION ENERGY": "zero",  #"CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
    "1b_cp": {
        "CP-CORRECTED INTERACTION ENERGY": "zero",  #"CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
    },
# TODO table defines the general qcvar as 0 even if 1-body qcvar not available. continue?
}


@pytest.fixture
def he_tetramer():
    a2 = 2 / constants.bohr2angstroms
    return Molecule(symbols=["He", "He", "He", "He"], fragments=[[0], [1], [2], [3]], geometry=[0, 0, 0, 0, 0, a2, 0, a2, 0, 0, a2, a2])


@pytest.mark.parametrize("levels", [
    # pattern 121
    #pytest.param({4: "c4-hf", 3: "c4-mp2", 1: "c4-ccsd"}, id="121-cfour_pure", marks=using("cfour")),
    ##pytest.param({4: "gms-hf", 3: "gms-mp2", 1: "gms-ccsd"}, id="121-gamess_pure", marks=using("gamess")),
    #pytest.param({4: "nwc-hf", 3: "nwc-mp2", 1: "nwc-ccsd"}, id="121-nwchem_pure", marks=using("nwchem")),
    pytest.param({4: "p4-hf", 3: "p4-mp2", 1: "p4-ccsd"}, id="121-psi4_pure", marks=using("psi4")),

    #pytest.param({4: "p4-hf", 3: "c4-mp2", 1: "c4-ccsd"}, id="121-cfour_psi4", marks=[using("cfour"), using("psi4")]),
    #pytest.param({4: "nwc-hf", 3: "nwc-mp2", 1: "p4-ccsd"}, id="121-nwchem_psi4", marks=[using("nwchem"), using("psi4")]),
    #pytest.param({4: "c4-hf", 3: "nwc-mp2", 1: "p4-ccsd"}, id="121-cfour_nwchem_psi4", marks=[using("cfour"), using("nwchem"), using("psi4")]),

    # pattern 22
    pytest.param({4: "p4-mp2", 2: "p4-ccsd"}, id="22-psi4_pure", marks=using("psi4")),
])
@pytest.mark.parametrize("mbe_keywords,anskey,bodykeys,calcinfo_nmbe", [
#    pytest.param(
#        {"bsse_type": ["nocp", "cp", "vmfc"]},
#        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv],
#        {"121": 65,
#         "22": 99},  #
#        id="4b_all"),
    pytest.param(  # ODD entry
        {"bsse_type": ["nocp", "cp"]},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if not k.startswith("VMFC-")],
        {"121": 61,  # cp(14md + 15lo) + nocp(14md + 15lo) + 4hi - 1lo 1234@1234
         "22": 49},  # cp(10hi + 15lo) + nocp(10hi + 15lo) - 1lo 1234@1234
        id="4b_nocpcp"),
#    pytest.param(
#        {"bsse_type": "nocp", "return_total_data": True, "supersystem_ie_only": True},
#        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
#        {"121": 5,
#         "22": 99},  #
#        id="4b_nocp_rtd_sio"),
#    pytest.param(
#        {"bsse_type": "nocp", "return_total_data": False, "supersystem_ie_only": True},
#        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
#        {"121": 5,
#         "22": 99},  #
#        id="4b_nocp_sio"),
#    pytest.param(
#        {"bsse_type": "cp", "return_total_data": True, "supersystem_ie_only": True},
#        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv if (k.startswith("CP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k))],
#        {"121": 9,
#         "22": 99},  #
#        id="4b_cp_rtd_sio"),
#    pytest.param(
#        {"bsse_type": "cp", "return_total_data": False, "supersystem_ie_only": True},
#        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
#        [k for k in he4_refs_conv if (k.startswith("CP-") and ("THROUGH 4-BODY" in k or "THROUGH 1-BODY" in k) and "TOTAL ENERGY" not in k)],
#        {"121": 5,
#         "22": 99},  #
#        id="4b_cp_sio"),
### TODO add vmfc. 3b nmbe=50
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-"))],
        {"121": 33,  # 4hi + 14md + 15lo vs. 15 for single-level
         "22": 25},  # 10hi + 15lo
        id="4b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-"))],
        {"121": 33,  # could be 29 TODO,  # 14md + 15lo vs. 15 for single-level
         "22": 25},  # 10hi + 15lo
        id="4b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True},
        "CP-CORRECTED TOTAL ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-"))],
        {"121": 33,  # 4hi(nocp) + 14md + 15lo vs. 19 for single-level,
         "22": 29},  # 10hi + 15lo + 4hi(nocp)
        id="4b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 4-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "TOTAL ENERGY" not in k)],
        {"121": 29, # 14md + 15lo vs. 15 for single-level,
         "22": 25},  # 10hi + 15lo
        id="4b_cp"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 3},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k))],
        {"121": 18,  # 4hi + 14md vs. 14 for single-level
         "22": 24},  # 10hi + 14lo
        id="3b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 3},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and "4-BODY" not in k)],
        {"121": 18,  # 4hi + 14md vs. 14 for single-level
         "22": 24},  # 10hi + 14lo
        id="3b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 3},
        "CP-CORRECTED TOTAL ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "4-BODY" not in k)],
        {"121": 18,  # 4hi + 14md vs. 18 for single-level  # bugfix: was 28
         "22": 28},  # 10hi + 14lo + 4hi(nocp)
        id="3b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 3},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 3-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and "4-BODY" not in k and "TOTAL ENERGY" not in k)],
        {"121": 14,  # 14md vs. 14 for single-level
         "22": 24},  # 10hi + 14lo
        id="3b_cp"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 2},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        {"121": 14,  # 4hi + 10md vs. 10 for single-level
         "22": 10},  # 10hi
        id="2b_nocp_rtd"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 2},
        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        {"121": 14,  # 4hi + 10md vs. 10 for single-level,
         "22": 10},  # 10hi
        id="2b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 2},
        "CP-CORRECTED TOTAL ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("4-BODY" not in k) and ("3-BODY" not in k))],
        {"121": 14,  # 4hi + 10md vs. 14 for single-level,
         "22": 14},  # 10hi + 4hi(nocp)
        id="2b_cp_rtd"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 2},
        "CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("4-BODY" not in k) and ("3-BODY" not in k) and "TOTAL ENERGY" not in k)],
        {"121": 10,  # 10md vs. 10 for single-level,
         "22": 10},  # 10hi
        id="2b_cp"),
    pytest.param(
        {"bsse_type": "nocp", "return_total_data": True, "max_nbody": 1},
        "NOCP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("1-BODY" in k))],
        {"121": 4,  # 4hi
         "22": 4},
        id="1b_nocp_rtd"),
## TODO fix 1b for rtd=F
##    pytest.param(
##        {"bsse_type": "nocp", "return_total_data": False, "max_nbody": 1},
##        "NOCP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
##        [k for k in he4_refs_conv if (k.startswith("NOCP-") and ("1-BODY" in k))],
##        {"121": 10,
##         "22": 99},  #
##        id="1b_nocp"),
    pytest.param(
        {"bsse_type": "cp", "return_total_data": True, "max_nbody": 1},
        "CP-CORRECTED TOTAL ENERGY THROUGH 1-BODY",
        [k for k in he4_refs_conv if (k.startswith("CP-") and ("1-BODY" in k))],
        {"121": 4,  # 4hi
         "22": 4},
        id="1b_cp_rtd"),
##    pytest.param(
##        {"bsse_type": "cp", "return_total_data": False, "max_nbody": 1},
##        "CP-CORRECTED INTERACTION ENERGY THROUGH 1-BODY",
##        [k for k in he4_refs_conv if (k.startswith("CP-") and ("1-BODY" in k) and "TOTAL ENERGY" not in k)],
##        {"121": 4,
##         "22": 99},  #
##        id="1b_cp"),
])
def test_nbody_he4_multi(levels, mbe_keywords, anskey, bodykeys, calcinfo_nmbe, he_tetramer, request, mbe_data_multilevel_631g):
    _inner = request.node.name.split("[")[1].split("]")[0]
    kwdsln, pattern, progln = _inner.split("-")

    levels = copy.deepcopy(levels)
    if pattern == "121":
        if mbe_keywords.get("max_nbody", 4) == 3:
            del levels[4]  # max_nbody and levels silently contradict w/o this
        elif mbe_keywords.get("max_nbody", 4) == 2:
            levels = {2: levels[3], 1: levels[1]}
        elif mbe_keywords.get("max_nbody", 4) == 1:
            del levels[4]
            del levels[3]
    elif pattern == "22":
        if mbe_keywords.get("max_nbody", 4) == 3:
            levels = {3: levels[4], 2: levels[2]}
        if mbe_keywords.get("max_nbody", 4) == 2:
            del levels[4]
        if mbe_keywords.get("max_nbody", 4) == 1:
            levels = {1: levels[2]}

    mbe_keywords = ManyBodyKeywords(levels=levels, **mbe_keywords)
    mbe_data_multilevel_631g["molecule"] = he_tetramer
    mbe_data_multilevel_631g["specification"]["keywords"] = mbe_keywords
    mbe_model = ManyBodyInput(**mbe_data_multilevel_631g)

    if "gamess" in progln:
        with pytest.raises(ValueError) as exe:
            qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
        assert "GAMESS+QCEngine can't handle ghost atoms yet" in str(exe.value)
        pytest.xfail("GAMESS can't do ghosts")

    ret = qcng.compute_procedure(mbe_model, "manybody", raise_error=True)
    print(f"MMMMMMM {request.node.name}")
    pprint.pprint(ret.model_dump(), width=200)

    refs = he4_refs_conv_multilevel_631g[pattern]
    ans = refs[anskey]
    ref_nmbe = calcinfo_nmbe[pattern]

    for qcv, ref in refs.items():
        skp = skprop(qcv)
        if qcv in bodykeys:
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[a] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[b] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[z] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv in sumdict["4b_all"]:
        skp = skprop(qcv)
        if qcv in sumdict[kwdsln]:
            refkey = sumdict[kwdsln][qcv]
            ref = 0.0 if refkey == "zero" else refs[refkey]
            assert compare_values(ref, ret.extras["qcvars"]["nbody"][qcv], atol=1.0e-8, label=f"[c] qcvars {qcv}")
            assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[d] skprop {skp}")
        else:
            assert qcv not in ret.extras["qcvars"]["nbody"], f"[y] {qcv=} wrongly present"
            assert getattr(ret.properties, skp) is None

    for qcv, ref in {
        "CURRENT ENERGY": ans,
    }.items():
        skp = skprop(qcv)
        assert compare_values(ref, ret.extras["qcvars"][qcv], atol=1.0e-8, label=f"[e] qcvars {qcv}")
        assert compare_values(ref, getattr(ret.properties, skp), atol=1.0e-8, label=f"[f] skprop {skp}")
    assert compare_values(ans, ret.return_result, atol=1.0e-8, label=f"[g] ret")
    assert ret.properties.calcinfo_nmbe == ref_nmbe, f"{ret.properties.calcinfo_nmbe=} != {ref_nmbe}"
