# This file sets up tests for fundamental QC methods for the purpose
#  of developing and maintaining harnesses to check that a QC program runs
#  through QCEngine, that it returns properties as expected for the method,
#  and that those results have the same answers as other programs. For
#  new programs, try to add the six rhf/uhf/rohf ae/fc combinations to
#  each method (separate "def test..."). Common things to account for are
#  sph/cart basis sets, fc/ae defaults, semicanonical orbitals.
#  Feel free to add test case variations to check keyword handling.

# Handling/Testing Errors
# -----------------------
# * (E1) Calculation that QC program just can't run.
#   example: GAMESS UHF CC
#   handle: add "error" tuple to test like `_q2`
#   result: runner executes through pytest.raises so test passes
#
# * (E2) Calculation that QC program can run but full properties data not harvestable.
#   example: Psi4 run with dfocc module doesn't split out MP2 spin components for RHF
#   handle: add logic in `contractual_mp2` (or appropriate fn) to forgive specific properties
#   result: runner checks that expected props are present and unexpected are absent.
#           even if full properties not available, QC code is working normally, so test passes
#
# * (E3) Testing instance parameters not evenly available to QC programs.
#   example: Odd basis like QZ2P not in Q-Chem's library
#   handle: add pytest.skip at most general position
#   result: test registers but skipped
#
# * (E4) Calculation the QC program thinks it can run but we know better -- `return_result` wrong.
#   example: NWChem ROHF MP2 that's missing singles contribution
#   handle: add "wrong" tuple to test like `_w1`
#   result: runner checks assertion fails at expected result and triggers pytest.xfail


import sys

import pytest
import qcelemental as qcel

import qcengine as qcng
from qcengine.programs.tests.standard_suite_cfour_basis import cfour_def2_svp_ecpdata, cfour_def2_svp_genbas
from qcengine.programs.tests.standard_suite_ref import std_molecules, std_refs
from qcengine.testing import using

from .standard_suite_runner import runner_asserter

_basis_keywords = ["cfour_basis", "basis"]


@pytest.fixture(params=["as_v2"])
def schema_versions_v2(request):
    # stdsuite runs QCSchema v2 in and out only. v1 and v1<->v2 conversion are covered by `schema_versions` tests
    #   elsewhere (e.g., test_canonical_*.py). the param id is needed, as the runner reads schema version from the test name.
    return qcel.models.v2, -1, qcel.models.v2


@pytest.fixture
def clsd_open_pmols():
    return {
        name[:-4]: qcel.models.v2.Molecule.from_data(smol, name=name[:-4])
        for name, smol in std_molecules.items()
        if name.endswith("-xyz")
    }


# <<<  Notes
#
# * `"xptd": {"fd": False}` is just a note that qcprog is doing internal findif, but the results are accurate
# * `"xptd": {"fd": True}` indicates internal findif, but probably 3-point, and results need a looser atol
# * `"xptd"` values may be per-dertype dicts, e.g., `"xptd": {"fd": {1: True}}`
# * `"marks"` is a per-dertype dict of marks applied at runtime, e.g., `"marks": {1: pytest.mark.xfail(reason=..., strict=True)}`.
#   selection marks like `long` must go on the pytest.param instead, since they're needed at collection.
#
# >>> Notes

# convergence settings for generating/matching references. qcdb applied e_convergence=10 & d_convergence=9 to all
#   calcs and translated per-program; here they're set natively per row.
_c4_tight = {
    "scf_conv": 12,
    "cc_conv": 12,
    "lineq_conv": 10,
    "bruck_conv": 7,  # 9 doesn't converge some dz
}
_p4_conv = {"e_convergence": 10, "d_convergence": 9}
_p4_fd = {"points": 5, "fd_project": False}  # tighten fd for a psi calc that's already doing it
_nwc_conv = {"scf__thresh": 1.0e-8}


# yapf: disable
# * tuple is (error_type, string_match_in_error_message, reason_for_human)
# * note that 2nd is regex matched, so raw strings and escape chars may be needed
# * numbering follows qcdb's stdsuite where tuples were ported from it (gaps are qcdb tuples for GAMESS or qcdb-only errors)
_q1 = (qcng.exceptions.InputError, "unknown SCFTYPE", "no ROHF reference for NWChem hand-coded MP2.")
_q2 = (qcng.exceptions.InputError, "CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF", "no UHF CC in GAMESS.")
_q3 = (qcng.exceptions.InputError, "ccsd: nopen is not zero", "no non-RHF reference for NWChem hand-coded CC.")
_q6 = (qcng.exceptions.InputError, r"Only RHF/UHF(/RKS|/RKS/UKS|) Hessians are currently implemented.", "no ROHF Hessian for Psi4 HF.")
_q7 = (qcng.exceptions.UnknownError, "CALCLEVEL LCCD not implemented for ROHF references.", "no ROHF LCCD in CFOUR.")
_q8 = (qcng.exceptions.InputError, r"Method=lccd is not available for any derivative level under conditions CC_TYPE=CONV, REFERENCE=ROHF, FREEZE_CORE=(TRUE|1|0|FALSE)", "no ROHF LCCD in Psi4")
_q10 = (qcng.exceptions.InputError, r"Invalid reference type (U|RO)HF != RHF for FNOCC energy", "no non-RHF LCCSD in Psi4.")
_q11 = (qcng.exceptions.InputError, r"Method=ccd is not available for any derivative level under conditions CC_TYPE=CONV, REFERENCE=(R|U|RO)HF, FREEZE_CORE=(TRUE|1|0|FALSE)", "no conv CCD in Psi4")
_q15 = (qcng.exceptions.InputError, r"Method=mp3 is not available for any derivative level under conditions MP_TYPE=CONV, REFERENCE=(U|RO)HF, FREEZE_CORE=(TRUE|1|0|FALSE)", "only detci for conv open-shell MP3 in Psi4, and it only runs rohf ZAPT")
_q17 = (qcng.exceptions.InputError, r"Method=ccd is not available for any derivative level under conditions CC_TYPE=CONV", "no conv CCD (so no CCD Hessian) in Psi4.")
_q18 = (qcng.exceptions.InputError, "Frozen core is not available for the CC gradients.", "no FC CC gradients in Psi4.")
_q19 = (qcng.exceptions.UnknownError, "UHF/ROHF gradients not implemented for this level of theory", "no non-RHF CCSDT-n or CCSDT gradients in CFOUR.")
_q20 = (qcng.exceptions.InputError, r"Method=ccsd\+t\(ccsd\) is not available for any derivative level", "no CCSD+T(CCSD) in Psi4.")
_q23 = (qcng.exceptions.InputError, r"Method not recognized: ccsdt-3", "no CCSDT-3 in NWChem.")
_q24 = (qcng.exceptions.InputError, r"Method not recognized: ccsdt-1a", "no CCSDT-1a in NWChem.")
_q25 = (qcng.exceptions.InputError, r"Method not recognized: ccsdt-1b", "no CCSDT-1b in NWChem.")
_q26 = (qcng.exceptions.InputError, r"(Method not recognized: ccsdt-2|Method=ccsdt-2 is not available for any derivative level)", "no CCSDT-2 in NWChem or Psi4.")
_q27 = (qcng.exceptions.InputError, r"Method not recognized: a-ccsd\(t\)", "no a-CCSD(T) in NWChem.")
_q28 = (qcng.exceptions.UnknownError, r"UHF has not been implemented for CCSD\(2\)_T", "no UHF a-CCSD(T) in CFOUR.")
_q29 = (qcng.exceptions.InputError, r"Method=a-ccsd\(t\) is not available for any derivative level under conditions CC_TYPE=CONV, REFERENCE=UHF, FREEZE_CORE=(TRUE|1|0|FALSE), QC_MODULE=\(auto\)", "no UHF a-CCSD(T) in Psi4 without MRCC")
_q31 = (qcng.exceptions.InputError, r"Method not recognized: mp4\(sdq\)", "no specialty MP4(SDQ) in NWChem.")
_q32 = (qcng.exceptions.InputError, r"Invalid reference type UHF != RHF for FNOCC energy", "no non-RHF MP4(SDQ) in Psi4.")
_q34 = (qcng.exceptions.InputError, r"Method=mp4 is not available for any derivative level under conditions MP_TYPE=CONV, REFERENCE=(U|RO)HF, FREEZE_CORE=(TRUE|1|0|FALSE)", "only detci for conv open-shell MP4 in Psi4, and it only runs rohf ZAPT")
_q35 = (qcng.exceptions.InputError, r"Method not recognized: (pbe|b3lyp|b3lyp5|b2plyp)", "no DFT in CFOUR.")
_q36 = (qcng.exceptions.InputError, "ROHF reference for DFT is not available.", "no ROHF DFT in Psi4.")
_q37 = (qcng.exceptions.InputError, r"Method not recognized: (pbe|b3lyp|b3lyp5)", "no DFT in CFOUR.")
_q38 = (qcng.exceptions.UnknownError, "RHF/UHF gradient calculations not possible for CALCLEVEL LCCD.", "no LCCD gradient in CFOUR.")
_q39 = (qcng.exceptions.InputError, r"Method=cisd is not available for any derivative level under conditions CI_TYPE=CONV, REFERENCE=UHF, FREEZE_CORE=(TRUE|1|0|FALSE)", "no UHF CI in Psi4")
_q41 = (qcng.exceptions.InputError, r"Invalid reference type UHF != RHF for FNOCC energy", "no non-RHF QCISD in Psi4.")
_q44 = (qcng.exceptions.InputError, r"Method not recognized: qcisd\(t\)", "no QCISD(T) in NWChem.")
_q45 = (qcng.exceptions.UnknownError, "non-Abelian symmetry not permitted", "temporary excuse of failure. I think NWChem has fixed upstream.")
_q46 = (qcng.exceptions.UnknownError, r"CALCLEVEL CC3 not implemented for ROHF reference", "no ROHF CC3 in CFOUR.")
_q47 = (qcng.exceptions.InputError, r"Method=(ccsdt-1a|ccsdt-1b|ccsdt-3|ccsdt) is not available for any derivative level", "no high-order CC in Psi4 without MRCC")
_q48 = (qcng.exceptions.InputError, r"Method=cc2 is not available for requested derivative level \(reqd=1 > avail=0\) under conditions CC_TYPE=CONV, REFERENCE=UHF, FREEZE_CORE=(TRUE|1|0|FALSE)", "no UHF CC2 gradients in Psi4")  # qcdb _q45

_w1 = ("MP2 CORRELATION ENERGY", "nonstandard answer: NWChem TCE MP2 doesn't report singles (affects ROHF)")
_w2 = ("CCSD CORRELATION ENERGY", "nonstandard answer: GAMESS CCSD ROHF FC energy")
_w3 = ("(T) CORRECTION ENERGY", "nonstandard answer: NWChem CCSD(T) ROHF AE/FC energy.")
_w4 = ("MP2 CORRELATION ENERGY", "nonstandard answer: NWChem TCE MP3 & MP4 doesn't report singles (affects ROHF), may be off by MP2 singles value.")
_w6 = ("CCSDTQ CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDTQ gradient mixed fc/ae parts.")
_w7 = ("CCSDT CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT & CCSDT(Q) gradient mixed fc/ae parts.")
_w8 = ("CCSD CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSD & A-CCSD(T) gradient mixed fc/ae parts.")
_w9 = ("CCD CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCD gradient mixed fc/ae parts.")
_w10 = ("MP2.5 CORRELATION ENERGY", "misdirected calc: CFOUR NCC MP3 gradient mixed fc/ae parts.")
_w11 = ("MP3 TOTAL GRADIENT", "nonstandard answer: CFOUR MP3 RHF FC doesn't match findif.")
_w12 = ("MP2 CORRELATION ENERGY", "nonstandard answer: CFOUR CCSD ROHF FC gradient right but energies wrong. (Paper says NYI.)")
_w13 = ("CCSDT-3 CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-3 gradient mixed fc/ae parts.")
_w14 = ("CCSDT-3 TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-3 FC doesn't match findif.")
_w15 = ("CCSDT-3 TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-3 FC doesn't match findif.")
_w16 = ("CCSDT-1A TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-1A FC doesn't match findif.")
_w17 = ("CCSDT-1A TOTAL HESSIAN", "nonstandard answer: CFOUR CCSDT-1A FC doesn't match findif.")
_w18 = ("CCSDT-1B TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-1B FC doesn't match findif.")
_w20 = ("CCSDT-2 TOTAL GRADIENT", "nonstandard answer: CFOUR CCSDT-2 FC doesn't match findif.")
_w22 = ("CCSDT-1A CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-1a gradient mixed fc/ae parts.")
_w23 = ("CCSDT-1B CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-1b gradient mixed fc/ae parts.")
_w24 = ("CCSDT-2 CORRELATION ENERGY", "misdirected calc: CFOUR NCC CCSDT-2 gradient mixed fc/ae parts.")
_w25 = ("A-CCSD(T) TOTAL GRADIENT", "nonstandard answer: CFOUR ECC A-CCSD(T) AE doesn't match findif.")
_w26 = ("MP2 TOTAL HESSIAN", "nonstandard answer: CFOUR MP2 ROHF hessian doesn't match findif. (Paper says NYI, but program does yield a value.)")
_w27 = ("MP2 CORRELATION ENERGY", "nonstandard answer: CFOUR MP2 ROHF hessian doesn't match findif, nor FC energy. (Paper says NYI, but program does yield a value.)")
_w28 = ("CISD CORRELATION ENERGY", "nonstandard answer: ROHF CISD vcc=tce!=guga=detci.")
_w31 = ("HF TOTAL HESSIAN", "nonstandard answer: CFOUR analytic Hessian omits ECP second-derivative integrals.")
_w32 = ("MP2 TOTAL HESSIAN", "nonstandard answer: CFOUR analytic Hessian omits ECP second-derivative integrals.")
# yapf: enable


# frozen-core orbital count where rows' hard-coded `dropmo: [1]` / `*__freeze: 1` (fine for 2nd-row systems) don't
#   apply. ECP systems: C 1s + I 4s4p (the 28-electron def2 ECP replaces 1s-3d), matching psi4's freeze_core
_fc_nfrzc = {
    ("def2-svp", "ch2i2"): 9,
    ("def2-svp", "ch2i"): 5,
}

def _system_keywords(qcprog, basis, method, driver):
    """Keywords a system needs that depend on more than one parametrize axis (program x basis x method x driver), so
    can't be written into a single pytest.param."""
    kw = {}
    if basis == "def2-svp" and qcprog == "psi4" and driver == "hessian" and method != "hf":
        # ECP systems ch2i2/ch2i: psi4 findif steps on I are mass-weighted to 0.005/sqrt(127) bohr, so findif-of-energies
        #   Hessians amplify energy noise. ch2i2 rhf fc mp2: e/d 10/9 3.9e-3 from cfour findif-of-analytic-grads; 12/11 5.8e-5
        kw.update({"e_convergence": 12, "d_convergence": 11})
    if basis == "def2-svp" and qcprog == "psi4" and driver == "gradient" and method in ["ccsd", "ccsd(t)"]:
        # ECP systems ch2i2/ch2i: psi4 has no fc (or rohf ccsd(t)) analytic CC gradients, so does findif-of-energies. at
        #   e/d 10/9 and default r_conv, ch2i rohf fc ccsd lands 1.2e-6 from refs. psi4 stdsuite tightens & uses Cartesian steps
        kw.update({"e_convergence": 12, "d_convergence": 11, "r_convergence": 11, "disp_mass_weight": False})
    return kw


# `AtomicInput.specification.extras` needed by a program for a basis (beyond keywords)
_basis_spec_extras = {
    ("cfour", "def2-svp"): {"cfour": {"genbas": cfour_def2_svp_genbas, "ecpdata": cfour_def2_svp_ecpdata}},
}


def _trans_key(qc, bas, key):
    # note that cfour isn't ready for keywords.basis instead of model.basis
    lkey = key.lower()

    # translate basis
    if lkey == "basis":
        return bas
    if bas.lower() == "cc-pvdz":
        if lkey == "cfour_basis":
            return "pvdz"
        elif lkey == "qcng_basis":
            return {
                "cfour": ("pvdz", {}),
                "gamess": ("ccd", {"contrl__ispher": 1}),
                "jaguar": (bas, {}),
                "nwchem": ("cc-pvdz", {"basis__spherical": True}),
                "psi4": (bas, {}),
                "qchem": ("cc-pvdz", {}),
            }[qc]
    elif bas.lower() == "aug-cc-pvdz":
        if lkey == "cfour_basis":
            return "aug-pvdz"
        elif lkey == "qcng_basis":
            return {
                "cfour": ("aug-pvdz", {}),
                "gamess": ("accd", {"contrl__ispher": 1}),
                "jaguar": ("cc-pVDZ++", {}),
                "nwchem": ("aug-cc-pvdz", {"basis__spherical": True}),
                "psi4": (bas, {}),
                "qchem": ("aug-cc-pvdz", {}),
            }[qc]
    elif bas.lower() == "def2-svp":
        # ECP systems ch2i2/ch2i: I carries the 28-electron def2 ECP. NWChem needs ECPs named per element, and I is the
        #   only ECP element in these systems. CFOUR's library lacks def2-SVP, so it comes through `_basis_spec_extras`
        if lkey == "qcng_basis":
            return {
                "cfour": ("def2-svp", {"ecp": "on"}),
                "gamess": (None, {}),  # untested
                "jaguar": (None, {}),  # untested
                "nwchem": ("def2-svp", {"basis__spherical": True, "ecp__I": "library def2-ecp"}),
                "psi4": (bas, {}),
                "qchem": (None, {}),  # untested
            }[qc]
    elif bas.lower() == "cfour-qz2p":
        if lkey == "cfour_basis":
            return "qz2p"
        elif lkey == "qcng_basis":
            return {
                "cfour": ("qz2p", {}),
                "gamess": (None, {}),  # not in GAMESS-US library
                "jaguar": (None, {}),  # not in Jaguar library
                "nwchem": (None, {}),  # not in NWChem library
                "psi4": (bas, {}),
                "qchem": (None, {}),  # not in Q-Chem library
            }[qc]

    sys.exit(1)


# http://patorjk.com/software/taag/#p=display&c=bash&f=Soft&t=MP3


#  ,--.  ,--.,------.    ,------.
#  |  '--'  ||  .---'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .--.  ||  `--,     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |  |  ||  |`       |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'  `--'`--'        `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                      `---' `---'
#  <<<  HF Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=pytest.mark.ecp),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {"scf_conv": 12},                                                                     }, id="hf  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                   }, id="hf  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "jaguar", "reference": "rhf",  "fcae": "ae", "keywords": {"dconv": 1.e-8, "econv": 1.e-9, "nops": 1},                                                                          }, id="hf  rhf ae: jaguar",     marks=using("jaguar")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv},                                                                        }, id="hf  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {"scf_type": "pk"},                                                                   }, id="hf  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "scf_conv": 12},                                                 }, id="hf  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "uhf",  "fcae": "ae", "keywords": {"contrl__scftyp": "uhf"},                                                            }, id="hf  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "jaguar", "reference": "uhf",  "fcae": "ae", "keywords": {"dconv": 1.e-8, "econv": 1.e-9, "iuhf": 1, "nops": 1},                                                               }, id="hf  uhf ae: jaguar",     marks=using("jaguar")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {"scf__uhf": True, "scf__thresh": 1.0e-8},                                            }, id="hf  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "scf_type": "pk"},                                               }, id="hf  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "scf_conv": 12},                                                }, id="hf rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rohf", "fcae": "ae", "keywords": {"contrl__scftyp": "rohf"},                                                           }, id="hf rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "jaguar", "reference": "rohf", "fcae": "ae", "keywords": {"dconv": 1.e-8, "econv": 1.e-9, "iuhf": 0, "nops": 1},                                                               }, id="hf rohf ae: jaguar",     marks=using("jaguar")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {"scf__rohf": True, "scf__thresh": 1.0e-8},                                           }, id="hf rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "scf_type": "pk"},                                              }, id="hf rohf ae: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_hf_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "hf")
    )


#
#  ,--.  ,--.,------.     ,----.                     ,--.,--.                 ,--.
#  |  '--'  ||  .---'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .--.  ||  `--,     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |  |  ||  |`       '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'  `--'`--'         `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  HF Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        #    pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=pytest.mark.ecp),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        ######## Are all possible ways of computing <method> working?

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {"scf_conv": 12},                                                                     }, id="hf  rhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                   }, id="hf  rhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "jaguar", "reference": "rhf",  "fcae": "ae", "keywords": {"dconv": 1.e-8, "econv": 1.e-9, "nops": 1},                                                                          }, id="hf  rhf ae: jaguar",     marks=using("jaguar")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {"scf__thresh": 1.e-8},                                                               }, id="hf  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {"scf_type": "pk"},                                                                   }, id="hf  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "scf_conv": 12},                                                 }, id="hf  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "uhf",  "fcae": "ae", "keywords": {"contrl__scftyp": "uhf"},                                                            }, id="hf  uhf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "jaguar", "reference": "uhf",  "fcae": "ae", "keywords": {"dconv": 1.e-8, "econv": 1.e-9, "iuhf": 1, "nops": 1},                                                               }, id="hf  uhf ae: jaguar",     marks=using("jaguar")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {"scf__uhf": True, "scf__thresh": 1.e-8},                                             }, id="hf  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "scf_type": "pk"},                                               }, id="hf  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "scf_conv": 12},                                                }, id="hf rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rohf", "fcae": "ae", "keywords": {"contrl__scftyp": "rohf"},                                                           }, id="hf rohf ae: gamess",     marks=using("gamess")),
        pytest.param({"call": "jaguar", "reference": "rohf", "fcae": "ae", "keywords": {"dconv": 1.e-8, "econv": 1.e-9, "iuhf": 0, "nops": 1},                                                               }, id="hf rohf ae: jaguar",     marks=using("jaguar")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {"scf__rohf": True, "scf__thresh": 1.e-8},                                            }, id="hf rohf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "scf_type": "pk"},                                              }, id="hf rohf ae: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_hf_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "hf")
    )


#
#  ,--.  ,--.,------.    ,--.  ,--.                     ,--.
#  |  '--'  ||  .---'    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  .--.  ||  `--,     |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  |  |  |  ||  |`       |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#  `--'  `--'`--'        `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  HF Hessian


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=pytest.mark.ecp),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {"scf_conv": 12},                                                                     "wrong": {"def2-svp": _w31}}, id="hf  rhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                              }, id="hf  rhf ae: gamess", marks=using("gamess")),
        pytest.param({"call": "jaguar", "reference": "rhf",  "fcae": "ae", "keywords": {"dconv": 1.e-8, "econv": 1.e-9, "nops": 1},                                                                                                     }, id="hf  rhf ae: jaguar", marks=using("jaguar")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {},                                                                                                              }, id="hf  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {"scf_type": "pk"},                                                                                              }, id="hf  rhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "scf_conv": 12},                                                 "wrong": {"def2-svp": _w31}}, id="hf  uhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "uhf",  "fcae": "ae", "keywords": {"contrl__scftyp": "uhf"},                                                                                       }, id="hf  uhf ae: gamess", marks=using("gamess")),
        pytest.param({"call": "jaguar", "reference": "uhf",  "fcae": "ae", "keywords": {"dconv": 1.e-8, "econv": 1.e-9, "iuhf": 1, "nops": 1},                                                                                          }, id="hf  uhf ae: jaguar", marks=using("jaguar")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {"scf__uhf": True},                                                                                              }, id="hf  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {"reference": "uhf", "scf_type": "pk"},                                                                          }, id="hf  uhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "scf_conv": 12},                                                "wrong": {"def2-svp": _w31}}, id="hf rohf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rohf", "fcae": "ae", "keywords": {"contrl__scftyp": "rohf", "scf__dirscf": True},                                                                 }, id="hf rohf ae: gamess", marks=using("gamess")),
        pytest.param({"call": "jaguar", "reference": "rohf", "fcae": "ae", "keywords": {"dconv": 1.e-8, "econv": 1.e-9, "iuhf": 0, "nops": 1},                                                                                          }, id="hf rohf ae: jaguar", marks=using("jaguar")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {"scf__rohf": True, "scf__thresh": 1.e-8},                                                                                             }, id="hf rohf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {"reference": "rohf", "scf_type": "pk"},                                                       "error": {2: _q6 }}, id="hf rohf ae: psi4",   marks=using("psi4")),
        # yapf: enable
    ],
)
def test_hf_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "hf")
    )


#
#
#  ,--.   ,--.,------.  ,---.     ,------.
#  |   `.'   ||  .--. ''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' | .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  <<<  MP2 Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=pytest.mark.ecp),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                   "keywords": {"scf_conv": 12},                                                                                                         }, id="mp2  rhf ae: cfour",            marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                   "keywords": {**_c4_tight},                                                                                                            }, id="mp2  rhf ae: cfour",            marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rhf",  "fcae": "ae",                                   "keywords": {"mp2__nacore": 0},                                                                                                       }, id="mp2  rhf ae: gamess",           marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv},                                                                                                            }, id="mp2  rhf ae: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",                                   "keywords": {**_nwc_conv, "qc_module": "directmp2"},                                                                                  }, id="mp2  rhf ae: nwchem-directmp2", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",                                   "keywords": {**_nwc_conv, "qc_module": "tce"},                                                              "error": {"cc-pvdz": _q45}}, id="mp2  rhf ae: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                   "keywords": {**_p4_conv, "qc_module": "fnocc", "mp2_type": "conv"},                                                                   }, id="mp2  rhf ae: psi4-fnocc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                   "keywords": {**_p4_conv, "qc_module": "occ", "mp2_type": "conv"},                                                                     }, id="mp2  rhf ae: psi4-occ",         marks=using("psi4")),
        pytest.param({"call": "qchem",  "reference": "rhf",  "fcae": "ae",                                   "keywords": {"N_FROZEN_CORE": 0},                                                                                                     }, id="mp2  rhf ae: qchem",            marks=using("qchem")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                   "keywords": {**_c4_tight, "dropmo": 1},                                                                                               }, id="mp2  rhf fc: cfour",            marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rhf",  "fcae": "fc",                                   "keywords": {},                                                                                                                       }, id="mp2  rhf fc: gamess",           marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "mp2__freeze": 1},                                                                                          }, id="mp2  rhf fc: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",                                   "keywords": {**_nwc_conv, "mp2__freeze": 1, "qc_module": "directmp2"},                                                                }, id="mp2  rhf fc: nwchem-directmp2", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",                                   "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                            "error": {"cc-pvdz": _q45}}, id="mp2  rhf fc: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",                                   "keywords": {**_nwc_conv, "mp2__freeze__core": 1},                                                                                    }, id="mp2  rhf fc: nwchem",           marks=using("nwchem")),
      # pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",                                   "keywords": {**_nwc_conv, "mp2__freeze__core__atomic": True},                                                                         }, id="mp2  rhf fc: nwchem",           marks=using("nwchem")),  # uses qcdb alias
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",                                   "keywords": {**_nwc_conv, "mp2__freeze__atomic": {"O": 1}},                                                                           }, id="mp2  rhf fc: nwchem",           marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                                   "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc", "mp2_type": "conv"},                                              }, id="mp2  rhf fc: psi4-fnocc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                                   "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "occ", "mp2_type": "conv"},                                                }, id="mp2  rhf fc: psi4-occ",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                                   "keywords": {**_c4_tight, "reference": "uhf"},                                                                                        }, id="mp2  uhf ae: cfour",            marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "uhf",  "fcae": "ae",                                   "keywords": {"contrl__scftyp": "uhf", "mp2__nacore": 0},                                                                              }, id="mp2  uhf ae: gamess",           marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "scf__uhf": True},                                                                                          }, id="mp2  uhf ae: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",                                   "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                      }, id="mp2  uhf ae: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                                   "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "occ", "mp2_type": "conv"},                                                 }, id="mp2  uhf ae: psi4-occ",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                                   "keywords": {**_c4_tight, "reference": "uhf", "dropmo": [1]},                                                                         }, id="mp2  uhf fc: cfour",            marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "uhf",  "fcae": "fc",                                   "keywords": {"contrl__scftyp": "uhf"},                                                                                                }, id="mp2  uhf fc: gamess",           marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "scf__uhf": True, "mp2__freeze": 1},                                                                        }, id="mp2  uhf fc: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",                                   "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                    }, id="mp2  uhf fc: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",                                   "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True, "qc_module": "occ", "mp2_type": "conv"},                            }, id="mp2  uhf fc: psi4-occ",         marks=using("psi4")),
        pytest.param({"call": "qchem",  "reference": "uhf",  "fcae": "fc",                                   "keywords": {"N_frozen_CORE": "fC"},                                                                                                  }, id="mp2  uhf fc: qchem",            marks=using("qchem")),

            # "cfour_occupation": [[3, 1, 1, 0], [3, 0, 1, 0]]
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                                   "keywords": {**_c4_tight, "reference": "rohf"},                                                                                       }, id="mp2 rohf ae: cfour",            marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rohf", "fcae": "ae",                                   "keywords": {"contrl__scftyp": "rohf", "mp2__nacore": 0, "mp2__ospt": "RMP"},                                                         }, id="mp2 rohf ae: gamess",           marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "scf__rohf": True},                                                                        "error": {0: _q1}}, id="mp2 rohf ae: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",                                   "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 0, "qc_module": "tce"},                                  "wrong": {0: _w1}}, id="mp2 rohf ae: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",                                   "keywords": {**_p4_conv, "reference": "rohf", "mp2_type": "conv", "qc_module": "occ"},                                                }, id="mp2 rohf ae: psi4-occ",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc",                                   "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1},                                                                          }, id="mp2 rohf fc: cfour",            marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rohf", "fcae": "fc",                                   "keywords": {"contrl__scftyp": "rohf", "mp2__ospt": "RMP"},                                                                           }, id="mp2 rohf fc: gamess",           marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "scf__rohf": True, "mp2__freeze": 1},                                                      "error": {0: _q1}}, id="mp2 rohf fc: nwchem-mp2grad",   marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc",                                   "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                  "wrong": {0: _w1}}, id="mp2 rohf fc: nwchem-tce",       marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc",                                   "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "qc_module": "occ", "mp2_type": "conv"},                           }, id="mp2 rohf fc: psi4-occ",         marks=using("psi4")),
        # yapf: enable
    ],
)
def test_mp2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "mp2")
    )


#
#  ,--.   ,--.,------.  ,---.      ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' | .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     '-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  MP2 Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=[pytest.mark.ecp, pytest.mark.long]),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                   "keywords": {**_c4_tight},                                                                                                                      }, id="mp2  rhf ae: cfour",          marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv},                                                                                                                      }, id="mp2  rhf ae: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},             "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                        "error": {"cc-pvdz": _q45}}, id="mp2  rhf ae: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                   "keywords": {**_p4_conv, "qc_module": "occ", "mp2_type": "conv"},                                                                               }, id="mp2  rhf ae: psi4",           marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                   "keywords": {**_c4_tight, "dropmo": 1},                                                                                                         }, id="mp2  rhf fc: cfour",          marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "mp2__freeze": 1},                                                                                                    }, id="mp2  rhf fc: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},             "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                                      "error": {"cc-pvdz": _q45}}, id="mp2  rhf fc: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False},            "keywords": {**_p4_conv, **_p4_fd, "freeze_core": True, "mp2_type": "conv"},                                                                    }, id="mp2  rhf fc: psi4",           marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                                   "keywords": {**_c4_tight, "reference": "uhf"},                                                                                                  }, id="mp2  uhf ae: cfour",          marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "scf__uhf": True},                                                                                                    }, id="mp2  uhf ae: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},             "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                                }, id="mp2  uhf ae: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                                   "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "occ", "mp2_type": "conv"},                                                           }, id="mp2  uhf ae: psi4",           marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                                   "keywords": {**_c4_tight, "reference": "uhf", "dropmo": [1]},                                                                                   }, id="mp2  uhf fc: cfour",          marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "scf__uhf": True, "mp2__freeze": 1},                                                                                  }, id="mp2  uhf fc: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},             "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                              }, id="mp2  uhf fc: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False},            "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "freeze_core": True, "mp2_type": "conv"},                                                }, id="mp2  uhf fc: psi4",           marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                                   "keywords": {**_c4_tight, "reference": "rohf"},                                                                                                 }, id="mp2 rohf ae: cfour",          marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "scf__rohf": True},                                                                                  "error": {1: _q1}}, id="mp2 rohf ae: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",                                   "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                              "wrong": {1: _w1}}, id="mp2 rohf ae: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "xptd": {"fd": False},            "keywords": {**_p4_conv, **_p4_fd, "reference": "rohf", "qc_module": "occ", "mp2_type": "conv"},                                                }, id="mp2 rohf ae: psi4",           marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc",                                   "keywords": {**_c4_tight, "reference": "rohf", "dropmo": [1]},                                                                                  }, id="mp2 rohf fc: cfour",          marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "xptd": {"qc_module": "mp2grad"}, "keywords": {**_nwc_conv, "scf__rohf": True, "mp2__freeze": 1},                                                                "error": {1: _q1}}, id="mp2 rohf fc: nwchem-mp2grad", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc",                                   "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                            "wrong": {1: _w1}}, id="mp2 rohf fc: nwchem-tce",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "xptd": {"fd": False},            "keywords": {**_p4_conv, **_p4_fd, "reference": "rohf", "freeze_core": True, "qc_module": "occ", "mp2_type": "conv"},                           }, id="mp2 rohf fc: psi4",           marks=[*using("psi4"), *using("psi4_112")]),
        # yapf: enable
    ],
)
def test_mp2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "mp2")
    )


#
#  ,--.   ,--.,------.  ,---.     ,--.  ,--.                     ,--.
#  |   `.'   ||  .--. ''.-.  \    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |'.'|  ||  '--' | .-' .'    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  |  |   |  ||  | --' /   '-.    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#  `--'   `--'`--'     '-----'    `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  MP2 Hessian


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=[pytest.mark.ecp, pytest.mark.long]),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight},                                                               "wrong": {"def2-svp": _w32}}, id="mp2  rhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv},                                                                                          }, id="mp2  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "mp2_type": "conv"},                                                             }, id="mp2  rhf ae: psi4",   marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1},                                                  "wrong": {"def2-svp": _w32}}, id="mp2  rhf fc: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "mp2__freeze": 1},                                                                        }, id="mp2  rhf fc: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "freeze_core": True, "mp2_type": "conv"},                                        }, id="mp2  rhf fc: psi4",   marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "uhf"},                                           "wrong": {"def2-svp": _w32}}, id="mp2  uhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True},                                                                        }, id="mp2  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "mp2_type": "conv"},                                         }, id="mp2  uhf ae: psi4",   marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "uhf", "dropmo": [1]},                            "wrong": {"def2-svp": _w32}}, id="mp2  uhf fc: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "mp2__freeze": 1},                                                      }, id="mp2  uhf fc: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "freeze_core": True, "mp2_type": "conv"},                    }, id="mp2  uhf fc: psi4",   marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "rohf"},                                 "wrong": {2: _w26, "def2-svp": _w32}}, id="mp2 rohf ae: cfour",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "rohf", "mp2_type": "conv"},                                        }, id="mp2 rohf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "rohf", "dropmo": [1]},                                    "wrong": {2: _w27}}, id="mp2 rohf fc: cfour",  marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "rohf", "freeze_core": True, "mp2_type": "conv"},                   }, id="mp2 rohf fc: psi4",   marks=using("psi4")),
        # for rohf, cfour analytic 2nd deriv (which cfour ppr doesn't admit to) doesn't match findif HbyE from psi4 not HbyG from p4c4
        # yapf: enable
    ],
)
def test_mp2_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "mp2")
    )


#
#  ,--.   ,--.,------. ,----.     ,------.
#  |   `.'   ||  .--. ''.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |  .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' /'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'     `----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  <<<  MP3 Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc doesn't compute mp3
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc"},                                                                                     }, id="mp3  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                     }, id="mp3  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, "qc_module": "tce"},                                                            "error": {"cc-pvdz": _q45}}, id="mp3  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "fnocc", "mp_type": "conv"},                                                                  }, id="mp3  rhf ae: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "occ", "mp_type": "conv"},                                                                    }, id="mp3  rhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc"},                                                                        }, id="mp3  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                                                        }, id="mp3  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                          "error": {"cc-pvdz": _q45}}, id="mp3  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc", "mp_type": "conv"},                                             }, id="mp3  rhf fc: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "occ", "mp_type": "conv"},                                               }, id="mp3  rhf fc: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc"},                                                                 }, id="mp3  uhf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                    }, id="mp3  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "occ", "mp_type": "conv"},                                                }, id="mp3  uhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc"},                                                    }, id="mp3  uhf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                  }, id="mp3  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True, "qc_module": "occ", "mp_type": "conv"},                           }, id="mp3  uhf fc: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf", "cc_program": "vcc"},                                                                }, id="mp3 rohf ae: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {"scf__rohf": True, "qc_module": "tce", "scf__thresh": 1.e-8},                                         "wrong": {0: _w4}}, id="mp3 rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, "reference": "rohf", "mp_type": "conv"},                                                 "error": {0: _q15}}, id="mp3 rohf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "cc_program": "vcc"},                                                   }, id="mp3 rohf fc: cfour",      marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                "wrong": {0: _w4}}, id="mp3 rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "mp_type": "conv"},                            "error": {0: _q15}}, id="mp3 rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_mp3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "mp3")
    )


#
#  ,--.   ,--.,------. ,----.      ,----.                     ,--.,--.                 ,--.
#  |   `.'   ||  .--. ''.-.  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |'.'|  ||  '--' |  .' <     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  |   |  ||  | --' /'-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'   `--'`--'     `----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  MP3 Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "vcc"},                                                                                                              }, id="mp3  rhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                                              }, id="mp3  rhf ae: cfour-ncc", marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                                     "error": {"cc-pvdz": _q45}}, id="mp3  rhf ae: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                        "keywords": {**_p4_conv, "qc_module": "occ", "mp_type": "conv"},                                                                                             }, id="mp3  rhf ae: psi4-occ",  marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc"},                                                                               "wrong": {1: _w11}}, id="mp3  rhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                                                               "wrong": {1: _w10}}, id="mp3  rhf fc: cfour-ncc", marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                                                   "error": {"cc-pvdz": _q45}}, id="mp3  rhf fc: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "freeze_core": True, "qc_module": "occ", "mp_type": "conv"},                                                              }, id="mp3  rhf fc: psi4-occ",  marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc"},                                                                                          }, id="mp3  uhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                                             }, id="mp3  uhf ae: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                        "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "occ", "mp_type": "conv"},                                                                         }, id="mp3  uhf ae: psi4-occ",  marks=using("psi4")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "mp_type": "conv", "function_kwargs": {"dertype": 0}},                                                }, id="mp3  uhf ae: psi4",      marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "uhf", "dropmo": [1]},                                                                                                }, id="mp3  uhf fc: cfour",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                                           }, id="mp3  uhf fc: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "freeze_core": True, "qc_module": "occ", "mp_type": "conv"},                                          }, id="mp3  uhf fc: psi4-occ",  marks=[*using("psi4"), *using("psi4_112")]),
        # DEBUG pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "freeze_core": True, "mp_type": "conv", "function_kwargs": {"dertype": 0}},                           }, id="mp3  uhf fc: psi4",      marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "rohf"},                                                                                                              }, id="mp3 rohf ae: cfour",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",                        "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                                           "wrong": {1: _w1}}, id="mp3 rohf ae: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",                        "keywords": {**_p4_conv, "reference": "rohf", "mp_type": "conv", "points": 5},                                                             "error": {1: _q15}}, id="mp3 rohf ae: psi4",      marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "rohf", "dropmo": [1]},                                                                                               }, id="mp3 rohf fc: cfour",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc",                        "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                                         "wrong": {1: _w1}}, id="mp3 rohf fc: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc",                        "keywords": {**_p4_conv, "reference": "rohf", "mp_type": "conv", "freeze_core": True},                                                     "error": {1: _q15}}, id="mp3 rohf fc: psi4",      marks=using("psi4")),
        # yapf: enable
    ],
)
def test_mp3_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "mp3")
    )


#
#  ,--.   ,--.,------.   ,---.  ,-. ,---.  ,------.   ,-----.   ,-.      ,------.
#  |   `.'   ||  .--. ' /    | / .''   .-' |  .-.  \ '  .-.  '  '. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |/  '  ||  | `.  `-. |  |  \  :|  | |  |   |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' '--|  ||  | .-'    ||  '--'  /'  '-'  '-. |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'        `--' \ '.`-----' `-------'  `-----'--'.' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                               `-'                             `-'                                    `---' `---'
#  <<<  MP4(SDQ) Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc doesn't compute mp4
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc"},                                                           }, id="mp4_sdq_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                           }, id="mp4_sdq_  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv},                                                              "error": {0: _q31}}, id="mp4_sdq_  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "mp_type": "conv", "qc_module": "fnocc"},                                        }, id="mp4_sdq_  rhf ae: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc"},                                              }, id="mp4_sdq_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                              }, id="mp4_sdq_  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc", "mp_type": "conv"},                   }, id="mp4_sdq_  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc"},                                       }, id="mp4_sdq_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, "reference": "uhf", "mp_type": "conv"},                        "error": {0: _q32}}, id="mp4_sdq_  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc"},                          }, id="mp4_sdq_  uhf fc: cfour-vcc",  marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf", "cc_program": "vcc"},                                      }, id="mp4_sdq_ rohf ae: cfour-vcc",  marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "cc_program": "vcc"},                         }, id="mp4_sdq_ rohf fc: cfour-vcc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_mp4_prsdq_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "mp4(sdq)")
    )


#
#  ,--.   ,--.,------.   ,---.    ,------.
#  |   `.'   ||  .--. ' /    |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |'.'|  ||  '--' |/  '  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |   |  ||  | --' '--|  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'   `--'`--'        `--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                               `---' `---'
#  <<<  MP4 Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc"},                                                                   }, id="mp4  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                   }, id="mp4  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, "qc_module": "tce"},                                          "error": {"cc-pvdz": _q45}}, id="mp4  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "mp_type": "conv", "qc_module": "fnocc"},                                                }, id="mp4  rhf ae: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc"},                                                      }, id="mp4  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                                      }, id="mp4  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                        "error": {"cc-pvdz": _q45}}, id="mp4  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "mp_type": "conv", "qc_module": "fnocc"},                           }, id="mp4  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc"},                                               }, id="mp4  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                  }, id="mp4  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, "reference": "uhf", "mp_type": "conv"},                                "error": {0: _q34}}, id="mp4  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc"},                                  }, id="mp4  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                }, id="mp4  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True, "mp_type": "conv"},           "error": {0: _q34}}, id="mp4  uhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf", "cc_program": "vcc"},                                              }, id="mp4 rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                "wrong": {0: _w4}}, id="mp4 rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, "reference": "rohf", "mp_type": "conv"},                               "error": {0: _q34}}, id="mp4 rohf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "cc_program": "vcc"},                                 }, id="mp4 rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},              "wrong": {0: _w4}}, id="mp4 rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "mp_type": "conv"},          "error": {0: _q34}}, id="mp4 rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_mp4_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "mp4")
    )


#
#  ,-------.  ,---.  ,------. ,--------. ,---.     ,------.
#  `--.   /  /  O  \ |  .--. ''--.  .--''.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#    /   /  |  .-.  ||  '--' |   |  |    .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#   /   `--.|  | |  ||  | --'    |  |   /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-------'`--' `--'`--'        `--'   '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                `---' `---'
#  <<<  ZAPT2 Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        # pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),  # larger basis sets too long, and detci ae can't allocate arrays
        # pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, "reference": "rohf", "mp_type": "conv", "qc_module": "detci"},                     }, id="zapt2 rohf ae: psi4-detci", marks=using("psi4")),

        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "qc_module": "detci", "mp_type": "conv"},}, id="zapt2 rohf fc: psi4-detci", marks=using("psi4")),
        # yapf: enable
    ],
)
def test_zapt2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "zapt2")
    )


#
#   ,-----.,--. ,---.  ,------.      ,------.
#  '  .--./|  |'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |`.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\|  |.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'`--'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                  `---' `---'
#  <<<  CISD Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc, ncc don't run
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                                                 }, id="cisd  rhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                                    "error": {"cc-pvdz": _q45}}, id="cisd  rhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv, "qc_module": "fnocc"},                                                                                                             }, id="cisd  rhf ae: psi4-fnocc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv, "qc_module": "detci"},                                                                                                             }, id="cisd  rhf ae: psi4-detci",    marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                                                                  }, id="cisd  rhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                                                  "error": {"cc-pvdz": _q45}}, id="cisd  rhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",               "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc"},                                                                                        }, id="cisd  rhf fc: psi4-fnocc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",               "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "detci"},                                                                                        }, id="cisd  rhf fc: psi4-detci",    marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                                                             }, id="cisd  uhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                                            }, id="cisd  uhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",               "keywords": {**_p4_conv, "reference": "uhf"},                                                                                             "error": {0: _q39}}, id="cisd  uhf ae: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",               "keywords": {**_c4_tight, "refereNCE": "UHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                                              }, id="cisd  uhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                                          }, id="cisd  uhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",               "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True},                                                                        "error": {0: _q39}}, id="cisd  uhf fc: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "orbitals": 0, "print": 2},                                                             }, id="cisd rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2},                                                                            }, id="cisd rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                                                           }, id="cisd rohf ae sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",               "keywords": {**_p4_conv, "reference": "rohf", "qc_module": "detci"},                                                                      "wrong": {0: _w28}}, id="cisd rohf ae   : psi4-detci", marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                                             }, id="cisd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "orbitals": 0, "print": 2},                                              }, id="cisd rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                                                         }, id="cisd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        # pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc",               "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "qc_module": "detci", "e_convergence": 8, "r_convergence": 7},         "error": {2: _q8 }}, id="cisd rohf fc: psi4",          marks=using("psi4")),
        # yapf: enable
    ],
)
def test_cisd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "cisd")
    )


#
#   ,-----.    ,-----.,--. ,---.  ,------.      ,------.
#  '  .-.  '  '  .--./|  |'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  |  |  |    |  |`.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '-.'  '--'\|  |.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'--' `-----'`--'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                             `---' `---'
#  <<<  QCISD Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ncc errors out
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                              }, id="qcisd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                          }, id="qcisd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, "qc_module": "tce"},                                                 "error": {"cc-pvdz": _q45}}, id="qcisd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "fnocc"},                                                                          }, id="qcisd  rhf ae: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                               }, id="qcisd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ecc"},                                                           }, id="qcisd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                               "error": {"cc-pvdz": _q45}}, id="qcisd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc"},                                                     }, id="qcisd  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                          }, id="qcisd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "ecc"},                                                      }, id="qcisd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                         }, id="qcisd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, "reference": "uhf"},                                                          "error": {0: _q41}}, id="qcisd  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "refereNCE": "UHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                           }, id="qcisd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "refereNCE": "UHF", "dropmo": [1], "cc_program": "ecc"},                                       }, id="qcisd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                       }, id="qcisd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True},                                     "error": {0: _q41}}, id="qcisd  uhf fc: psi4",       marks=using("psi4")),
        # * no ROHF in Cfour, so hard to establish a ref
        # yapf: enable
    ],
)
def test_qcisd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "qcisd")
    )


#   ,-----.    ,-----.,--. ,---.  ,------.    ,-.,--------.,-.      ,------.
#  '  .-.  '  '  .--./|  |'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  | |  |  |  |    |  |`.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '-'  '-.'  '--'\|  |.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----'--' `-----'`--'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                             `-'          `-'                                    `---' `---'
#  <<<  QCISD(T) Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                      }, id="qcisd_t_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv},                                                                     "error": {0: _q44}}, id="qcisd_t_  rhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "fnocc"},                                                                  }, id="qcisd_t_  rhf ae: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                       }, id="qcisd_t_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc"},                                             }, id="qcisd_t_  rhf fc: psi4-fnocc", marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                  }, id="qcisd_t_  uhf ae: cfour-vcc",  marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "refereNCE": "UHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                   }, id="qcisd_t_  uhf fc: cfour-vcc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_qcisd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "qcisd(t)")
    )


#
#  ,------. ,-----.,--.    ,------.
#  |  .---''  .--./|  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  `--, |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  |`   '  '--'\|  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'     `-----'`--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                        `---' `---'
#  <<<  FCI Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "detci"},                                          }, id="fci  rhf ae: psi4-detci", marks=using("psi4")),  # 6m

        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "detci"},                     }, id="fci  rhf fc: psi4-detci", marks=using("psi4")),

        # TOO LONG pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, "freeze_core": False, "reference": "rohf", "qc_module": "detci"},}, id="fci rohf ae: psi4-detci", marks=using("psi4")),

        # pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "reference": "rohf", "qc_module": "detci"},}, id="fci rohf fc: psi4-detci", marks=using("psi4")),
        # yapf: enable
    ],
)
def test_fci_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "fci")
    )


#
#  ,--.    ,-----. ,-----.,------.      ,------.
#  |  |   '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |   |  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--.'  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-----' `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                     `---' `---'
#  <<<  LCCD Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc doesn't run
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                              }, id="lccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                          }, id="lccd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                 "error": {"cc-pvdz": _q45}}, id="lccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "fnocc"},                                                                                          }, id="lccd  rhf ae: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "occ"},                                                                                            }, id="lccd  rhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                                               }, id="lccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ncc"},                                                                           }, id="lccd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                               "error": {"cc-pvdz": _q45}}, id="lccd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc"},                                                                     }, id="lccd  rhf fc: psi4-fnocc", marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "occ"},                                                                       }, id="lccd  rhf fc: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                                          }, id="lccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                         }, id="lccd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "occ"},                                                                        }, id="lccd  uhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "refereNCE": "UHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                           }, id="lccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                       }, id="lccd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True, "qc_module": "occ"},                                                   }, id="lccd  uhf fc: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2},                                       "error": {0: _q7 }}, id="lccd rohf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                                        }, id="lccd rohf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, "reference": "rohf"},                                                                         "error": {0: _q8 }}, id="lccd rohf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "orbitals": 0, "print": 2},         "error": {0: _q7 }}, id="lccd rohf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                                      }, id="lccd rohf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "e_convergence": 8, "r_convergence": 7},            "error": {0: _q8 }}, id="lccd rohf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_lccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "lccd")
    )


#
#  ,--.    ,-----. ,-----.,------.       ,----.                     ,--.,--.                 ,--.
#  |  |   '  .--./'  .--./|  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |   |  |    |  |    |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--.'  '--'\'  '--'\|  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `-----' `-----' `-----'`-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  LCCD Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc doesn't run
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                              "error": {1: _q38}}, id="lccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "ncc"},                                                          "error": {1: _q38}}, id="lccd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "qc_module": "tce"},                                                   "error": {"cc-pvdz": _q45}}, id="lccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "qc_module": "fnocc"},                                                                  }, id="lccd  rhf ae: psi4-fnocc", marks=[*using("psi4"), *using("psi4_112")]),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                        "keywords": {**_p4_conv, "qc_module": "occ"},                                                                              }, id="lccd  rhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                               "error": {1: _q38}}, id="lccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ncc"},                                           "error": {1: _q38}}, id="lccd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                 "error": {"cc-pvdz": _q45}}, id="lccd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "freeze_core": True, "qc_module": "fnocc"},                                             }, id="lccd  rhf fc: psi4-fnocc", marks=[*using("psi4"), *using("psi4_112")]),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "freeze_core": True, "qc_module": "occ"},                                               }, id="lccd  rhf fc: psi4-occ",   marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                          "error": {1: _q38}}, id="lccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                           }, id="lccd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                        "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "occ"},                                                          }, id="lccd  uhf ae: psi4-occ",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "refereNCE": "UHF", "dropmo": [1], "cc_program": "vcc", "print": 2},           "error": {1: _q38}}, id="lccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                         }, id="lccd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "freeze_core": True, "qc_module": "occ"},                           }, id="lccd  uhf fc: psi4-occ",   marks=[*using("psi4"), *using("psi4_112")]),
        # yapf: enable
    ],
)
def test_lccd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "lccd")
    )


#
#  ,--.    ,-----. ,-----. ,---.  ,------.      ,------.
#  |  |   '  .--./'  .--./'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |   |  |    |  |    `.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--.'  '--'\'  '--'\.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `-----' `-----' `-----'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                             `---' `---'
#  <<<  LCCSD Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                              }, id="lccsd  rhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                          }, id="lccsd  rhf ae: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                 "error": {"cc-pvdz": _q45}}, id="lccsd  rhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv, "qc_module": "fnocc"},                                                                                          }, id="lccsd  rhf ae: psi4-fnocc",    marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                                               }, id="lccsd  rhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ncc"},                                                                           }, id="lccsd  rhf fc: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "qc_module": "tce", "tce__freeze": 1},                                               "error": {"cc-pvdz": _q45}}, id="lccsd  rhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",               "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc"},                                                                     }, id="lccsd  rhf fc: psi4-fnocc",    marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                                          }, id="lccsd  uhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                         }, id="lccsd  uhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",               "keywords": {**_p4_conv, "reference": "uhf"},                                                                          "error": {0: _q10}}, id="lccsd  uhf ae: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",               "keywords": {**_c4_tight, "REFerence": "UHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                           }, id="lccsd  uhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                       }, id="lccsd  uhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",               "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True},                                                     "error": {0: _q10}}, id="lccsd  uhf fc: psi4",          marks=using("psi4")),

        # Sum 2021: uncertain rohf reference values: cfour == nwc-tce, but those != p4n. cfour says rohf nyi. Sum 2022: switching from p4n to vcc ref values so c4 and nwc show as right, not `"wrong": {0: _w29}`
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2},                                                         }, id="lccsd rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "REFerence": "roHF", "orbitals": 0, "cc_program": "vcc", "print": 2},                                          }, id="lccsd rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_nwc_conv, "qc_module": "tce", "scf__rohf": True},                                                                        }, id="lccsd rohf ae sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",               "keywords": {**_p4_conv, "reference": "rohf"},                                                                         "error": {0: _q10}}, id="lccsd rohf ae   : psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                          }, id="lccsd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "orbitals": 0, "cc_program": "vcc", "print": 2},                           }, id="lccsd rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                                      }, id="lccsd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc",               "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True},                                                    "error": {0: _q10}}, id="lccsd rohf fc   : psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_lccsd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "lccsd")
    )


#
#   ,-----. ,-----.,------.      ,------.
#  '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                              `---' `---'
#  <<<  CCD Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc doesn't run
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                              }, id="ccd  rhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                          }, id="ccd  rhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                          }, id="ccd  rhf ae: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                 "error": {"cc-pvdz": _q45}}, id="ccd  rhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv},                                                                                              "error": {0: _q11}}, id="ccd  rhf ae: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                                               }, id="ccd  rhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ecc"},                                                                           }, id="ccd  rhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ncc"},                                                                           }, id="ccd  rhf fc: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                               "error": {"cc-pvdz": _q45}}, id="ccd  rhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",               "keywords": {**_p4_conv, "freeze_core": True},                                                                         "error": {0: _q11}}, id="ccd  rhf fc: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                                          }, id="ccd  uhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "ecc"},                                                                      }, id="ccd  uhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                         }, id="ccd  uhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",               "keywords": {**_p4_conv, "reference": "uhf"},                                                                          "error": {0: _q11}}, id="ccd  uhf ae: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",               "keywords": {**_c4_tight, "REFerence": "uHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                           }, id="ccd  uhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",               "keywords": {**_c4_tight, "REFerence": "uHF", "dropmo": [1], "cc_program": "ecc"},                                                       }, id="ccd  uhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                       }, id="ccd  uhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",               "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True},                                                     "error": {0: _q11}}, id="ccd  uhf fc: psi4",          marks=using("psi4")),

        # rohf vcc = tce, but cfour paper disavows rohf, so I'm suspicious
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "orbitals": 0, "print": 2},                                          }, id="ccd rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2},                                                         }, id="ccd rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                                        }, id="ccd rohf ae sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",               "keywords": {**_p4_conv, "reference": "rohf"},                                                                         "error": {0: _q11}}, id="ccd rohf ae   : psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "orbitals": 0, "cc_program": "vcc", "print": 2},                           }, id="ccd rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                          }, id="ccd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                                      }, id="ccd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc",               "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True},                                                    "error": {0: _q11}}, id="ccd rohf fc   : psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccd")
    )


#
#   ,-----. ,-----.,------.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./|  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\|  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCD Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                              }, id="ccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                          }, id="ccd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                          }, id="ccd  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True}, "keywords": {**_nwc_conv, "qc_module": "tce"},                                                 "error": {"cc-pvdz": _q45}}, id="ccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                       "keywords": {**_p4_conv, "cc_type": "conv"},                                                           "error": {1: _q11}}, id="ccd  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                       "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc", "print": 2},                                                 }, id="ccd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                       "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                                             }, id="ccd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                       "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                            "wrong": {1: _w9}}, id="ccd  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True}, "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                               "error": {"cc-pvdz": _q45}}, id="ccd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                       "keywords": {**_p4_conv, "freeze_core": True, "cc_type": "conv"},                                      "error": {1: _q11}}, id="ccd  rhf fc: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc", "print": 2},                                          }, id="ccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                                      }, id="ccd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True}, "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                         }, id="ccd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                       "keywords": {**_p4_conv, "reference": "uhf", "cc_type": "conv"},                                       "error": {1: _q11}}, id="ccd  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                       "keywords": {**_c4_tight, "reference": "uhf", "dropmo": [1], "cc_program": "vcc", "print": 2},                           }, id="ccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                       "keywords": {**_c4_tight, "reference": "uhf", "dropmo": [1], "cc_program": "ecc"},                                       }, id="ccd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True}, "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                       }, id="ccd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",                       "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True, "cc_type": "conv"},                  "error": {1: _q11}}, id="ccd  uhf fc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccd")
    )


#
#   ,-----. ,-----.,------.      ,--.  ,--.                     ,--.
#  '  .--./'  .--./|  .-.  \     |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    |  |  \  :    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\|  '--'  /    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-------'     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  CCD Hessian


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ncc errors
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                              }, id="ccd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                          }, id="ccd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True}, "keywords": {**_nwc_conv, "qc_module": "tce"},                                                 "error": {"cc-pvdz": _q45}}, id="ccd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                       "keywords": {**_p4_conv, "cc_type": "conv"},                                                           "error": {2: _q17}}, id="ccd  rhf ae: psi4",       marks=using("psi4")),

        # FC: vcc errors for analytic hess
        # pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                       "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc", "print": 2},                                                 }, id="ccd  rhf fc: cfour-vcc",  marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc", "print": 2},                                          }, id="ccd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                       "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                                      }, id="ccd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True}, "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                         }, id="ccd  uhf ae: nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                       "keywords": {**_p4_conv, "reference": "uhf", "cc_type": "conv"},                                       "error": {2: _q17}}, id="ccd  uhf ae: psi4",       marks=using("psi4")),
        # FC: vcc errors for analytic hess
        # pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                       "keywords": {**_c4_tight, "reference": "uhf", "dropmo": [1], "cc_program": "vcc", "print": 2},                           }, id="ccd  uhf fc: cfour-vcc",  marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccd_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "ccd")
    )


#
#  ,-----.   ,-----. ,-----.,------.      ,------.
#  |  |) /_ '  .--./'  .--./|  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \|  |    |  |    |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' /'  '--'\'  '--'\|  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------'  `-----' `-----'`-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                       `---' `---'
#  <<<  BCCD Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc, ncc doesn't run
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                    }, id="bccd  rhf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv},                                                                      }, id="bccd  rhf ae: psi4",         marks=using("psi4")),

        # cfour fc Brueckner hangs!
        # pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                     }, id="bccd  rhf fc: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",               "keywords": {**_p4_conv, "freeze_core": True},                                                 }, id="bccd  rhf fc: psi4",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                }, id="bccd  uhf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",               "keywords": {**_p4_conv, "reference": "uhf"},                                                  }, id="bccd  uhf ae: psi4",         marks=using("psi4")),

        # pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",               "keywords": {**_c4_tight, "REFerence": "uHF", "dropmo": [1], "cc_program": "vcc", "print": 2}, }, id="bccd  uhf fc: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",               "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True},                             }, id="bccd  uhf fc: psi4",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2},               }, id="bccd rohf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",               "keywords": {**_p4_conv, "reference": "rohf"},                                                 }, id="bccd rohf ae: psi4",         marks=using("psi4")),  # NOTE: sc

        # pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "print": 2},}, id="bccd rohf fc sc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True},                            }, id="bccd rohf fc sc: psi4",      marks=using("psi4")),
        # yapf: enable
    ],
)
def test_bccd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "bccd")
    )


#
#   ,-----. ,-----. ,---.     ,------.
#  '  .--./'  .--./'.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |     .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\/   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----''-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                           `---' `---'
#  <<<  CC2 Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc, ncc doesn't run
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                               }, id="cc2  rhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "qc_module": "tce"},                                                  "error": {"cc-pvdz": _q45}}, id="cc2  rhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv},                                                                                                 }, id="cc2  rhf ae: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                                }, id="cc2  rhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                "error": {"cc-pvdz": _q45}}, id="cc2  rhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",               "keywords": {**_p4_conv, "freeze_core": True},                                                                            }, id="cc2  rhf fc: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                           }, id="cc2  uhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                          }, id="cc2  uhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",               "keywords": {**_p4_conv, "reference": "uhf"},                                                                             }, id="cc2  uhf ae: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",               "keywords": {**_c4_tight, "REFerence": "uHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                            }, id="cc2  uhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                        }, id="cc2  uhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",               "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True},                                                        }, id="cc2  uhf fc: psi4",          marks=using("psi4")),

        # rohf ae vcc != tce == psi4, and cfour paper disavows rohf, so I'm suspicious
        # pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2},                                          }, id="cc2 rohf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",               "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                         }, id="cc2 rohf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",               "keywords": {**_p4_conv, "reference": "rohf"},                                                                            }, id="cc2 rohf ae: psi4",          marks=using("psi4")),  # NOTE: sc

        # nwc slightly different
        # pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                           }, id="cc2 rohf fc sc: cfour-vcc",  marks=using("cfour")),
        # pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                       }, id="cc2 rohf fc   : nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True},                                                       }, id="cc2 rohf fc sc: psi4",       marks=using("psi4")),
        # yapf: enable
    ],
)
def test_cc2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "cc2")
    )


#
#   ,-----. ,-----. ,---.      ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |     .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\/   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----''-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CC2 Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                                                  }, id="cc2  rhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "cc_type": "conv", "qc_module": "ccenergy"},                                                                                        }, id="cc2  rhf ae: psi4-cc",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cc_program": "vcc", "dropmo": 1, "print": 2},                                                                                     }, id="cc2  rhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "cc_type": "conv", "qc_module": "ccenergy"},                                                 "error": {1: _q18}}, id="cc2  rhf fc: psi4-cc",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "freeze_core": True, "cc_type": "conv", "qc_module": "ccenergy", "function_kwargs": {"dertype": 0}, "points": 5},                   }, id="cc2  rhf fc: psi4-cc",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc", "print": 2},                                                                              }, id="cc2  uhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, "reference": "uhf", "cc_type": "conv", "function_kwargs": {"dertype": 1}},                                        "error": {1: _q48}}, id="cc2  uhf ae: psi4-cc",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, "reference": "uhf", "cc_type": "conv", "function_kwargs": {"dertype": 0}, "points": 5},                                             }, id="cc2  uhf ae: psi4-cc",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc", "print": 2},                                                                 }, id="cc2  uhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "keywords": {**_p4_conv, "reference": "uhf", "cc_type": "conv", "freeze_core": True, "function_kwargs": {"dertype": 1}},                   "error": {1: _q48}}, id="cc2  uhf fc: psi4",      marks=using("psi4")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "keywords": {**_p4_conv, "reference": "uhf", "cc_type": "conv", "freeze_core": True, "function_kwargs": {"dertype": 0}, "points": 5},                        }, id="cc2  uhf fc: psi4",      marks=using("psi4")),
        # yapf: enable
    ],
)
def test_cc2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "cc2")
    )


#
#   ,-----. ,-----. ,---.  ,------.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                      `---' `---'
#  <<<  CCSD Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=pytest.mark.ecp),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                                                                           }, id="ccsd  rhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                                                                       }, id="ccsd  rhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                                                                       }, id="ccsd  rhf ae: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rhf",  "fcae": "ae",                                            "keywords": {"ccinp__ncore": 0},                                                                                                                                                      }, id="ccsd  rhf ae: gamess",        marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                                                              "error": {"cc-pvdz": _q45}}, id="ccsd  rhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv},                                                                                                                                                            }, id="ccsd  rhf ae: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "qc_module": "ccenergy"},                                                                                                                                    }, id="ccsd  rhf ae: psi4-cc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "qc_module": "fnocc"},                                                                                                                                       }, id="ccsd  rhf ae: psi4-fnocc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "qc_module": "mrcc"},                                                                                                                                        }, id="ccsd  rhf ae: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                                                                                            }, id="ccsd  rhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ecc"},                                                                                                                        }, id="ccsd  rhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ncc"},                                                                                                                        }, id="ccsd  rhf fc: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rhf",  "fcae": "fc",                                            "keywords": {},                                                                                                                                                                       }, id="ccsd  rhf fc: gamess",        marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                                                                            "error": {"cc-pvdz": _q45}}, id="ccsd  rhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv, "ccsd__freeze": 1},                                                                                                                                         }, id="ccsd  rhf fc: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "ccenergy"},                                                                                                               }, id="ccsd  rhf fc: psi4-cc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc"},                                                                                                                  }, id="ccsd  rhf fc: psi4-fnocc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "mrcc"},                                                                                                                   }, id="ccsd  rhf fc: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        # "cfour_occupation": [[3, 1, 1, 0], [3, 0, 1, 0]]
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                                                                                       }, id="ccsd  uhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "ecc"},                                                                                                                   }, id="ccsd  uhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "uhf",  "fcae": "ae",                                            "keywords": {"contrl__scftyp": "uhf", "ccinp__ncore": 0},                                                                                                           "error": {0: _q2 }}, id="ccsd  uhf ae: gamess",        marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                                                                      }, id="ccsd  uhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv, "scf__uhf": True},                                                                                                                        "error": {0: _q3 }}, id="ccsd  uhf ae: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "ccenergy"},                                                                                                                }, id="ccsd  uhf ae: psi4-cc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "mrcc"},                                                                                                                    }, id="ccsd  uhf ae: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "REFerence": "UHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                                                                        }, id="ccsd  uhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "REFerence": "UHF", "dropmo": [1], "cc_program": "ecc"},                                                                                                    }, id="ccsd  uhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "uhf",  "fcae": "fc",                                            "keywords": {"contrl__scftyp": "uhf"},                                                                                                                              "error": {0: _q2 }}, id="ccsd  uhf fc: gamess",        marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                                                                    }, id="ccsd  uhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv, "scf__uhf": True, "ccsd__freeze": 1},                                                                                                     "error": {0: _q3 }}, id="ccsd  uhf fc: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "freeze_core": True, "reference": "uhf", "qc_module": "ccenergy"},                                                                                           }, id="ccsd  uhf fc: psi4-cc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "freeze_core": True, "reference": "uhf", "qc_module": "mrcc"},                                                                                               }, id="ccsd  uhf fc: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sc",                              "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2},                                                                                                      }, id="ccsd rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sc",                              "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "ecc"},                                                                                                                  }, id="ccsd rohf ae sc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sd",                              "keywords": {**_c4_tight, "REFerence": "roHF", "orbitals": 0, "cc_program": "vcc", "print": 2},                                                                                       }, id="ccsd rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sd",                              "keywords": {**_c4_tight, "REFerence": "roHF", "orbitals": 0, "cc_program": "ecc"},                                                                                                   }, id="ccsd rohf ae sd: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rohf", "fcae": "ae",                                            "keywords": {"contrl__scftyp": "rohf", "ccinp__ncore": 0, "ccinp__maxcc": 50},                                                                                                        }, id="ccsd rohf ae: gamess",        marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",                                            "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                                                                                     }, id="ccsd rohf ae   : nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv, "scf__rohf": True},                                                                                                                       "error": {0: _q3 }}, id="ccsd rohf ae   : nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "sdsc": "sd",                              "keywords": {**_p4_conv, "reference": "rohf", "qc_module": "ccenergy"},                                                                                                               }, id="ccsd rohf ae sd: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "sdsc": "sc",                              "keywords": {**_p4_conv, "reference": "rohf", "qc_module": "mrcc"},                                                                                                                   }, id="ccsd rohf ae sc: psi4-mrcc",  marks=[*using("psi4"), *using("mrcc")]),

        # at least two formulations: standard (sd) or semicanonical (sc) orbitals. most implementations choose one.
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sd",                              "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "orbitals": 0, "cc_program": "vcc", "print": 2},                                                                        }, id="ccsd rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sd",                              "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "orbitals": 0, "cc_program": "ecc"},                                                                                    }, id="ccsd rohf fc sd: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                                                                       }, id="ccsd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "ecc"},                                                                                                   }, id="ccsd rohf fc sc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "gamess", "reference": "rohf", "fcae": "fc",                                            "keywords": {"contrl__scftyp": "rohf", "ccinp__iconv": 9, "scf__conv": 9},                                                                                          "wrong": {0: _w2 }}, id="ccsd rohf fc: gamess",        marks=using("gamess")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "sdsc": "sd",                              "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                                                                                   }, id="ccsd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv, "scf__rohf": True, "ccsd__freeze": 1},                                                                                                    "error": {0: _q3 }}, id="ccsd rohf fc   : nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sd",                              "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "qc_module": "ccenergy", "e_convergence": 8, "r_convergence": 7},                                                  }, id="ccsd rohf fc sd: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "qc_module": "ccenergy", "e_convergence": 8, "r_convergence": 7, "semicanonical": True},                           }, id="ccsd rohf fc sc: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "qc_module": "mrcc", "e_convergence": 8, "r_convergence": 7},                                                      }, id="ccsd rohf fc sc: psi4-mrcc",  marks=[*using("psi4"), *using("mrcc")]),
        # yapf: enable
    ],
)
def test_ccsd_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsd")
    )


# @pytest.mark.parametrize("mode", ["driver", "sandwich"])
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable

        ## local
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "vcc"},      "keywords": {"dropmo": 1, "print": 2},                                             }, id="ccsd  rhf fc: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "vcc"},      "keywords": {"dropmo": 1, "reference": "uhf", "print": 2},                         }, id="ccsd  uhf fc: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "vcc"},      "keywords": {"dropmo": 1, "reference": "rohf", "orbitals": 0, "print": 2},         }, id="ccsd rohf fc: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "vcc"},      "keywords": {"print": 2},                                                          }, id="ccsd  rhf ae: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "vcc"},      "keywords": {"reference": "uhf", "print": 2},                                      }, id="ccsd  uhf ae: cfour dd1",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "xptd": {"qc_module": "vcc"},      "keywords": {"reference": "rohf", "print": 2},                                     }, id="ccsd rohf ae: cfour dd1",  marks=using("cfour")),


        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "cc"},       "keywords": {**_nwc_conv, "ccsd__freeze": 1},                                      }, id="ccsd  rhf fc: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {**_nwc_conv, "tce__freeze": 1, "scf__uhf": True, "qc_module": "tce"}, }, id="ccsd  uhf fc: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {**_nwc_conv, "tce__freeze": 1, "scf__rohf": True, "qc_module": "tce"},}, id="ccsd rohf fc: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "cc"},       "keywords": {**_nwc_conv},                                                         }, id="ccsd  rhf ae: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                   }, id="ccsd  uhf ae: nwchem dd1", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "tce"},      "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                  }, id="ccsd rohf ae: nwchem dd1", marks=using("nwchem")),

        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {**_p4_conv, "freeze_core": True},                                     }, id="ccsd  rhf fc: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {**_p4_conv, "freeze_core": True, "reference": "uhf"},                 }, id="ccsd  uhf fc: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {**_p4_conv, "freeze_core": True, "reference": "rohf"},                }, id="ccsd rohf fc: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {**_p4_conv},                                                          }, id="ccsd  rhf ae: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {**_p4_conv, "reference": "uhf"},                                      }, id="ccsd  uhf ae: psi4 dd1",   marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "sdsc": "sd", "xptd": {"qc_module": "ccenergy"}, "keywords": {**_p4_conv, "reference": "rohf"},                                     }, id="ccsd rohf ae: psi4 dd1",   marks=using("psi4")),

        ## translated



        # yapf: enable
    ],
)
def test_ccsd_energy_default(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsd")
    )


#
#   ,-----. ,-----. ,---.  ,------.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSD Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=[pytest.mark.ecp, pytest.mark.long]),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                                                                          }, id="ccsd  rhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                                                                      }, id="ccsd  rhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                                                                      }, id="ccsd  rhf ae: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",               "xptd": {"fd": True},        "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                                                             "error": {"cc-pvdz": _q45}}, id="ccsd  rhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",               "xptd": {"fd": True},        "keywords": {**_nwc_conv},                                                                                                                                                           }, id="ccsd  rhf ae: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "cc_type": "conv", "qc_module": "ccenergy"},                                                                                                                }, id="ccsd  rhf ae: psi4-cc",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "cc_program": "vcc", "dropmo": 1, "print": 2},                                                                                                             }, id="ccsd  rhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "cc_program": "ecc", "dropmo": 1},                                                                                                                         }, id="ccsd  rhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "cc_program": "ncc", "dropmo": 1},                                                                                                        "wrong": {1: _w8}}, id="ccsd  rhf fc: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",               "xptd": {"fd": True},        "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                                                                           "error": {"cc-pvdz": _q45}}, id="ccsd  rhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",               "xptd": {"fd": True},        "keywords": {**_nwc_conv, "ccsd__freeze": 1},                                                                                                                                        }, id="ccsd  rhf fc: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "freeze_core": True, "cc_type": "conv", "qc_module": "ccenergy"},                                                                         "error": {1: _q18}}, id="ccsd  rhf fc: psi4-cc",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc", "print": 2},                                                                                                      }, id="ccsd  uhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                                                                                                  }, id="ccsd  uhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",               "xptd": {"fd": True},        "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                                                                     }, id="ccsd  uhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv, "scf__uhf": True},                                                                                                                        "error": {1: _q3}}, id="ccsd  uhf ae: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "reference": "uhf", "cc_type": "conv"},                                                                                                                     }, id="ccsd  uhf ae: psi4-cc",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc", "print": 2},                                                                                         }, id="ccsd  uhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "ecc"},                                                                                                     }, id="ccsd  uhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",               "xptd": {"fd": True},        "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                                                                   }, id="ccsd  uhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_nwc_conv, "scf__uhf": True, "ccsd__freeze": 1},                                                                                                     "error": {1: _q3}}, id="ccsd  uhf fc: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "reference": "uhf", "cc_type": "conv", "freeze_core": True},                                                                              "error": {1: _q18}}, id="ccsd  uhf fc: psi4",          marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                                            "keywords": {**_c4_tight, "reference": "rohf", "cc_program": "vcc", "print": 2},                                                                                                     }, id="ccsd rohf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                                            "keywords": {**_c4_tight, "reference": "rohf", "cc_program": "ecc", "print": 2},                                                                                                     }, id="ccsd rohf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",               "xptd": {"fd": True},        "keywords": {**_nwc_conv, "qc_module": "tce", "scf__rohf": True},                                                                                                                    }, id="ccsd rohf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",                                            "keywords": {**_nwc_conv, "scf__rohf": True},                                                                                                                       "error": {1: _q3}}, id="ccsd rohf ae: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "sdsc": "sd",                              "keywords": {**_p4_conv, "reference": "rohf", "qc_module": "ccenergy"},                                                                                                              }, id="ccsd rohf ae: psi4-cc",       marks=using("psi4")),

        # * c. Summer 2021 vcc and ecc yield correct gradients w/pert_orb=0 but not at the same time with correct energies (orbitals=0). c. Summer 2022: previous refers to "sd" targets. "sc" targets work w/ or w/o pert_orb=0
        # * note that ecc not recc for rohf and no rohf gradients on list in cfour paper
        # * below (for both ecc and vcc) no longer "wrong" for sc
            # pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc",                                            "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "cc_program": "vcc", "print": 2, "pert_orb": 0},                                                       "wrong": {1: _w12}}, id="ccsd rohf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "cc_program": "vcc", "print": 2},                                                                                        }, id="ccsd rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "cc_program": "ecc", "print": 2},                                                                                        }, id="ccsd rohf fc sc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"fd": True},        "keywords": {**_nwc_conv, "tce__freeze": 1, "scf__rohf": True, "qc_module": "tce"},                                                                                                  }, id="ccsd rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv, "ccsd__freeze": 1, "scf__rohf": True},                                                                                                    "error": {1: _q3}}, id="ccsd rohf fc   : nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc",                                            "keywords": {**_p4_conv, "reference": "rohf", "cc_type": "conv", "freeze_core": True},                                                                             "error": {1: _q18}}, id="ccsd rohf fc   : psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sd", "xptd": {"fd": False},       "keywords": {**_p4_conv, "reference": "rohf", "cc_type": "conv", "freeze_core": True, **_p4_fd, "function_kwargs": {"dertype": 0}},                                                  }, id="ccsd rohf fc sd: psi4-cc",    marks=[*using("psi4"), *using("psi4_112")]),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc", "xptd": {"fd": False},       "keywords": {**_p4_conv, "reference": "rohf", "cc_type": "conv", "freeze_core": True, **_p4_fd, "function_kwargs": {"dertype": 0}, "semicanonical": True},                           }, id="ccsd rohf fc sc: psi4-cc",    marks=[*using("psi4"), *using("psi4_112")]),
        # yapf: enable
    ],
)
def test_ccsd_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccsd")
    )


# PERT_ORB
# specifies the type of perturbed orbitals used in energy derivative calculations. STANDARD means that the gradient formulation assumes that the perturbed orbitals are not those in which the (perturbed) Fock matrix is diagonal. CANONICAL means that the perturbed orbitals are assumed to be canonical. This keyword is set automatically to CANONICAL in derivative calculations with methods which include triple excitations (MBPT[4]/MP4, CCSD+T[CCSD], CCSD[T], QCISD[T] and all iterative schemes like CCSDT-n and CC3) apart from CCSDT. IJ_CANONICAL requests a canonical perturbed-orbital treatment only for the occupied-occupied block of the unperturbed density matrix in analytic derivative calculations.
# For testing purpose, it is possible to force the use standard perturbed orbitals even in case of iterative triple excitations via the option FORCE_STANDA.
# Note also that in case of unrelaxed derivatives standard orbitals must be used.
# Default: STANDARD for all methods without triples (except CCSDT), CANONICAL for all methods with triples in case of relaxed derivatives.


#
#   ,-----. ,-----. ,---.  ,------.      ,--.  ,--.                     ,--.
#  '  .--./'  .--./'   .-' |  .-.  \     |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    `.  `-. |  |  \  :    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\.-'    ||  '--'  /    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-----' `-------'     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  CCSD Hessian


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ncc errors
        # * nwchem-cc Error in Communication
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                                                  }, id="ccsd  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "ecc", "print": 2},                                                                                                  }, id="ccsd  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                                     "error": {"cc-pvdz": _q45}}, id="ccsd  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "cc_type": "conv"},                                                                                                       }, id="ccsd  rhf ae: psi4",       marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc", "print": 2},                                                                                     }, id="ccsd  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc", "print": 2},                                                                                     }, id="ccsd  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                                                   "error": {"cc-pvdz": _q45}}, id="ccsd  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "freeze_core": True, "cc_type": "conv", "function_kwargs": {"dertype": 0}},                                               }, id="ccsd  rhf fc: psi4",       marks=[*using("psi4"), *using("psi4_112")]),  # FD SWITCH

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc", "print": 2},                                                                              }, id="ccsd  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                                                                          }, id="ccsd  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                                             }, id="ccsd  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                        "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "cc_type": "conv"},                                                                                   }, id="ccsd  uhf ae: psi4",       marks=[*using("psi4"), *using("psi4_112")]),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc", "print": 2},                                                                 }, id="ccsd  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "ecc", "print": 2},                                                                 }, id="ccsd  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                                           }, id="ccsd  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "freeze_core": True, "cc_type": "conv", "function_kwargs": {"dertype": 0}},                           }, id="ccsd  uhf fc: psi4",       marks=[*using("psi4"), *using("psi4_112")]),
        # FD SWITCH
        # yapf: enable
    ],
)
def test_ccsd_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "ccsd")
    )


#                                        ,--.
#   ,-----. ,-----. ,---.  ,------.      |  |    ,--------. ,-. ,-----. ,-----. ,---.  ,------. ,-.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \ ,---    ---.'--.  .--'/ .''  .--./'  .--./'   .-' |  .-.  \'. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :'---    ---'   |  |  |  | |  |    |  |    `.  `-. |  |  \  :|  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /    |  |       |  |  |  | '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'     `--'       `--'   \ '. `-----' `-----'`-----' `-------'.' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                           `-'                                 `-'                                    `---' `---'
#  <<<  CCSD+T(CCSD) Energy aka CCSD[T] Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                              "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                              }, id="ccsd_t_ccsd_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                              "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                          }, id="ccsd_t_ccsd_  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                              "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                          }, id="ccsd_t_ccsd_  rhf ae: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",                              "keywords": {**_nwc_conv, "qc_module": "tce"},                                                 "error": {"cc-pvdz": _q45}}, id="ccsd_t_ccsd_  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv},                                                                                               }, id="ccsd_t_ccsd_  rhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                              "keywords": {**_p4_conv},                                                                              "error": {0: _q20}}, id="ccsd_t_ccsd_  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                              "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                               }, id="ccsd_t_ccsd_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                              "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ecc"},                                                           }, id="ccsd_t_ccsd_  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                              "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ncc"},                                                           }, id="ccsd_t_ccsd_  rhf fc: cfour-ncc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",                              "keywords": {**_nwc_conv, "qc_module": "tce", "tce__freeze": 1},                               "error": {"cc-pvdz": _q45}}, id="ccsd_t_ccsd_  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv, "ccsd__freeze": 1},                                                                            }, id="ccsd_t_ccsd_  rhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                              "keywords": {**_p4_conv, "freeze_core": True},                                                         "error": {0: _q20}}, id="ccsd_t_ccsd_  rhf fc: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                              "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                          }, id="ccsd_t_ccsd_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                              "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "ecc"},                                                      }, id="ccsd_t_ccsd_  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",                              "keywords": {**_nwc_conv, "qc_module": "tce", "scf__uhf": True},                                                         }, id="ccsd_t_ccsd_  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",                              "keywords": {**_nwc_conv, "scf__uhf": True},                                                            "error": {0: _q3}}, id="ccsd_t_ccsd_  uhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                              "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "ccenergy"},                                 "error": {0: _q20}}, id="ccsd_t_ccsd_  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                              "keywords": {**_c4_tight, "dropmo": [1], "REFerence": "UHF", "cc_program": "vcc", "print": 2},                           }, id="ccsd_t_ccsd_  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                              "keywords": {**_c4_tight, "dropmo": [1], "REFerence": "UHF", "cc_program": "ecc"},                                       }, id="ccsd_t_ccsd_  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",                              "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce", "scf__uhf": True},                                       }, id="ccsd_t_ccsd_  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",                              "keywords": {**_nwc_conv, "ccsd__freeze": 1, "scf__uhf": True},                                         "error": {0: _q3}}, id="ccsd_t_ccsd_  uhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",                              "keywords": {**_p4_conv, "freeze_core": True, "reference": "uhf"},                                     "error": {0: _q20}}, id="ccsd_t_ccsd_  uhf fc: psi4-cc",    marks=using("psi4")),
        # yapf: enable
    ],
)
def test_ccsdpt_prccsd_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsd+t(ccsd)")
    )


#
#   ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                      `-'          `-'                                    `---' `---'
#  <<<  CCSD(T) Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=pytest.mark.ecp),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                                                     }, id="ccsd_t_  rhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                                                 }, id="ccsd_t_  rhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                                                 }, id="ccsd_t_  rhf ae: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                                        "error": {"cc-pvdz": _q45}}, id="ccsd_t_  rhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv},                                                                                                                                      }, id="ccsd_t_  rhf ae: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "qc_module": "ccenergy"},                                                                                                              }, id="ccsd_t_  rhf ae: psi4-cc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "qc_module": "fnocc"},                                                                                                                 }, id="ccsd_t_  rhf ae: psi4-fnocc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "qc_module": "mrcc"},                                                                                                                  }, id="ccsd_t_  rhf ae: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                                                                      }, id="ccsd_t_  rhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ecc"},                                                                                                  }, id="ccsd_t_  rhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "ncc"},                                                                                                  }, id="ccsd_t_  rhf fc: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                                                      "error": {"cc-pvdz": _q45}}, id="ccsd_t_  rhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",               "xptd": {"qc_module": "cc"}, "keywords": {**_nwc_conv, "ccsd__freeze": 1},                                                                                                                   }, id="ccsd_t_  rhf fc: nwchem",        marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "ccenergy"},                                                                                         }, id="ccsd_t_  rhf fc: psi4-cc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "fnocc"},                                                                                            }, id="ccsd_t_  rhf fc: psi4-fnocc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "mrcc"},                                                                                             }, id="ccsd_t_  rhf fc: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                                                                 }, id="ccsd_t_  uhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "ecc"},                                                                                             }, id="ccsd_t_  uhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                                                }, id="ccsd_t_  uhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_nwc_conv, "scf__uhf": True},                                                                                                   "error": {0: _q3}}, id="ccsd_t_  uhf ae: nwchem",        marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "ccenergy"},                                                                                          }, id="ccsd_t_  uhf ae: psi4-cc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                                            "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "mrcc"},                                                                                              }, id="ccsd_t_  uhf ae: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "REFerence": "UHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                                                  }, id="ccsd_t_  uhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_c4_tight, "REFerence": "UHF", "dropmo": [1], "cc_program": "ecc"},                                                                              }, id="ccsd_t_  uhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                                              }, id="ccsd_t_  uhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_nwc_conv, "scf__uhf": True, "ccsd__freeze": 1},                                                                                "error": {0: _q3}}, id="ccsd_t_  uhf fc: nwchem",        marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True, "qc_module": "ccenergy"},                                                                     }, id="ccsd_t_  uhf fc: psi4-cc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",                                            "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True, "qc_module": "mrcc"},                                                                         }, id="ccsd_t_  uhf fc: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                                            "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2},                                                                                }, id="ccsd_t_ rohf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                                            "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "ecc"},                                                                                            }, id="ccsd_t_ rohf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",                                            "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                                              "wrong": {0: _w3}}, id="ccsd_t_ rohf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",                                            "keywords": {**_nwc_conv, "scf__rohf": True},                                                                                                  "error": {0: _q3}}, id="ccsd_t_ rohf ae: nwchem",        marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",                                            "keywords": {**_p4_conv, "reference": "rohf", "qc_module": "ccenergy"},                                                                                         }, id="ccsd_t_ rohf ae: psi4-cc",       marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",                                            "keywords": {**_p4_conv, "reference": "rohf", "qc_module": "mrcc"},                                                                                             }, id="ccsd_t_ rohf ae: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        ## Sum 2021: can get all passing c4 = p4 = mrcc but at cost of CCSD  not matching plain CCSD. Sum 2022: CCSD often has default standard (sd) or semicanonical (sc), while CCSD(T) generally only sc that were matching Sum 2021. sd allows (T) + default CCSD.
        ## TODO express sdsc variants in table. CCSD, too
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sd",                              "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "orbitals": 0, "cc_program": "vcc", "print": 2},                                                  }, id="ccsd_t_ rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sd",                              "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "orbitals": 0, "cc_program": "ecc"},                                                              }, id="ccsd_t_ rohf fc sd: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                                                                 }, id="ccsd_t_ rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "ecc"},                                                                             }, id="ccsd_t_ rohf fc sc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "sdsc": "sd",                              "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                                            "wrong": {0: _w3}}, id="ccsd_t_ rohf fc   : nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc",                                            "keywords": {**_nwc_conv, "scf__rohf": True, "ccsd__freeze": 1},                                                                               "error": {0: _q3}}, id="ccsd_t_ rohf fc   : nwchem",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "qc_module": "ccenergy", "e_convergence": 10, "r_convergence": 9},                           }, id="ccsd_t_ rohf fc sc: psi4-cc",    marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc",                              "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "qc_module": "mrcc", "e_convergence": 10, "r_convergence": 9},                               }, id="ccsd_t_ rohf fc sc: psi4-mrcc",  marks=[*using("psi4"), *using("mrcc")]),
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsd(t)")
    )


#
#   ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'  \ '.   `--'   .' /      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#                                      `-'          `-'
#  <<<  CCSD(T) Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
        pytest.param("def2-svp", ["ch2i2", "ch2i", "ch2i"], id="def2ecp", marks=[pytest.mark.ecp, pytest.mark.long]),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                                                         }, id="ccsd_t_  rhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                                                     }, id="ccsd_t_  rhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                                                     }, id="ccsd_t_  rhf ae: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                                            "error": {"cc-pvdz": _q45}}, id="ccsd_t_  rhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv},                                                                                                                                          }, id="ccsd_t_  rhf ae: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",                        "keywords": {**_p4_conv, "qc_module": "ccenergy"},                                                                                                                  }, id="ccsd_t_  rhf ae: psi4-cc",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc", "print": 2},                                                                                            }, id="ccsd_t_  rhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                                                                                        }, id="ccsd_t_  rhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                                                                       "wrong": {1: _w8}}, id="ccsd_t_  rhf fc: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                                                          "error": {"cc-pvdz": _q45}}, id="ccsd_t_  rhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "ccsd__freeze": 1},                                                                                                                       }, id="ccsd_t_  rhf fc: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",                        "keywords": {**_p4_conv, "freeze_core": True, "qc_module": "ccenergy"},                                                                           "error": {1: _q18}}, id="ccsd_t_  rhf fc: psi4-cc",       marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc", "print": 2},                                                                                     }, id="ccsd_t_  uhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                                                                                 }, id="ccsd_t_  uhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                                                    }, id="ccsd_t_  uhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",                        "keywords": {**_nwc_conv, "scf__uhf": True},                                                                                                       "error": {1: _q3}}, id="ccsd_t_  uhf ae: nwchem-cc",     marks=using("nwchem")),
        # VERY LONG pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",                        "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "ccenergy"},                                                                                   }, id="ccsd_t_  uhf ae: psi4-cc",    marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc", "print": 2},                                                                        }, id="ccsd_t_  uhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "ecc"},                                                                                    }, id="ccsd_t_  uhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                                                  }, id="ccsd_t_  uhf fc: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",                        "keywords": {**_nwc_conv, "scf__uhf": True, "ccsd__freeze": 1},                                                                                    "error": {1: _q3}}, id="ccsd_t_  uhf fc: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "freeze_core": True},                                                                                        }, id="ccsd_t_  uhf fc: psi4",          marks=[*using("psi4"), *using("psi4_112")]),  # psi4 falls back to findif

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "rohf", "cc_program": "vcc", "print": 2},                                                                                    }, id="ccsd_t_ rohf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "rohf", "cc_program": "ecc"},                                                                                                }, id="ccsd_t_ rohf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                                                  "wrong": {1: _w3}}, id="ccsd_t_ rohf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae",                        "keywords": {**_nwc_conv, "scf__rohf": True},                                                                                                      "error": {1: _q3}}, id="ccsd_t_ rohf ae: nwchem-cc",     marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_conv, "reference": "rohf", "qc_module": "ccenergy", "points": 5},                                                                                }, id="ccsd_t_ rohf ae: psi4-cc",       marks=[*using("psi4"), *using("psi4_112")]),
        # Skip ROHF FC until energy resolved

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "cc_program": "vcc", "print": 2},                                                                       }, id="ccsd_t_ rohf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "cc_program": "ecc"},                                                                                   }, id="ccsd_t_ rohf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "rohf", "freeze_core": True, "qc_module": "ccenergy", "function_kwargs": {"dertype": 0}},                           }, id="ccsd_t_ rohf fc: psi4-ccenergy", marks=[*using("psi4"), *using("psi4_112")]),
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccsd(t)")
    )


#
#   ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.      ,--.  ,--.                     ,--.
#  '  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-----' `-------'  \ '.   `--'   .' /     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#                                      `-'          `-'
#  <<<  CCSD(T) Hessian


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ncc errors
        # * nwchem-cc Error in Communication
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                                                    }, id="ccsd_t_  rhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                                }, id="ccsd_t_  rhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                       "error": {"cc-pvdz": _q45}}, id="ccsd_t_  rhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "qc_module": "ccenergy", "function_kwargs": {"dertype": 0}},                                                }, id="ccsd_t_  rhf ae: psi4-cc",    marks=[*using("psi4"), *using("psi4_112")]),  # FD SWITCH

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc", "print": 2},                                                                       }, id="ccsd_t_  rhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                                                                   }, id="ccsd_t_  rhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                                     "error": {"cc-pvdz": _q45}}, id="ccsd_t_  rhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "freeze_core": True, "qc_module": "ccenergy", "function_kwargs": {"dertype": 0}},                           }, id="ccsd_t_  rhf fc: psi4-cc",    marks=[*using("psi4"), *using("psi4_112")]),  # FD SWITCH

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc", "print": 2},                                                                }, id="ccsd_t_  uhf ae: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                                                            }, id="ccsd_t_  uhf ae: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                                               }, id="ccsd_t_  uhf ae: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",                        "keywords": {**_nwc_conv, "scf__uhf": True},                                                                                  "error": {2: _q3}}, id="ccsd_t_  uhf ae: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "qc_module": "ccenergy", "function_kwargs": {"dertype": 0}},                            }, id="ccsd_t_  uhf ae: psi4-cc",    marks=[*using("psi4"), *using("psi4_112")]),  # FD SWITCH

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc", "print": 2},                                                   }, id="ccsd_t_  uhf fc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "ecc"},                                                               }, id="ccsd_t_  uhf fc: cfour-ecc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                                             }, id="ccsd_t_  uhf fc: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",                        "keywords": {**_nwc_conv, "scf__uhf": True, "ccsd__freeze": 1},                                                               "error": {2: _q3}}, id="ccsd_t_  uhf fc: nwchem-cc",  marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc", "xptd": {"fd": False}, "keywords": {**_p4_conv, **_p4_fd, "reference": "uhf", "freeze_core": True, "function_kwargs": {"dertype": 0}},                                }, id="ccsd_t_  uhf fc: psi4",       marks=[*using("psi4"), *using("psi4_112")]),  # FD SWITCH
        # yapf: enable
    ],
)
def test_ccsd_prt_pr_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "ccsd(t)")
    )


#
#                  ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.      ,------.
#   ,--,--.,-----.'  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  ' ,-.  |'-----'|  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  \ '-'  |       '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `--`--'        `-----' `-----'`-----' `-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                     `-'          `-'                                    `---' `---'
#  <<<  a-CCSD(T) Energy aka Lambda-CCSD(T) aka CCSD(aT) aka CCSD(T)_L


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # vcc stops prematurely
        #pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc"},                                           }, id="a-ccsd_t_  rhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                           }, id="a-ccsd_t_  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                           }, id="a-ccsd_t_  rhf ae: cfour-ncc", marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv},                                              "error": {0: _q27}}, id="a-ccsd_t_  rhf ae: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "ccenergy"},                                        }, id="a-ccsd_t_  rhf ae: psi4-cc",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                              }, id="a-ccsd_t_  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                              }, id="a-ccsd_t_  rhf fc: cfour-ncc", marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, "qc_module": "ccenergy", "freeze_core": True},                   }, id="a-ccsd_t_  rhf fc: psi4-cc",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},     "error": {0: _q28}}, id="a-ccsd_t_  uhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, "reference": "uhf"},                           "error": {0: _q29}}, id="a-ccsd_t_  uhf ae: psi4",      marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, "reference": "uhf", "qc_module": "mrcc"},                        }, id="a-ccsd_t_  uhf ae: psi4-mrcc", marks=[*using("psi4"), *using("mrcc")]),
        # yapf: enable
    ],
)
def test_accsd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "a-ccsd(t)")
    )


#                  ,-----. ,-----. ,---.  ,------.    ,-.,--------.,-.       ,----.                     ,--.,--.                 ,--.
#   ,--,--.,-----.'  .--./'  .--./'   .-' |  .-.  \  / .''--.  .--''. \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  ' ,-.  |'-----'|  |    |  |    `.  `-. |  |  \  :|  |    |  |    |  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  \ '-'  |       '  '--'\'  '--'\.-'    ||  '--'  /|  |    |  |    |  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `--`--'        `-----' `-----'`-----' `-------'  \ '.   `--'   .' /      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#                                                     `-'          `-'
#  <<<  a-CCSD(T) Gradient aka Lambda-CCSD(T) aka CCSD(aT) aka CCSD(T)_L


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},              "wrong": {1: _w25}}, id="a-ccsd_t_  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                }, id="a-ccsd_t_  rhf ae: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"}, "wrong": {1: _w25}}, id="a-ccsd_t_  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},  "wrong": {1: _w8}}, id="a-ccsd_t_  rhf fc: cfour-ncc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_accsd_prt_pr_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "a-ccsd(t)")
    )


#
#  ,-----.   ,-----. ,-----.,------.    ,-.,--------.,-.      ,------.
#  |  |) /_ '  .--./'  .--./|  .-.  \  / .''--.  .--''. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \|  |    |  |    |  |  \  :|  |    |  |    |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' /'  '--'\'  '--'\|  '--'  /|  |    |  |    |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------'  `-----' `-----'`-------'  \ '.   `--'   .' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                       `-'          `-'                                    `---' `---'
#  <<<  BCCD(T) Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2, "orbitals": 1},                     }, id="bccd_t_  rhf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv},                                                                      }, id="bccd_t_  rhf ae: psi4",         marks=using("psi4")),

        # cfour fc Brueckner hangs!
        # pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                     }, id="bccd_t_  rhf fc: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",               "keywords": {**_p4_conv, "freeze_core": True},                                                 }, id="bccd_t_  rhf fc: psi4",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2, "orbitals": 1}, }, id="bccd_t_  uhf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",               "keywords": {**_p4_conv, "reference": "uhf"},                                                  }, id="bccd_t_  uhf ae: psi4",         marks=using("psi4")),

        # pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",               "keywords": {**_c4_tight, "REFerence": "uHF", "dropmo": [1], "cc_program": "vcc", "print": 2}, }, id="bccd_t_  uhf fc: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",               "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True},                             }, id="bccd_t_  uhf fc: psi4",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2, "orbitals": 1},}, id="bccd_t_ rohf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",               "keywords": {**_p4_conv, "reference": "rohf"},                                                 }, id="bccd_t_ rohf ae: psi4",         marks=using("psi4")),  # NOTE: sc

        # pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "print": 2},}, id="bccd_t_ rohf fc sc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True},                            }, id="bccd_t_ rohf fc sc: psi4",      marks=using("psi4")),

        # cfour: rhf & uhf fail w/o "cfour_orbitals": 1
        # yapf: enable
    ],
)
def test_bccd_prt_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "bccd(t)")
    )


#
#   ,-----. ,-----.,----.     ,------.
#  '  .--./'  .--./'.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |      .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\/'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                           `---' `---'
#  <<<  CC3 Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc, ncc doesn't run
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "vcc", "print": 2},                                                       }, id="cc3  rhf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv},                                                                                         }, id="cc3  rhf ae: psi4",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": [1], "cc_program": "vcc", "print": 2},                                        }, id="cc3  rhf fc: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc",               "keywords": {**_p4_conv, "freeze_core": True},                                                                    }, id="cc3  rhf fc: psi4",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "UHF", "cc_program": "vcc", "print": 2},                                   }, id="cc3  uhf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae",               "keywords": {**_p4_conv, "reference": "uhf"},                                                                     }, id="cc3  uhf ae: psi4",         marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",               "keywords": {**_c4_tight, "REFerence": "uHF", "dropmo": [1], "cc_program": "vcc", "print": 2},                    }, id="cc3  uhf fc: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "fc",               "keywords": {**_p4_conv, "reference": "uhf", "freeze_core": True},                                                }, id="cc3  uhf fc: psi4",         marks=using("psi4")),

        # cfour paper disavows cc3 rohf
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",               "keywords": {**_c4_tight, "REFerence": "roHF", "cc_program": "vcc", "print": 2},                "error": {0: _q46}}, id="cc3 rohf ae: cfour-vcc",    marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae",               "keywords": {**_p4_conv, "reference": "rohf"},                                                                    }, id="cc3 rohf ae: psi4",         marks=using("psi4")),  # NOTE: sc

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "REFerence": "roHF", "dropmo": [1], "cc_program": "vcc", "print": 2}, "error": {0: _q46}}, id="cc3 rohf fc sc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True},                                               }, id="cc3 rohf fc sc: psi4",      marks=using("psi4")),
        # yapf: enable
    ],
)
def test_cc3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "cc3")
    )


#
#   ,-----. ,-----.,----.      ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'.-.  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |      .' <     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\/'-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CC3 Gradient


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.            ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   | ,--,--.    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  |' ,-.  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  |\ '-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'           `--' `--`--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                                   `---' `---'
#  <<<  CCSDT-1a Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc"},                                                    }, id="ccsdt-1a  rhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                    }, id="ccsdt-1a  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                    }, id="ccsdt-1a  rhf ae: cfour-ncc", marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv},                                                       "error": {0: _q24}}, id="ccsdt-1a  rhf ae: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv},                                                        "error": {0: _q47}}, id="ccsdt-1a  rhf ae: psi4",      marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "mrcc"},                                                     }, id="ccsdt-1a  rhf ae: psi4-mrcc", marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc"},                                       }, id="ccsdt-1a  rhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                       }, id="ccsdt-1a  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                       }, id="ccsdt-1a  rhf fc: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc"},                                }, id="ccsdt-1a  uhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                }, id="ccsdt-1a  uhf ae: cfour-ecc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc"},                   }, id="ccsdt-1a  uhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "ecc"},                   }, id="ccsdt-1a  uhf fc: cfour-ecc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1a_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsdt-1a")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.             ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   | ,--,--.    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  |' ,-.  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  |\ '-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'           `--' `--`--'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSDT-1a Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * vcc turns into ecc
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                   }, id="ccsdt-1a  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                   }, id="ccsdt-1a  rhf ae: cfour-ncc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "function_kwargs": {"dertype": 0}},                                        }, id="ccsdt-1a  rhf ae: psi4-mrcc",),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                                    "wrong": {1: _w16}}, id="ccsdt-1a  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                                    "wrong": {1: _w22}}, id="ccsdt-1a  rhf fc: cfour-ncc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "freeze_core": True, "function_kwargs": {"dertype": 0}},                   }, id="ccsdt-1a  rhf fc: psi4-mrcc",),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                  "error": {1: _q19}}, id="ccsdt-1a  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                 "error": {1: _q19}}, id="ccsdt-1a rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1a_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccsdt-1a")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.            ,--.  ,--.                     ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   | ,--,--.    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  |' ,-.  |    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  |\ '-'  |    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-----' `-------'   `--'           `--' `--`--'    `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  CCSDT-1a Hessian


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                   }, id="ccsdt-1a  rhf ae: cfour-ecc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "function_kwargs": {"dertype": 0}},                                        }, id="ccsdt-1a  rhf ae: psi4-mrcc",),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                                    "wrong": {2: _w17}}, id="ccsdt-1a  rhf fc: cfour-ecc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "freeze_core": True, "function_kwargs": {"dertype": 0}},                   }, id="ccsdt-1a  rhf fc: psi4-mrcc",),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                  "error": {2: _q19}}, id="ccsdt-1a  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                 "error": {2: _q19}}, id="ccsdt-1a rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1a_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "ccsdt-1a")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.,--.       ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   ||  |-.     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  || .-. '    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  || `-' |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'           `--' `---'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                                  `---' `---'
#  <<<  CCSDT-1b Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc"},                                                    }, id="ccsdt-1b  rhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                    }, id="ccsdt-1b  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                    }, id="ccsdt-1b  rhf ae: cfour-ncc", marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv},                                                       "error": {0: _q25}}, id="ccsdt-1b  rhf ae: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv},                                                        "error": {0: _q47}}, id="ccsdt-1b  rhf ae: psi4",      marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "mrcc"},                                                     }, id="ccsdt-1b  rhf ae: psi4-mrcc", marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc"},                                       }, id="ccsdt-1b  rhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                       }, id="ccsdt-1b  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                       }, id="ccsdt-1b  rhf fc: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc"},                                }, id="ccsdt-1b  uhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                }, id="ccsdt-1b  uhf ae: cfour-ecc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc"},                   }, id="ccsdt-1b  uhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "ecc"},                   }, id="ccsdt-1b  uhf fc: cfour-ecc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1b_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsdt-1b")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.,--.        ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   ||  |-.     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  || .-. '    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  || `-' |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'           `--' `---'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSDT-1b Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * vcc turns into ecc
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                   }, id="ccsdt-1b  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                   }, id="ccsdt-1b  rhf ae: cfour-ncc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "function_kwargs": {"dertype": 0}},                                        }, id="ccsdt-1b  rhf ae: psi4-mrcc",),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                                    "wrong": {1: _w18}}, id="ccsdt-1b  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                                    "wrong": {1: _w23}}, id="ccsdt-1b  rhf fc: cfour-ncc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "freeze_core": True, "function_kwargs": {"dertype": 0}},                   }, id="ccsdt-1b  rhf fc: psi4-mrcc",),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                  "error": {1: _q19}}, id="ccsdt-1b  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                 "error": {1: _q19}}, id="ccsdt-1b rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1b_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccsdt-1b")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,--.,--.       ,--.  ,--.                     ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----./   ||  |-.     |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'`|  || .-. '    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |           |  || `-' |    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-----' `-------'   `--'           `--' `---'     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  CCSDT-1b Hessian


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * adz ae the two debugs don't match to 1e-5 - went with p4c4 to pass
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                   }, id="ccsdt-1b  rhf ae: cfour-ecc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "function_kwargs": {"dertype": 0}},                                        }, id="ccsdt-1b  rhf ae: psi4-mrcc",),

        # * fc is wrong wrt the two debugs for adz - commented b/c not filling in ref for other bases
        # pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "cc_program": "ecc", "dropmo": 1},                                                    "wrong": {2: _w17}}, id="ccsdt-1b  rhf fc: cfour-ecc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "freeze_core": True, "function_kwargs": {"dertype": 0}},                   }, id="ccsdt-1b  rhf fc: psi4-mrcc",),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                  "error": {2: _q19}}, id="ccsdt-1b  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                 "error": {2: _q19}}, id="ccsdt-1b rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt1b_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "ccsdt-1b")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,---.     ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  \    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----' .-' .'    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /   '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'          '-----'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                             `---' `---'
#  <<<  CCSDT-2 Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc"},                                                    }, id="ccsdt-2  rhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                    }, id="ccsdt-2  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                    }, id="ccsdt-2  rhf ae: cfour-ncc", marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv},                                                       "error": {0: _q26}}, id="ccsdt-2  rhf ae: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv},                                                        "error": {0: _q26}}, id="ccsdt-2  rhf ae: psi4",      marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc"},                                       }, id="ccsdt-2  rhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                       }, id="ccsdt-2  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                       }, id="ccsdt-2  rhf fc: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc"},                                }, id="ccsdt-2  uhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                }, id="ccsdt-2  uhf ae: cfour-ecc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc"},                   }, id="ccsdt-2  uhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "ecc"},                   }, id="ccsdt-2  uhf fc: cfour-ecc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt2_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsdt-2")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,---.      ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  \    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----' .-' .'    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /   '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'          '-----'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSDT-2 Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ccsdt-2 not available in mrcc
        # * vcc turns into ecc
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                }, id="ccsdt-2  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                }, id="ccsdt-2  rhf ae: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"}, "wrong": {1: _w20}}, id="ccsdt-2  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"}, "wrong": {1: _w24}}, id="ccsdt-2  rhf fc: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},               "error": {1: _q19}}, id="ccsdt-2  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},              "error": {1: _q19}}, id="ccsdt-2 rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt2_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccsdt-2")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.        ,---.     ,--.  ,--.                     ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  \    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----' .-' .'    |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /   '-.    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-----' `-------'   `--'          '-----'    `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  CCSDT-2 Hessian


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ccsdt-2 not available in mrcc
        # * adz has 7e-6 sized discrepancies from debug - adding false fd marking to pass it
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "xptd": {"fd": False}, "keywords": {**_c4_tight, "cc_program": "ecc"},                                }, id="ccsdt-2  rhf ae: cfour-ecc", marks=using("cfour")),  # TODO FD

        # * no use testing hess if grad wrong
        # pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",                        "keywords": {**_c4_tight, "cc_program": "ecc", "dropmo": 1}, "wrong": {2: _w17}}, id="ccsdt-2  rhf fc: cfour-ecc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "uhf"},               "error": {2: _q19}}, id="ccsdt-2  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae",                        "keywords": {**_c4_tight, "reference": "rohf"},              "error": {2: _q19}}, id="ccsdt-2 rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt2_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "ccsdt-2")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.       ,----.     ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  |    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'  .' <     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /'-'  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'          `----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                             `---' `---'
#  <<<  CCSDT-3 Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "vcc"},                                                    }, id="ccsdt-3  rhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                    }, id="ccsdt-3  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                    }, id="ccsdt-3  rhf ae: cfour-ncc", marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv},                                                       "error": {0: _q23}}, id="ccsdt-3  rhf ae: nwchem",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv},                                                        "error": {0: _q47}}, id="ccsdt-3  rhf ae: psi4",      marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, "qc_module": "mrcc"},                                                     }, id="ccsdt-3  rhf ae: psi4-mrcc", marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc"},                                       }, id="ccsdt-3  rhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                       }, id="ccsdt-3  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                       }, id="ccsdt-3  rhf fc: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc"},                                }, id="ccsdt-3  uhf ae: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "ecc"},                                }, id="ccsdt-3  uhf ae: cfour-ecc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc"},                   }, id="ccsdt-3  uhf fc: cfour-vcc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc", "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "ecc"},                   }, id="ccsdt-3  uhf fc: cfour-ecc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt3_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsdt-3")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.       ,----.      ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  |    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'  .' <     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /'-'  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'          `----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSDT-3 Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * vcc turns into ecc
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                   }, id="ccsdt-3  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                                   }, id="ccsdt-3  rhf ae: cfour-ncc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "function_kwargs": {"dertype": 0}},                                        }, id="ccsdt-3  rhf ae: psi4-mrcc",),

        # p4c4 by psi findif or mrcc by psi findif differ from ecc analytic by 1e-6
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                                    "wrong": {1: _w14}}, id="ccsdt-3  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                                    "wrong": {1: _w13}}, id="ccsdt-3  rhf fc: cfour-ncc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "freeze_core": True, "function_kwargs": {"dertype": 0}},                   }, id="ccsdt-3  rhf fc: psi4-mrcc",),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                  "error": {1: _q19}}, id="ccsdt-3  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                 "error": {1: _q19}}, id="ccsdt-3 rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt3_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccsdt-3")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.       ,----.     ,--.  ,--.                     ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--',-----.'.-.  |    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    `.  `-. |  |  \  :  |  |   '-----'  .' <     |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |          /'-'  |    |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-----' `-------'   `--'          `----'     `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  CCSDT-3 Hessian


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                                   }, id="ccsdt-3  rhf ae: cfour-ecc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "function_kwargs": {"dertype": 0}},                                        }, id="ccsdt-3  rhf ae: psi4-mrcc",),

        # p4c4 by psi findif or mrcc by psi findif differ from ecc analytic by 1e-4
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                                    "wrong": {2: _w15}}, id="ccsdt-3  rhf fc: cfour-ecc", marks=using("cfour")),
        # DEBUG pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "fc", "keywords": {**_p4_conv, **_p4_fd, "reference": "rhf", "freeze_core": True, "function_kwargs": {"dertype": 0}},                   }, id="ccsdt-3  rhf fc: psi4-mrcc",),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                  "error": {2: _q19}}, id="ccsdt-3  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                 "error": {2: _q19}}, id="ccsdt-3 rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt3_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "ccsdt-3")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.    ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |       |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |       |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'       `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                               `---' `---'
#  <<<  CCSDT Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "vcc"},                                                                            }, id="ccsdt  rhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "ecc"},                                                                            }, id="ccsdt  rhf ae: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae",               "keywords": {**_c4_tight, "cc_program": "ncc"},                                                                            }, id="ccsdt  rhf ae: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "qc_module": "tce"},                                                   "error": {"cc-pvdz": _q45}}, id="ccsdt  rhf ae: nwchem-tce",    marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv},                                                                                "error": {0: _q47}}, id="ccsdt  rhf ae: psi4",          marks=using("psi4")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae",               "keywords": {**_p4_conv, "qc_module": "mrcc"},                                                                             }, id="ccsdt  rhf ae: psi4-mrcc",     marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "vcc"},                                                               }, id="ccsdt  rhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                                                               }, id="ccsdt  rhf fc: cfour-ecc",     marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc",               "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},                                                               }, id="ccsdt  rhf fc: cfour-ncc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "tce__freeze": 1, "qc_module": "tce"},                                 "error": {"cc-pvdz": _q45}}, id="ccsdt  rhf fc: nwchem-tce",    marks=using("nwchem")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae",               "keywords": {**_c4_tight, "reference": "uhf", "cc_program": "vcc"},                                                        }, id="ccsdt  uhf ae: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae",               "keywords": {**_nwc_conv, "scf__uhf": True, "qc_module": "tce"},                                                           }, id="ccsdt  uhf ae: nwchem-tce",    marks=using("nwchem")),

        # cfour uhf/rohf ecc does not converge, ncc does not run
        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "fc",               "keywords": {**_c4_tight, "reference": "uhf", "dropmo": 1, "cc_program": "vcc"},                                           }, id="ccsdt  uhf fc: cfour-vcc",     marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc",               "keywords": {**_nwc_conv, "scf__uhf": True, "tce__freeze": 1, "qc_module": "tce"},                                         }, id="ccsdt  uhf fc: nwchem-tce",    marks=using("nwchem")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_c4_tight, "reference": "rohf", "cc_program": "vcc"},                                                       }, id="ccsdt rohf ae sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_c4_tight, "reference": "rohf", "orbitals": 0, "cc_program": "vcc"},                                        }, id="ccsdt rohf ae sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "sdsc": "sd", "keywords": {**_nwc_conv, "scf__rohf": True, "qc_module": "tce"},                                                          }, id="ccsdt rohf ae sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "sdsc": "sc", "keywords": {**_p4_conv, "reference": "rohf", "qc_module": "mrcc"},                                                        }, id="ccsdt rohf ae sc: psi4-mrcc",  marks=[*using("psi4"), *using("mrcc")]),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "cc_program": "vcc"},                                          }, id="ccsdt rohf fc sc: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_c4_tight, "reference": "rohf", "dropmo": 1, "orbitals": 0, "cc_program": "vcc"},                           }, id="ccsdt rohf fc sd: cfour-vcc",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "fc", "sdsc": "sd", "keywords": {**_nwc_conv, "scf__rohf": True, "tce__freeze": 1, "qc_module": "tce"},                                        }, id="ccsdt rohf fc sd: nwchem-tce", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "fc", "sdsc": "sc", "keywords": {**_p4_conv, "reference": "rohf", "freeze_core": True, "qc_module": "mrcc"},                                   }, id="ccsdt rohf fc sc: psi4-mrcc",  marks=[*using("psi4"), *using("mrcc")]),
        # yapf: enable
    ],
)
def test_ccsdt_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsdt")
    )


#   ,-----. ,-----. ,---.  ,------. ,--------.     ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |       |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |       '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'        `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSDT Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * vcc turns into ecc
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                }, id="ccsdt  rhf ae: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ncc"},                                }, id="ccsdt  rhf ae: cfour-ncc", marks=using("cfour")),
        # VLONG pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "xptd": {"fd": True},  "keywords": {**_nwc_conv, "qc_module": "tce"},                                                                                 }, id="ccsdt  rhf ae: nwchem-tce", marks=using("nwchem")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                   }, id="ccsdt  rhf fc: cfour-ecc", marks=using("cfour")),
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ncc"},  "wrong": {1: _w7}}, id="ccsdt  rhf fc: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},               "error": {1: _q19}}, id="ccsdt  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},              "error": {1: _q19}}, id="ccsdt rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccsdt")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------.    ,--.  ,--.                     ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'    |  '--'  | ,---.  ,---.  ,---. `--' ,--,--.,--,--,
#  |  |    |  |    `.  `-. |  |  \  :  |  |       |  .--.  || .-. :(  .-' (  .-' ,--.' ,-.  ||      \
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |       |  |  |  |\   --..-'  `).-'  `)|  |\ '-'  ||  ||  |
#   `-----' `-----'`-----' `-------'   `--'       `--'  `--' `----'`----' `----' `--' `--`--'`--''--'
#
#  <<<  CCSDT Hessian


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(2, id="hes2"),
        # pytest.param(1, id="hes1", marks=pytest.mark.long),
        # pytest.param(0, id="hes0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # * ncc errors
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight, "cc_program": "ecc"},                                }, id="ccsdt  rhf ae: cfour-ecc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {**_c4_tight, "dropmo": 1, "cc_program": "ecc"},                   }, id="ccsdt  rhf fc: cfour-ecc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},               "error": {2: _q19}}, id="ccsdt  uhf ae: cfour",     marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},              "error": {2: _q19}}, id="ccsdt rohf ae: cfour",     marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt_hessian_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "hessian", "ccsdt")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------. ,-. ,-----.   ,-.      ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'/ .''  .-.  '  '. \     |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |  |  | |  | |  |   |  |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |  |  | '  '-'  '-. |  |    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'   \ '. `-----'--'.' /     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                              `-'           `-'                                    `---' `---'
#  <<<  CCSDT(Q) Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {"cc_program": "ncc"},             }, id="ccsdt_q_  rhf ae: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {"dropmo": 1, "cc_program": "ncc"},}, id="ccsdt_q_  rhf fc: cfour-ncc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt_prq_pr_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsdt(q)")
    )


#   ,-----. ,-----. ,---.  ,------. ,--------. ,-. ,-----.   ,-.       ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--'/ .''  .-.  '  '. \     '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |  |  | |  | |  |   |  |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |  |  | '  '-'  '-. |  |    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'   \ '. `-----'--'.' /      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#                                              `-'           `-'
#  <<<  CCSDT(Q) Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {"cc_program": "ncc"},                                }, id="ccsdt_q_  rhf ae: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {"dropmo": 1, "cc_program": "ncc"}, "wrong": {1: _w7 }}, id="ccsdt_q_  rhf fc: cfour-ncc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdt_prq_pr_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccsdt(q)")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------. ,-----.       ,------.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--''  .-.  '      |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   |  | |  |      |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |   '  '-'  '-.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#   `-----' `-----'`-----' `-------'   `--'    `-----'--'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                          `---' `---'
#  <<<  CCSDTQ Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz"),
        # pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # ecc/vcc sent to cfour, but it switches to ncc
        # adz only converges to cc_conv=9
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {"cc_program": "ncc"},             }, id="ccsdtq  rhf ae: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {"dropmo": 1, "cc_program": "ncc"},}, id="ccsdtq  rhf fc: cfour-ncc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdtq_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "ccsdtq")
    )


#
#   ,-----. ,-----. ,---.  ,------. ,--------. ,-----.        ,----.                     ,--.,--.                 ,--.
#  '  .--./'  .--./'   .-' |  .-.  \'--.  .--''  .-.  '      '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  |    |  |    `.  `-. |  |  \  :  |  |   |  | |  |      |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  '  '--'\'  '--'\.-'    ||  '--'  /  |  |   '  '-'  '-.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#   `-----' `-----'`-----' `-------'   `--'    `-----'--'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  CCSDTQ Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        # vcc/ecc error out
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {"cc_program": "ncc"},                                }, id="ccsdtq  rhf ae: cfour-ncc", marks=using("cfour")),

        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "fc", "keywords": {"dropmo": 1, "cc_program": "ncc"}, "wrong": {1: _w6 }}, id="ccsdtq  rhf fc: cfour-ncc", marks=using("cfour")),
        # yapf: enable
    ],
)
def test_ccsdtq_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "ccsdtq")
    )


#
#  ,------. ,-----.  ,------.    ,------.
#  |  .--. '|  |) /_ |  .---'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  '--' ||  .-.  \|  `--,     |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  | --' |  '--' /|  `---.    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `--'     `------' `------'    `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                              `---' `---'
#  <<<  PBE Energy

_psi_grid = {"dft_radial_points": 99, "dft_spherical_points": 590}
_nwc_grid = {"dft__grid__lebedev": [99, 14]}


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                                             "error": {0: _q35}}, id="pbe  rhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid},                                                                                                                                  }, id="pbe  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "scf_type": "pk"},                                                                                                                 }, id="pbe  rhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                                                         "error": {0: _q35}}, id="pbe  uhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__uhf": True},                                                                                                                }, id="pbe  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "uhf", "scf_type": "pk"},                                                                                             }, id="pbe  uhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                                                        "error": {0: _q35}}, id="pbe rohf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True},                                                                                           }, id="pbe rohf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "rohf", "scf_type": "pk"},                                                                          "error": {0: _q36}}, id="pbe rohf ae: psi4",   marks=using("psi4")),
        # pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True, "dft__convergence__energy": 1.e-9, "dft__convergence__gradient": 1.e-7},                   }, id="pbe rohf ae: nwchem", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_pbe_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "pbe")
    )


#
#  ,------. ,-----.  ,------.     ,----.                     ,--.,--.                 ,--.
#  |  .--. '|  |) /_ |  .---'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  '--' ||  .-.  \|  `--,     |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  | --' |  '--' /|  `---.    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `--'     `------' `------'     `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  PBE Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                                             "error": {1: _q37}}, id="pbe  rhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid},                                                                                                                                  }, id="pbe  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "scf_type": "pk"},                                                                                                                 }, id="pbe  rhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                                                         "error": {1: _q37}}, id="pbe  uhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__uhf": True},                                                                                                                }, id="pbe  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "uhf", "scf_type": "pk"},                                                                                             }, id="pbe  uhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                                                        "error": {1: _q37}}, id="pbe rohf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True, "dft__convergence__gradient": 1.e-5},                                                      }, id="pbe rohf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "rohf", "scf_type": "pk"},                                                                          "error": {1: _q36}}, id="pbe rohf ae: psi4",   marks=using("psi4")),
        # pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True, "dft__convergence__energy": 1.e-9, "dft__convergence__gradient": 1.e-7},                   }, id="pbe rohf ae: nwchem", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_pbe_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "pbe")
    )


#
#  ,-----.  ,----. ,--.,--.   ,--.,------.     ,------.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' //'-'  ||  '--. |  |   |  | --'     |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------' `----' `-----' `--'   `--'         `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                            `---' `---'
#  <<<  B3LYP Energy


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                                             "error": {0: _q35}}, id="b3lyp  rhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid},                                                                                                                                  }, id="b3lyp  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "scf_type": "pk"},                                                                                                                 }, id="b3lyp  rhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                                                         "error": {0: _q35}}, id="b3lyp  uhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__uhf": True},                                                                                                                }, id="b3lyp  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "uhf", "scf_type": "pk"},                                                                                             }, id="b3lyp  uhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                                                        "error": {0: _q35}}, id="b3lyp rohf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True},                                                                                           }, id="b3lyp rohf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "rohf", "scf_type": "pk"},                                                                          "error": {0: _q36}}, id="b3lyp rohf ae: psi4",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True, "dft__convergence__energy": 1.e-9, "dft__convergence__gradient": 1.e-7},                   }, id="b3lyp rohf ae: nwchem", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b3lyp_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "b3lyp")
    )


#
#  ,-----.  ,----. ,--.,--.   ,--.,------. ,-----.    ,------.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '|  .--'    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |'--. `\    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' //'-'  ||  '--. |  |   |  | --' .--'  /    |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------' `----' `-----' `--'   `--'     `----'     `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                   `---' `---'
#  <<<  B3LYP5 Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                                             "error": {0: _q35}}, id="b3lyp5  rhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid},                                                                                                                                  }, id="b3lyp5  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "scf_type": "pk"},                                                                                                                 }, id="b3lyp5  rhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                                                         "error": {0: _q35}}, id="b3lyp5  uhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__uhf": True},                                                                                                                }, id="b3lyp5  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "uhf", "scf_type": "pk"},                                                                                             }, id="b3lyp5  uhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                                                        "error": {0: _q35}}, id="b3lyp5 rohf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True},                                                                                           }, id="b3lyp5 rohf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "rohf", "scf_type": "pk"},                                                                          "error": {0: _q36}}, id="b3lyp5 rohf ae: psi4",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True, "dft__convergence__energy": 1.e-9, "dft__convergence__gradient": 1.e-7},                   }, id="b3lyp5 rohf ae: nwchem", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b3lyp5_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "b3lyp5")
    )


#
#  ,-----.  ,----. ,--.,--.   ,--.,------.      ,----.                     ,--.,--.                 ,--.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--' //'-'  ||  '--. |  |   |  | --'     '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `------' `----' `-----' `--'   `--'          `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  B3LYP Gradient


@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                                             "error": {1: _q37}}, id="b3lyp  rhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid},                                                                                                                                  }, id="b3lyp  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "scf_type": "pk"},                                                                                                                 }, id="b3lyp  rhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                                                         "error": {1: _q37}}, id="b3lyp  uhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__uhf": True},                                                                                                                }, id="b3lyp  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "uhf", "scf_type": "pk"},                                                                                             }, id="b3lyp  uhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                                                        "error": {1: _q37}}, id="b3lyp rohf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True, "dft__convergence__gradient": 1.e-5},                                                      }, id="b3lyp rohf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "rohf", "scf_type": "pk"},                                                                          "error": {1: _q36}}, id="b3lyp rohf ae: psi4",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True, "dft__convergence__energy": 1.e-9, "dft__convergence__gradient": 1.e-7},                   }, id="b3lyp rohf ae: nwchem", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b3lyp_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "b3lyp")
    )


#
#  ,-----.  ,----. ,--.,--.   ,--.,------. ,-----.     ,----.                     ,--.,--.                 ,--.
#  |  |) /_ '.-.  ||  | \  `.'  / |  .--. '|  .--'    '  .-./   ,--.--. ,--,--. ,-|  |`--' ,---. ,--,--, ,-'  '-.
#  |  .-.  \  .' < |  |  '.    /  |  '--' |'--. `\    |  | .---.|  .--'' ,-.  |' .-. |,--.| .-. :|      \'-.  .-'
#  |  '--' //'-'  ||  '--. |  |   |  | --' .--'  /    '  '--'  ||  |   \ '-'  |\ `-' ||  |\   --.|  ||  |  |  |
#  `------' `----' `-----' `--'   `--'     `----'      `------' `--'    `--`--' `---' `--' `----'`--''--'  `--'
#
#  <<<  B3LYP5 Gradient


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(1, id="grd1"),
        # pytest.param(0, id="grd0", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                                                                                             "error": {1: _q37}}, id="b3lyp5  rhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid},                                                                                                                                  }, id="b3lyp5  rhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "scf_type": "pk"},                                                                                                                 }, id="b3lyp5  rhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "uhf",  "fcae": "ae", "keywords": {**_c4_tight, "reference": "uhf"},                                                                                                         "error": {1: _q37}}, id="b3lyp5  uhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__uhf": True},                                                                                                                }, id="b3lyp5  uhf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "uhf", "scf_type": "pk"},                                                                                             }, id="b3lyp5  uhf ae: psi4",   marks=using("psi4")),

        pytest.param({"call": "cfour",  "reference": "rohf", "fcae": "ae", "keywords": {**_c4_tight, "reference": "rohf"},                                                                                                        "error": {1: _q37}}, id="b3lyp5 rohf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True, "dft__convergence__gradient": 1.e-5},                                                      }, id="b3lyp5 rohf ae: nwchem", marks=using("nwchem")),
        pytest.param({"call": "psi4",   "reference": "rohf", "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "rohf", "scf_type": "pk"},                                                                          "error": {1: _q36}}, id="b3lyp5 rohf ae: psi4",   marks=using("psi4")),
        # DEBUG pytest.param({"call": "nwchem", "reference": "rohf", "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__rohf": True, "dft__rodft": True, "dft__convergence__energy": 1.e-9, "dft__convergence__gradient": 1.e-7},                   }, id="b3lyp5 rohf ae: nwchem", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b3lyp5_gradient_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "gradient", "b3lyp5")
    )


#
#  ,-----.   ,---. ,------. ,--.,--.   ,--.,------.     ,------.
#  |  |) /_ '.-.  \|  .--. '|  | \  `.'  / |  .--. '    |  .---',--,--,  ,---. ,--.--. ,---.,--. ,--.
#  |  .-.  \ .-' .'|  '--' ||  |  '.    /  |  '--' |    |  `--, |      \| .-. :|  .--'| .-. |\  '  /
#  |  '--' //   '-.|  | --' |  '--. |  |   |  | --'     |  `---.|  ||  |\   --.|  |   ' '-' ' \   '
#  `------' '-----'`--'     `-----' `--'   `--'         `------'`--''--' `----'`--'   .`-  /.-'  /
#                                                                                     `---' `---'
#  <<<  B2PLYP Energy


@pytest.mark.nonroutine
@pytest.mark.parametrize(
    "dertype",
    [
        pytest.param(0, id="ene0"),
    ],
)
@pytest.mark.parametrize(
    "basis, subjects",
    [
        pytest.param("cc-pvdz", ["hf", "bh3p", "bh3p"], id="dz"),
        pytest.param("aug-cc-pvdz", ["h2o", "nh2", "nh2"], id="adz", marks=pytest.mark.long),
        pytest.param("cfour-qz2p", ["h2o", "nh2", "nh2"], id="qz2p", marks=pytest.mark.long),
    ],
)
@pytest.mark.parametrize(
    "inp",
    [
        # yapf: disable
        pytest.param({"call": "cfour",  "reference": "rhf",  "fcae": "ae", "keywords": {**_c4_tight},                                                  "error": {0: _q35}}, id="b2plyp  rhf ae: cfour",  marks=using("cfour")),
        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid},                                                       }, id="b2plyp  rhf ae: nwchem", marks=using("nwchem")),
        # ONLY DF-MP2 pytest.param({"call": "psi4",   "reference": "rhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "scf_type": "pk"},                                                               }, id="b2plyp  rhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "nwchem", "reference": "rhf",  "fcae": "fc", "keywords": {**_nwc_conv, **_nwc_grid, "mp2__freeze": 1},                                     }, id="b2plyp  rhf fc: nwchem", marks=using("nwchem")),

        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "ae", "keywords": {**_nwc_conv, **_nwc_grid, "scf__uhf": True},                                     }, id="b2plyp  uhf ae: nwchem", marks=using("nwchem")),
        # ONLY DF-MP2 pytest.param({"call": "psi4",   "reference": "uhf",  "fcae": "ae", "keywords": {**_p4_conv, **_psi_grid, "reference": "uhf", "scf_type": "pk"},                                           }, id="b2plyp  uhf ae: psi4",       marks=using("psi4")),

        pytest.param({"call": "nwchem", "reference": "uhf",  "fcae": "fc", "keywords": {**_nwc_conv, **_nwc_grid, "scf__uhf": True, "mp2__freeze": 1},                   }, id="b2plyp  uhf fc: nwchem", marks=using("nwchem")),
        # yapf: enable
    ],
)
def test_b2plyp_energy_module(inp, dertype, basis, subjects, clsd_open_pmols, request, schema_versions_v2):
    models, retver, _ = schema_versions_v2
    runner_asserter(
        *_processor(inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, "energy", "b2plyp")
    )


###################

# Notes from test_sp_mp2_rhf_full
# TODO Molpro has frozen-core on by default. For this to pass need keyword frozen_core = False
# pytest.param('molpro', 'aug-cc-pvdz', {}, marks=using("molpro")),

# Notes from test_sp_mp2_uhf_fc
# TODO Molpro needs a new keyword for unrestricted MP2 (otherwise RMP2 by default) and needs symmetry c1
# pytest.param('molpro', 'aug-cc-pvdz', {"reference": "unrestricted"}, marks=using("molpro")),

# Notes from sp_ccsd_rhf_ull
# pytest.param("qchem", "aug-cc-pvdz", {"N_FROZEN_CORE": 0}, marks=using("qchem")),
# TODO Molpro has frozen-core on by default. For this to pass need new keyword frozen_core = False
# pytest.param('molpro', 'aug-cc-pvdz', {}, marks=using("molpro")),

# Notes from sp_ccsd_rohf_full
# pytest.param("qchem", "AUG-CC-PVDZ", {"UNRESTRICTED": False}, marks=using("qchem")),
# TODO Molpro has frozen-core on by default. For this to pass need new keyword frozen_core = False
# pytest.param('molpro', 'aug-cc-pvdz', {}, marks=using("molpro")),

# def test_sp_ccsd_rohf_fc(method, keywords, nh2):
#    # differing std/semicanonical
#    ## from Cfour
#    #osccsdcorl = -0.1563275
#    #ssccsdcorl = -0.0365048
#    ##ssccsdcorl = -0.036502859383024
#
#    ## from Psi4. psi & nwchem agree to 6
#    osccsd_corl = -0.152968752464362
#    ssccsd_corl = -0.036502859383024
#
#    ## from gamess.
#    #ccsdcorl = -0.1928447371
#    #ccsdtot = -55.7775819971
#    # from Cfour
#
#    scf_tot = -55.5847372600528120
#    mp2_tot = -55.760613403041812
#    ccsd_tot = -55.7775634749542241
#
#    # gms CCSD correlation disagrees with Cfour, Psi4, and NWChem by ~2.e-4
#    # hf is in agreement across all programs
#    if method.startswith('gms'):
#        atol = 2.e-5


def _processor(
    inp, dertype, basis, subjects, clsd_open_pmols, request, models, retver, driver, method, *, scramble=None, frame=""
):
    qcprog = inp["call"]
    tnm = request.node.name
    suffix = "-fixed" if frame == "fixed" else ""
    subject = clsd_open_pmols[subjects[std_refs.index(inp["reference"])] + suffix]

    inpcopy = {k: v for k, v in inp.items() if k not in ["error", "wrong", "xptd", "marks"]}
    if inp.get("error", False) and inp["error"].get(dertype, False):
        inpcopy["error"] = inp["error"][dertype]
    if inp.get("wrong", False) and inp["wrong"].get(dertype, False):
        inpcopy["wrong"] = inp["wrong"][dertype]
    if inp.get("error", False) and inp["error"].get(basis, False):
        inpcopy["error"] = inp["error"][basis]
    if inp.get("wrong", False) and inp["wrong"].get(basis, False):
        inpcopy["wrong"] = inp["wrong"][basis]
    if inp.get("marks", False) and inp["marks"].get(dertype, False):
        request.node.add_marker(inp["marks"][dertype])

    # xptd values are plain (applies to all dertypes) or a per-dertype dict. keep falsy plain values, since
    #   `"fd": False` (accurate internal findif) is meaningful to the runner
    xptd = {}
    for k, v in inp.get("xptd", {}).items():
        if isinstance(v, dict):
            if dertype in v:
                xptd[k] = v[dertype]
        else:
            xptd[k] = v
    inpcopy["xptd"] = xptd

    inpcopy["keywords"] = {
        k: (_trans_key(qcprog, basis, k) if v == "<>" else v) for k, v in inpcopy["keywords"].items()
    }
    inpcopy["driver"] = driver
    if not any([k.lower() in _basis_keywords for k in inpcopy["keywords"]]):
        inpcopy["basis"], basis_extras = _trans_key(qcprog, basis, "qcng_basis")
        inpcopy["keywords"].update(basis_extras)
    inpcopy["keywords"].update(_system_keywords(qcprog, basis, method, driver))

    if inp["fcae"] == "fc" and (basis, subject.name) in _fc_nfrzc:
        nfrzc = _fc_nfrzc[(basis, subject.name)]
        for k, v in inpcopy["keywords"].items():
            lk = k.lower()
            if lk == "dropmo":
                inpcopy["keywords"][k] = list(range(1, nfrzc + 1))
            elif lk.endswith("__freeze") or lk.endswith("__freeze__core"):
                if isinstance(v, dict):
                    pytest.skip(f"row's per-element {k} is specific to its usual systems")
                inpcopy["keywords"][k] = nfrzc
            elif lk.endswith("__freeze__atomic"):
                pytest.skip(f"row's per-element {k} is specific to its usual systems")

    if (qcprog, basis) in _basis_spec_extras:
        inpcopy["extras"] = _basis_spec_extras[(qcprog, basis)]
    inpcopy["scf_type"] = "pk"
    inpcopy["corl_type"] = "conv"
    inpcopy["qc_module"] = "-".join(
        [
            qcprog,
            inp["keywords"].get("qc_module", inp["keywords"].get("cc_program", inp["keywords"].get("mp2__code", ""))),
        ]
    ).strip("-")

    return inpcopy, subject, method, basis, tnm, scramble, frame, models, retver
