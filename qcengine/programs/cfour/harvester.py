import re
from decimal import Decimal

import numpy as np
import qcelemental as qcel
from qcelemental.models import Molecule
from qcelemental.molparse import regex

from ..util import PreservingDict, load_hessian


def harvest_output(outtext):
    """Function to separate portions of a CFOUR output file *outtest*,
    divided by xjoda.

    """
    pass_psivar = []
    pass_coord = []
    pass_grad = []

    # for outpass in re.split(r'--invoking executable xjoda', outtext, re.MULTILINE):
    for outpass in re.split(r"JODA beginning optimization cycle", outtext, re.MULTILINE):
        psivar, c4coord, c4grad, version, error = harvest_outfile_pass(outpass)
        pass_psivar.append(psivar)
        pass_coord.append(c4coord)
        pass_grad.append(c4grad)

        # print('\n\nXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n')
        # print(outpass)
        # print(psivar, c4coord, c4grad, version, error)
        # print('\n\nxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n')

    retindx = -1 if pass_coord[-1] else -2

    #    print '    <<<  C4 PSIVAR  >>>'
    #    for item in pass_psivar[retindx]:
    #        print('       %30s %16.8f' % (item, pass_psivar[retindx][item]))
    #    print '    <<<  C4 COORD   >>>'
    #    for item in pass_coord[retindx]:
    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
    #    print '    <<<   C4 GRAD   >>>'
    #    for item in pass_grad[retindx]:
    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

    return pass_psivar[retindx], pass_coord[retindx], pass_grad[retindx], version, error


def harvest_outfile_pass(outtext):
    """Function to read CFOUR output file *outtext* and parse important
    quantum chemical information from it in

    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None
    version = ""
    error = ""

    #    TODO: BCC
    #          CI
    #          QCISD(T)
    #          other ROHF tests
    #          vcc/ecc

    NUMBER = r"(?x:" + regex.NUMBER + ")"

    # Process version
    mobj = re.search(r"^\s*" + r"Version" + r"\s+" + r"(?P<version>[\w.]+)" + r"\s*$", outtext, re.MULTILINE)
    if mobj:
        print("matched version")
        version = mobj.group("version")

    # Process NRE
    mobj = re.search(
        r"^\s+" + r"(?:Nuclear repulsion energy :)" + r"\s+" + NUMBER + r"\s+a\.u\.\s*$", outtext, re.MULTILINE
    )
    if mobj:
        print("matched nre")
        psivar["NUCLEAR REPULSION ENERGY"] = mobj.group(1)

    # Process SCF
    mobj = re.search(r"^\s+" + r"(?:E\(SCF\))" + r"\s+=\s+" + NUMBER + r"\s+a\.u\.\s*$", outtext, re.MULTILINE)
    if mobj:
        print("matched scf1")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)

    mobj = re.search(r"^\s+" + r"(?:E\(SCF\)=)" + r"\s+" + NUMBER + r"\s+" + NUMBER + r"\s*$", outtext, re.MULTILINE)
    if mobj:
        print("matched scf2")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)

    if "SCF TOTAL ENERGY" not in psivar:
        # can be too greedy and match across scf cycles
        mobj = re.search(
            # fmt: off
            r'^\s+' + r'(?:SCF has converged.)' + r'\s*$' +
            r'(?:.*?)' +
            r'^\s+' + r'(?:\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
            # fmt: on
            outtext,
            re.MULTILINE | re.DOTALL,
        )
        if mobj:
            print("matched scf3")
            psivar["SCF TOTAL ENERGY"] = mobj.group(1)

    mobj = re.search(r"^\s+" + r"(?:E\(ROHF\)=)" + r"\s+" + NUMBER + r"\s+" + NUMBER + r"\s*$", outtext, re.MULTILINE)
    if mobj:
        print("matched scf4")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)

    # Process MP2
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched mp2r")
        psivar["MP2 SAME-SPIN CORRELATION ENERGY"] = 2 * Decimal(mobj.group(1))
        psivar["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(2)
        psivar["MP2 CORRELATION ENERGY"] = 2 * Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar["MP2 TOTAL ENERGY"] = mobj.group(4)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched mp2u")
        psivar["MP2 SAME-SPIN CORRELATION ENERGY"] = Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(3)
        psivar["MP2 CORRELATION ENERGY"] = Decimal(mobj.group(1)) + Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar["MP2 TOTAL ENERGY"] = mobj.group(5)

    mobj = re.search(
        # particularly, want to avoid capture when following line present:
        #  "MP2 energies are correct only for semicanonical orbitals."
        # fmt: off
        r'Singles contribution will be calculated.' + r'\s*' +
        r'^\s+' + r'-*' + r'\s*' +
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(SINGLE\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched mp2ro")
        psivar["MP2 SAME-SPIN CORRELATION ENERGY"] = Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = mobj.group(4)
        psivar["MP2 SINGLES ENERGY"] = mobj.group(5)
        psivar["MP2 CORRELATION ENERGY"] = (
            Decimal(mobj.group(5)) + Decimal(mobj.group(2)) + Decimal(mobj.group(3)) + Decimal(mobj.group(4))
        )
        psivar["MP2 TOTAL ENERGY"] = mobj.group(7)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + r'(?P<sgl>' + NUMBER + r')' + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + r'(?P<dbl>' + NUMBER + r')' + r'\s+' +
                                                r'(?P<mp2tot>' + NUMBER + r')' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched mp2ro2")
        # psivar['MP2 SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        # psivar['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = mobj.group(3)
        psivar["MP2 SINGLES ENERGY"] = mobj.group("sgl")
        psivar["MP2 CORRELATION ENERGY"] = Decimal(mobj.group("sgl")) + Decimal(mobj.group("dbl"))
        psivar["MP2 TOTAL ENERGY"] = mobj.group("mp2tot")

    # Process MP3
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp3r")
        dmp2 = Decimal(mobj.group(1))
        dmp3 = Decimal(mobj.group(3))
        psivar["MP2 CORRELATION ENERGY"] = dmp2
        psivar["MP2 TOTAL ENERGY"] = mobj.group(2)
        psivar["MP3 CORRELATION ENERGY"] = dmp2 + dmp3
        psivar["MP3 TOTAL ENERGY"] = mobj.group(4)
        psivar["MP2.5 CORRELATION ENERGY"] = dmp2 + Decimal("0.500000000000") * dmp3
        psivar["MP2.5 TOTAL ENERGY"] = psivar["MP2.5 CORRELATION ENERGY"] + psivar["SCF TOTAL ENERGY"]
        psivar["MP3 SINGLES ENERGY"] = Decimal("0.0")
        psivar["MP3 DOUBLES ENERGY"] = dmp2 + dmp3

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp3ro")
        dmp2 = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp3 = Decimal(mobj.group(5)) + Decimal(mobj.group(7))
        psivar["MP3 CORRELATION ENERGY"] = dmp2 + dmp3
        psivar["MP3 TOTAL ENERGY"] = mobj.group(8)
        psivar["MP2.5 CORRELATION ENERGY"] = dmp2 + Decimal("0.500000000000") * dmp3
        psivar["MP2.5 TOTAL ENERGY"] = psivar["MP2.5 CORRELATION ENERGY"] + psivar["SCF TOTAL ENERGY"]
        psivar["MP3 SINGLES ENERGY"] = Decimal(mobj.group(1)) + Decimal(mobj.group(5))
        psivar["MP3 DOUBLES ENERGY"] = Decimal(mobj.group(3)) + Decimal(mobj.group(7))

    # Process MP4
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:Q-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp4r")
        dmp2 = Decimal(mobj.group(1))
        dmp3 = Decimal(mobj.group(3))
        dmp4sdq = Decimal(mobj.group(5)) + Decimal(mobj.group(7)) + Decimal(mobj.group(9))
        psivar["MP2 CORRELATION ENERGY"] = dmp2
        psivar["MP2 TOTAL ENERGY"] = mobj.group(2)
        psivar["MP3 CORRELATION ENERGY"] = dmp2 + dmp3
        psivar["MP3 TOTAL ENERGY"] = mobj.group(4)
        psivar["MP2.5 CORRELATION ENERGY"] = dmp2 + Decimal("0.500000000000") * dmp3
        psivar["MP2.5 TOTAL ENERGY"] = psivar["MP2.5 CORRELATION ENERGY"] + psivar["SCF TOTAL ENERGY"]
        psivar["MP4(SDQ) CORRELATION ENERGY"] = dmp2 + dmp3 + dmp4sdq
        psivar["MP4(SDQ) TOTAL ENERGY"] = mobj.group(10)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:L-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:NL-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp4ro")
        dmp2 = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp3 = Decimal(mobj.group(5)) + Decimal(mobj.group(7))
        dmp4sdq = Decimal(mobj.group(9)) + Decimal(mobj.group(11))
        psivar["MP2 CORRELATION ENERGY"] = dmp2
        psivar["MP2 TOTAL ENERGY"] = mobj.group(4)
        psivar["MP3 CORRELATION ENERGY"] = dmp2 + dmp3
        psivar["MP3 TOTAL ENERGY"] = mobj.group(8)
        psivar["MP2.5 CORRELATION ENERGY"] = dmp2 + Decimal("0.500000000000") * dmp3
        psivar["MP2.5 TOTAL ENERGY"] = psivar["MP2.5 CORRELATION ENERGY"] + psivar["SCF TOTAL ENERGY"]
        psivar["MP4(SDQ) CORRELATION ENERGY"] = dmp2 + dmp3 + dmp4sdq
        psivar["MP4(SDQ) TOTAL ENERGY"] = mobj.group(12)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:D-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:Q-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp4tr")
        dmp4sdq = Decimal(mobj.group(1)) + Decimal(mobj.group(3)) + Decimal(mobj.group(5))
        dmp4t = Decimal(mobj.group(7))
        psivar["MP4(SDQ) CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"] + dmp4sdq
        psivar["MP4(SDQ) TOTAL ENERGY"] = mobj.group(6)
        psivar["MP4(T) CORRECTION ENERGY"] = dmp4t
        psivar["MP4(SDTQ) CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"] + dmp4sdq + dmp4t
        psivar["MP4(SDTQ) TOTAL ENERGY"] = mobj.group(8)
        psivar["MP4 CORRELATION ENERGY"] = psivar["MP4(SDTQ) CORRELATION ENERGY"]
        psivar["MP4 TOTAL ENERGY"] = psivar["MP4(SDTQ) TOTAL ENERGY"]

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:L-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:NL-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:WT12-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched mp4tro")
        dmp4sdq = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp4t = Decimal(mobj.group(5)) + Decimal(mobj.group(7))  # TODO: WT12 with T, not SDQ?
        psivar["MP4(SDQ) CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"] + dmp4sdq
        psivar["MP4(SDQ) TOTAL ENERGY"] = mobj.group(4)
        psivar["MP4(T) CORRECTION ENERGY"] = dmp4t
        psivar["MP4(SDTQ) CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"] + dmp4sdq + dmp4t
        psivar["MP4(SDTQ) TOTAL ENERGY"] = mobj.group(8)
        psivar["MP4 CORRELATION ENERGY"] = psivar["MP4(SDTQ) CORRELATION ENERGY"]
        psivar["MP4 TOTAL ENERGY"] = psivar["MP4(SDTQ) TOTAL ENERGY"]

    # Process CC Iterations
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+DIIS\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s*(?:A miracle (?:has come|come) to pass. The CC iterations have converged.)\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched cc with full %s iterating %s" % (mobj.group("fullCC"), mobj.group("iterCC")))
        psivar["%s CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(3)
        psivar["%s TOTAL ENERGY" % (mobj.group("iterCC"))] = mobj.group(4)

        mobj3 = re.search(r"SCF reference function:  RHF", outtext)
        if mobj3:
            psivar[f"{mobj.group('iterCC')} DOUBLES ENERGY"] = mobj.group(3)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:\d+)' + r'\s+' + r'(?P<corl>' + NUMBER + r')\s+' +
                  #NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                  #NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER +  r'\s+' + NUMBER +  r'\s+' + NUMBER + r'\s*' +
                  NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER +  r'\s+' + r'(' + NUMBER +  r'\s+' + NUMBER + r')?' + r'\s*' +
        r'^\s*' +
        r'^\s*' + r'(?:\w+ iterations converged .*?)' +
        r'^\s*' +
        r'^\s*' + r'(?:Total (?P<iterCC>\w+) energy:)' + r'\s+' + r'(?P<tot>' + NUMBER + r')\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ncc cc iter")
        # looks like ncc is rhf-only
        psivar["{} CORRELATION ENERGY".format(mobj.group("iterCC"))] = mobj.group("corl")
        psivar["{} DOUBLES ENERGY".format(mobj.group("iterCC"))] = mobj.group("corl")
        psivar["{} TOTAL ENERGY".format(mobj.group("iterCC"))] = mobj.group("tot")

    # Process CC(T)
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\))' + r'\s+=\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\(T\)\))' + r'\s+=\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) vcc")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)
        psivar["CCSD TOTAL ENERGY"] = mobj.group(2)
        psivar["(T) CORRECTION ENERGY"] = Decimal(mobj.group(3)) - Decimal(mobj.group(2))
        psivar["CCSD(T) CORRELATION ENERGY"] = Decimal(mobj.group(3)) - Decimal(mobj.group(1))
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group(3)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E\(CCSD\))' + r'\s+=\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\(T\)\))' + r'\s+=\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) vcc v2")
        psivar["CCSD TOTAL ENERGY"] = mobj.group(1)
        psivar["(T) CORRECTION ENERGY"] = Decimal(mobj.group(2)) - Decimal(mobj.group(1))
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group(2)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s*' + NUMBER + r'\s+a\.u\.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:Total perturbative triples energy:)' + r'\s+' + NUMBER + r'\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) ecc")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)
        psivar["CCSD TOTAL ENERGY"] = mobj.group(2)
        psivar["(T) CORRECTION ENERGY"] = mobj.group(3)
        psivar["CCSD(T) CORRELATION ENERGY"] = Decimal(mobj.group(4)) - Decimal(mobj.group(1))
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group(4)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:HF-SCF energy)' + r'\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E4T \+ E5ST)' + r'\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) ecc v2")
        psivar["SCF TOTAL ENERGY"] = mobj.group(1)
        psivar["CCSD TOTAL ENERGY"] = mobj.group(2)
        psivar["(T) CORRECTION ENERGY"] = mobj.group(3)
        psivar["CCSD(T) CORRELATION ENERGY"] = Decimal(mobj.group(4)) - Decimal(mobj.group(1))
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group(4)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) lamb")
        psivar["CCSD TOTAL ENERGY"] = mobj.group(1)
        psivar["(T) CORRECTION ENERGY"] = Decimal(mobj.group(2)) - Decimal(mobj.group(1))
        psivar["CCSD(T) CORRELATION ENERGY"] = Decimal(mobj.group(2)) - psivar["SCF TOTAL ENERGY"]
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group(2)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?:CCSD\(T\) contribution:)\s+' + r'(?P<tcorr>' + NUMBER + ')' + r'\s*'
        r'^\s*' + r'(?:CCSD\[T\] contribution:)\s+' + r'(?P<bkttcorr>' + NUMBER + ')' + r'\s*'
        r'^\s*' + r'(?:Total CCSD\(T\) energy:)\s+' + r'(?P<ttot>' + NUMBER + ')' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) ncc")
        psivar["(T) CORRECTION ENERGY"] = mobj.group("tcorr")
        psivar["[T] CORRECTION ENERGY"] = mobj.group("bkttcorr")
        psivar["CCSD(T) TOTAL ENERGY"] = mobj.group("ttot")

    mobj = re.search(
        # fmt: off
        r'^\s*' + r'(?:CCSD\[T\] correlation energy:)\s+' + r'(?P<bkttcorr>' + NUMBER + ')' + r'\s*'
        r'^\s*' + r'(?:CCSD\(T\) correlation energy:)\s+' + r'(?P<tcorr>' + NUMBER + ')' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched ccsd(t) ncc v2")
        psivar["(T) CORRECTION ENERGY"] = mobj.group("tcorr")
        psivar["[T] CORRECTION ENERGY"] = mobj.group("bkttcorr")
        psivar["CCSD(T) TOTAL ENERGY"] = psivar["(T) CORRECTION ENERGY"] + psivar["CCSD TOTAL ENERGY"]
        psivar["CCSD(T) CORRELATION ENERGY"] = psivar["(T) CORRECTION ENERGY"] + psivar["CCSD CORRELATION ENERGY"]

    # Process DBOC
    mobj = re.search(
        # fmt: off
        r'^\s*' + r'(?:The total diagonal Born-Oppenheimer correction \(DBOC\) is:)\s+' +
        r'(?P<dboc>' + NUMBER + ')' + r'\s*a.u.\s*',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:
        print("matched dboc ecc")
        psivar["CCSD DBOC ENERGY"] = mobj.group("dboc")

    # Process SCS-CC
    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s*' + r'(?:@CCENRG-I, Correlation energies.)' + r'\s+(?:ECCAA)\s+' + NUMBER + r'\s*' +
        r'^\s+(?:ECCBB)\s+' + NUMBER + r'\s*' +
        r'^\s+(?:ECCAB)\s+' + NUMBER + r'\s*' +
        r'^\s+(?:Total)\s+' + NUMBER + r'\s*',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:  # PRINT=2 to get SCS-CC components
        print("matched scscc")
        if float(mobj.group(4)) == 0.0:
            ss = 2 * Decimal(mobj.group(3))
        else:
            ss = Decimal(mobj.group(3)) + Decimal(mobj.group(4))

        if not (
            re.search(r"executable xvcc finished", outtext)
            and re.search(r"The reference state is a ROHF wave function.", outtext)
        ):
            psivar["%s SAME-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = ss
        psivar["%s OPPOSITE-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(5)
        psivar["%s CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(6)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'Amplitude equations converged in' + r'\s*\d+\s*' + r'iterations.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'The AA contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The BB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The AB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The total correlation energy is\s+' + NUMBER + r'\s+a.u.\s*' +
        r'(?:.*?)' +
        # r'^\s+' + r'The CC iterations have converged.' + r'\s*$',
        r'^\s+' + r'(?:A miracle come to pass. )?' + r'The CC iterations have converged.' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:  # PRINT=2 to get SCS components
        print("matched scscc2")
        mobj3 = re.search(r"The reference state is a ROHF wave function.", outtext)
        mobj4 = re.search(r"executable xvcc finished", outtext)
        if mobj4:  # vcc
            psivar["%s OPPOSITE-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(5)
            if not mobj3:
                psivar[f'{mobj.group("iterCC")} SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(3)) + Decimal(
                    mobj.group(4)
                )
        else:  # ecc
            psivar[f'{mobj.group("iterCC")} SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(3)) + Decimal(
                mobj.group(4)
            )
            if not mobj3:
                psivar["%s OPPOSITE-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(5)
        psivar["%s CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(6)

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'Amplitude equations converged in' + r'\s*\d+\s*' + r'iterations.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'The AA contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The AB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The total correlation energy is\s+' + NUMBER + r'\s+a.u.\s*' +
        r'(?:.*?)' +
        # r'^\s+' + r'The CC iterations have converged.' + r'\s*$',
        r'^\s+' + r'(?:A miracle come to pass. )?' + r'The CC iterations have converged.' + r'\s*$',
        # fmt: on
        outtext,
        re.MULTILINE | re.DOTALL,
    )
    if mobj:  # PRINT=2 to get SCS components
        print("matched scscc rhf", mobj.groups())
        psivar["%s SAME-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = 2 * Decimal(mobj.group(3))
        psivar["%s OPPOSITE-SPIN CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(4)
        psivar["%s CORRELATION ENERGY" % (mobj.group("iterCC"))] = mobj.group(5)

    # Process gradient
    mobj = re.search(
        # fmt: off
        r'\s+' + r'Molecular gradient' + r'\s*' +
        r'\s+' + r'------------------' + r'\s*' +
        r'\s+' + r'\n' +
        r'(?:(?:\s+[A-Z]+\s*#\d+\s+[xyz]\s+[-+]?\d+\.\d+\s*\n)+)' +  # optional, it seems
        r'\n\n' +  # optional, it seems
        r'((?:\s+[A-Z]+\s*#\d+\s+\d?\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +
        r'\n\n' +
        r'\s+' + 'Molecular gradient norm',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched molgrad")
        atoms = []
        psivar_grad = []
        for line in mobj.group(1).splitlines():
            lline = line.split()
            atoms.append(lline[0])
            # psivar_gradient.append([Decimal(lline[-3]), Decimal(lline[-2]), Decimal(lline[-1])])
            psivar_grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])

    # Process geometry
    mobj = re.search(
        # fmt: off
        # r'\s+(?:-+)\s*' +
        # r'^\s+' + r'Z-matrix   Atomic            Coordinates (in bohr)' + r'\s*' +
        r'^\s+' + r'Symbol    Number           X              Y              Z' + r'\s*' +
        r'^\s+(?:-+)\s*' +
        r'((?:\s+[A-Z]+\s+[0-9]+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +
        r'^\s+(?:-+)\s*',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched geom")
        molxyz = "%d bohr\n\n" % len(mobj.group(1).splitlines())
        for line in mobj.group(1).splitlines():
            lline = line.split()
            molxyz += "%s %16s %16s %16s\n" % (lline[0], lline[-3], lline[-2], lline[-1])
        # Rather a dinky Molecule as no ghost, charge, or multiplicity
        psivar_coord = Molecule(
            validate=False,
            **qcel.molparse.to_schema(
                qcel.molparse.from_string(molxyz, dtype="xyz+", fix_com=True, fix_orientation=True)["qm"], dtype=2
            ),
        )

    # Process atom geometry
    mobj = re.search(r"^\s+" + r"@GETXYZ-I,     1 atoms read from ZMAT." + r"\s*$", outtext, re.MULTILINE)
    mobj2 = re.search(
        r"^([A-Z]+)#1" + r"\s+" + NUMBER + r"\s+" + NUMBER + r"\s+" + NUMBER + r"\s*$", outtext, re.MULTILINE
    )
    if mobj and mobj2:
        print("matched atom2")  # unsavory for when atom never printed except for basis file
        # Dinky Molecule
        molxyz = "1 bohr\n\n%s 0.0 0.0 0.0\n" % (mobj2.group(1))
        psivar_coord = Molecule(
            validate=False,
            **qcel.molparse.to_schema(
                qcel.molparse.from_string(molxyz, dtype="xyz+", fix_com=True, fix_orientation=True)["qm"], dtype=2
            ),
        )

    mobj = re.search(
        # fmt: off
        r'^\s+' + r'@GETXYZ-I,     1 atoms read from ZMAT.' + r'\s*' +
        r'^\s+' + r'[0-9]+\s+([A-Z]+)\s+[0-9]+\s+' + NUMBER + r'\s*',
        # fmt: on
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched atom")
        # Dinky Molecule
        molxyz = "1 bohr\n\n%s 0.0 0.0 0.0\n" % (mobj.group(1))
        psivar_coord = Molecule(
            validate=False,
            **qcel.molparse.to_schema(
                qcel.molparse.from_string(molxyz, dtype="xyz+", fix_com=True, fix_orientation=True)["qm"], dtype=2
            ),
        )

    # Process error codes
    mobj = re.search(
        r"^\s*" + r"--executable " + r"(\w+)" + r" finished with status" + r"\s+" + r"([1-9][0-9]*)",
        outtext,
        re.MULTILINE,
    )
    if mobj:
        print("matched error")
        # psivar['CFOUR ERROR CODE'] = mobj.group(2)
        if int(mobj.group(2)) != 0:
            error += "--executable {} finished with status {}".format(mobj.group(1), mobj.group(2))

    # Process CURRENT energies (TODO: needs better way)
    if "SCF TOTAL ENERGY" in psivar:
        psivar["CURRENT REFERENCE ENERGY"] = psivar["SCF TOTAL ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["SCF TOTAL ENERGY"]
        psivar["HF TOTAL ENERGY"] = psivar["SCF TOTAL ENERGY"]

    if "MP2 TOTAL ENERGY" in psivar and "MP2 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["MP2 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["MP2 TOTAL ENERGY"]

    if "MP3 TOTAL ENERGY" in psivar and "MP3 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["MP3 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["MP3 TOTAL ENERGY"]

    if "MP4 TOTAL ENERGY" in psivar and "MP4 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["MP4 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["MP4 TOTAL ENERGY"]

    #    if ('%s TOTAL ENERGY' % (mobj.group('fullCC')) in psivar) and \
    #       ('%s CORRELATION ENERGY' % (mobj.group('fullCC')) in psivar):
    #        psivar['CURRENT CORRELATION ENERGY'] = psivar['%s CORRELATION ENERGY' % (mobj.group('fullCC')]
    #        psivar['CURRENT ENERGY'] = psivar['%s TOTAL ENERGY' % (mobj.group('fullCC')]

    if "CC2 TOTAL ENERGY" in psivar and "CC2 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CC2 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CC2 TOTAL ENERGY"]

    if "CCSD TOTAL ENERGY" in psivar and "CCSD CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSD CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSD TOTAL ENERGY"]

    if "CCSD(T) TOTAL ENERGY" in psivar and "CCSD(T) CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSD(T) CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSD(T) TOTAL ENERGY"]

    if "CC3 TOTAL ENERGY" in psivar and "CC3 CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CC3 CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CC3 TOTAL ENERGY"]

    if "CCSDT TOTAL ENERGY" in psivar and "CCSDT CORRELATION ENERGY" in psivar:
        psivar["CURRENT CORRELATION ENERGY"] = psivar["CCSDT CORRELATION ENERGY"]
        psivar["CURRENT ENERGY"] = psivar["CCSDT TOTAL ENERGY"]

    return psivar, psivar_coord, psivar_grad, version, error


def harvest(p4Mol, c4out, **largs):
    """Parses all the pieces of output from Cfour: the stdout in
    *c4out* and the contents of various scratch files like GRD stored
    in their namesake keys in *largs*. Since all Cfour output uses
    its own orientation and atom ordering for the given molecule,
    a qcdb.Molecule *p4Mol*, if supplied, is used to transform the
    Cfour output back into consistency with *p4Mol*.

    """
    # Collect results from output file and subsidiary files
    outPsivar, outMol, outGrad, version, error = harvest_output(c4out)

    if largs.get("GRD"):
        grdMol, grdGrad = harvest_GRD(largs["GRD"])
    else:
        grdMol, grdGrad = None, None

    if largs.get("FCMFINAL"):
        fcmHess = load_hessian(largs["FCMFINAL"], dtype="fcmfinal")
        if np.count_nonzero(fcmHess) == 0:
            fcmHess = None
    else:
        fcmHess = None

    if largs.get("DIPOL"):
        dipolDip = harvest_DIPOL(largs["DIPOL"])
    else:
        dipolDip = None

    # Reconcile the coordinate information: several cases
    #   Case                            p4Mol   GRD      Check consistency           Apply orientation?     ReturnMol (1-19-2014)
    #   sp with mol thru cfour {}       None    None              outMol             N.C.                   outMol
    #   opt with mol thru cfour {}      None    grdMol            outMol && grdMol   N.C.                   grdMol
    #   sp with mol thru molecule {}    p4Mol   None     p4Mol && outMol             p4Mol <-- outMol       p4Mol (same as input arg)
    #   opt with mol thru molecule {}   p4Mol   grdMol   p4Mol && outMol && grdMol   p4Mol <-- grdMol       p4Mol (same as input arg)

    if outMol:
        if grdMol:
            if abs(outMol.nuclear_repulsion_energy() - grdMol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValueError(
                    """Cfour outfile (NRE: %f) inconsistent with Cfour GRD (NRE: %f)."""
                    % (outMol.nuclear_repulsion_energy(), grdMol.nuclear_repulsion_energy())
                )
        if p4Mol:
            if abs(outMol.nuclear_repulsion_energy() - p4Mol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValueError(
                    """Cfour outfile (NRE: %f) inconsistent with Psi4 input (NRE: %f)."""
                    % (outMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy())
                )
    else:
        raise ValueError("""No coordinate information extracted from Cfour output.""")

    #    print '    <<<   [1] P4-MOL   >>>'
    #    if p4Mol:
    #        p4Mol.print_out_in_bohr()
    #    print '    <<<   [2] C4-OUT-MOL   >>>'
    #    if outMol:
    #        outMol.print_out_in_bohr()
    #    print '    <<<   [3] C4-GRD-MOL   >>>'
    #    if grdMol:
    #        grdMol.print_out_in_bohr()

    # Set up array reorientation object
    if p4Mol and grdMol:
        amol, data = grdMol.align(p4Mol, atoms_map=False, mols_align=True, verbose=0)
        mill = data["mill"]

        oriCoord = mill.align_coordinates(grdMol.geometry)  # (np_out=True))
        oriGrad = mill.align_gradient(np.array(grdGrad))
        if dipolDip is None:
            oriDip = None
        else:
            oriDip = mill.align_vector(np.array(dipolDip))

        if fcmHess is None:
            oriHess = None
        else:
            oriHess = mill.align_hessian(np.array(fcmHess))

        # p4c4 = OrientMols(p4Mol, grdMol)
        # oriCoord = p4c4.transform_coordinates2(grdMol)
        # oriGrad = p4c4.transform_gradient(grdGrad)
        # oriDip = None if dipolDip is None else p4c4.transform_vector(dipolDip)

    elif p4Mol and outMol:
        # TODO watch out - haven't seen atom_map=False yet
        amol, data = outMol.align(p4Mol, atoms_map=True, mols_align=True, verbose=0)
        mill = data["mill"]

        oriCoord = mill.align_coordinates(outMol.geometry)  # (np_out=True))
        oriGrad = None
        oriHess = None  # I don't think we ever get FCMFINAL w/o GRAD
        if dipolDip is None:
            oriDip = None
        else:
            oriDip = mill.align_vector(np.array(dipolDip))
        # p4c4 = OrientMols(p4Mol, outMol)
        # oriCoord = p4c4.transform_coordinates2(outMol)
        # oriGrad = None
        # oriDip = None if dipolDip is None else p4c4.transform_vector(dipolDip)

    elif outMol:
        oriGrad = None
        oriHess = None
        oriDip = None if dipolDip is None else dipolDip

    #    print p4c4
    #    print '    <<<   [4] C4-ORI-MOL   >>>'
    #    if oriCoord is not None:
    #        for item in oriCoord:
    #            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
    #
    #    print '    <<<   [1] C4-GRD-GRAD   >>>'
    #    if grdGrad is not None:
    #        for item in grdGrad:
    #            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
    #    print '    <<<   [2] C4-ORI-GRAD   >>>'
    #    if oriGrad is not None:
    #        for item in oriGrad:
    #            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

    retMol = None if p4Mol else grdMol

    if oriDip is not None:
        outPsivar["CURRENT DIPOLE"] = oriDip
        oriDip *= qcel.constants.dipmom_au2debye
        # outPsivar["CURRENT DIPOLE X"] = oriDip[0]
        # outPsivar["CURRENT DIPOLE Y"] = oriDip[1]
        # outPsivar["CURRENT DIPOLE Z"] = oriDip[2]
        # outPsivar['CURRENT DIPOLE X'] = str(oriDip[0] * psi_dipmom_au2debye)
        # outPsivar['CURRENT DIPOLE Y'] = str(oriDip[1] * psi_dipmom_au2debye)
        # outPsivar['CURRENT DIPOLE Z'] = str(oriDip[2] * psi_dipmom_au2debye)

    if oriGrad is not None:
        retGrad = oriGrad
    elif grdGrad is not None:
        retGrad = grdGrad
    else:
        retGrad = None

    if oriHess is not None:
        retHess = oriHess
    else:
        retHess = None

    # if oriCoord is not None:
    #     retCoord = oriCoord
    # else:
    #     retCoord = None

    return outPsivar, retHess, retGrad, retMol, version, error


def harvest_GRD(grd):
    """Parses the contents *grd* of the Cfour GRD file into the gradient
    array and coordinate information. The coordinate info is converted
    into a rather dinky Molecule (no charge, multiplicity, or fragment),
    but this is these coordinates that govern the reading of molecule
    orientation by Cfour. Return qcel.models.Molecule and gradient array.

    """
    grd = grd.splitlines()
    Nat = int(grd[0].split()[0])
    molxyz = f"{Nat} bohr\n\n"

    grad = []
    for at in range(Nat):
        mline = grd[at + 1].split()
        el = "GH" if int(float(mline[0])) == 0 else qcel.periodictable.to_E(int(float(mline[0])))
        molxyz += "%s %16s %16s %16s\n" % (el, mline[-3], mline[-2], mline[-1])
        lline = grd[at + 1 + Nat].split()
        grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])
    mol = Molecule(
        validate=False,
        **qcel.molparse.to_schema(
            qcel.molparse.from_string(molxyz, dtype="xyz+", fix_com=True, fix_orientation=True)["qm"], dtype=2
        ),
    )

    return mol, grad


def harvest_DIPOL(dipol):
    """Parses the contents *dipol* of the Cfour DIPOL file into a dipol vector.

    """
    dipol = dipol.splitlines()
    lline = dipol[0].split()
    dip = [float(lline[0]), float(lline[1]), float(lline[2])]

    # return None if empty else dip
    return dip
