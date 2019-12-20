from typing import Any, Dict, Tuple

from qcengine.exceptions import InputError

# List of XC functionals known to NWChem
_xc_functionals = [
    "acm",
    "b3lyp",
    "beckehandh",
    "pbe0",
    "becke97",
    "becke97-1",
    "becke97-2",
    "becke97-3",
    "becke97-d",
    "becke98",
    "hcth",
    "hcth120",
    "hcth147",
    "hcth407",
    "becke97gga1",
    "hcth407p",
    "mpw91",
    "mpw1k",
    "xft97",
    "cft97",
    "ft97",
    "op",
    "bop",
    "pbeop",
    "xpkzb99",
    "cpkzb99",
    "xtpss03",
    "ctpss03",
    "xctpssh",
    "b1b95",
    "bb1k",
    "mpw1b95",
    "mpwb1k",
    "pw6b95",
    "pwb6k",
    "m05",
    "m05-2x",
    "vs98",
    "m06",
    "m06-hf",
    "m06-L",
    "m06-2x",
    "HFexch",
    "becke88",
    "xperdew91",
    "xpbe96",
    "gill96",
    "lyp",
    "perdew81",
    "perdew86",
    "perdew91",
    "cpbe96",
    "pw91lda",
    "slater",
    "vwn_1",
    "vwn_2",
    "vwn_3",
    "vwn_4",
    "vwn_5",
    "vwn_1_rpa",
    "xtpss03",
    "ctpss03",
    "bc95",
    "xpw6b95",
    "xpwb6k",
    "xm05",
    "xm05-2x",
    "cpw6b95",
    "cpwb6k",
    "cm05",
    "cm05-2x",
    "xvs98",
    "cvs98",
    "xm06-L",
    "xm06-hf",
    "xm06",
    "xm06-2x",
    "cm06-L",
    "cm06-hf",
    "cm06",
    "cm06-2x",
]


def muster_modelchem(method: str, derint: int, use_tce: bool) -> Tuple[str, Dict[str, Any]]:
    """Converts the QC method into NWChem keywords

    Args:
        method (str): Name of the QC method to use
        derint (str): Index of the run type
        use_tce (bool): Whether to use the Tensor Contraction Engine
    Returns:
        (str): Task command for NWChem
        (dict): Any options for NWChem
    """

    # Standardize the method name
    method = method.lower()
    opts = {}

    # Map the run type to
    runtyp = {
        0: "energy",
        1: "gradient",
        2: "hessian",
        # 'properties': 'prop',
    }[derint]

    # Write out the theory directive
    if method == "nwchem":
        mdccmd = ""

    elif method in ["scf", "hf"]:
        mdccmd = f"task scf {runtyp}\n\n"

    elif method == "mp2":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__mp2"] = True
        else:
            mdccmd = f"task mp2 {runtyp}\n\n"

    elif method == "mp3":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__mp3"] = True

    elif method == "mp4":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__mp4"] = True

    elif method == "ccd":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccd"] = True

    elif method == "ccsd":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccsd"] = True
        else:
            mdccmd = f"task ccsd {runtyp}\n\n"

    elif method == "ccsdt":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccsdt"] = True
        else:
            mdccmd = f"task ccsdt {runtyp}\n\n"

    elif method == "ccsd(t)":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__ccsd(t)"] = True
        else:
            mdccmd = f"task ccsd(t) {runtyp}\n\n"

    elif method == "tddft":
        mdccmd = f"task tddft {runtyp}\n\n"

    elif method in ["sodft", "direct_mp2", "rimp2", "mcscf", "selci", "md", "pspw", "band"]:
        raise InputError(f'Method "{method}" not yet supported by QCEngine')

    elif method == "tce":
        raise InputError(
            f"Do not specify TCE as a method. Instead specify the desired method " f'as a keyword and "qc_module=True".'
        )

    elif method in _xc_functionals:
        opts["dft__xc"] = method
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__"] = "dft"
        elseimport qcelemental as qcel
import qcengine as qcng
sdimerF='/home/adrian/projects/qca/benzeneDimer/Sdimer-004.xyz'
pddimerF='/home/adrian/projects/qca/benzeneDimer/PDdimer-051.xyz'
SDIMER_MOL = qcel.models.Molecule.from_file(filename=sdimerF)
PDDIMER_MOL = qcel.models.Molecule.from_file(filename=pddimerF)

KW = {'dft__xc': 'xc vwn_1_rpa 0.19 lyp 0.81 HFexch 0.20  slater 0.80 becke88 nonlocal 0.72',
      'dft__direct': True,
      # 'dft__disp__vdw': 3,
      'dft__grid': 'fine',
      'dft__tolerances__tight': True,
      'dft__tolerances__acccoul': 10,
      'dft__tolerances__radius': 30,
      'dft__convergence__density': 1e-6}

blyp = {'dft__xc': 'becke88 lyp',  # blyp
        'dft__direct': True,
        'dft__grid': 'fine',
        'dft__tolerances__tight': True,
        'dft__tolerances__radius': 30,
        'dft__tolerances__acccoul': 10,
        'dft__convergence__density': 1e-6}

b3lyp = {'dft__xc': 'b3lyp',
         'dft__direct': True,
         # 'dft__disp__vdw': 3,
         'dft__grid': 'fine',
         'dft__tolerances__tight': True,
         'dft__tolerances__acccoul': 10,
         'dft__tolerances__radius': 30,
         'dft__convergence__density': 1e-6}


inp2 = qcel.models.AtomicInput(
    molecule=PDDIMER_MOL,
    driver="energy",
    model={"method": "dft", "basis": 'sto-3g'},
    keywords=blyp
)
inp3 = qcel.models.AtomicInput(
    molecule=PDDIMER_MOL,
    driver="energy",
    model={"method": "dft", "basis": 'sto-3g'},
    keywords=b3lyp
)

# r1 = qcng.compute(inp1, 'nwchem')
r2 = qcng.compute(inp2, 'nwchem')
r3 = qcng.compute(inp3, 'nwchem')
            mdccmd = f"task dft {runtyp}\n\n"

    elif method == "dft":
        if use_tce:
            mdccmd = f"task tce {runtyp}\n\n"
            opts["tce__"] = "dft"
        else:
            mdccmd = f"task dft {runtyp}\n\n"

    else:
        raise InputError(f"Method not recognized: {method}")

    return mdccmd, opts
