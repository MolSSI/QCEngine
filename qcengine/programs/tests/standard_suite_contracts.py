# This file consolidates the expectations of what properties/QCVariables
#  a method should produce and what in practice a QC program does produce.

__all__ = [
    "contractual_current",
    "contractual_hf",
    "contractual_mp2",
    "contractual_mp2p5",
    "contractual_mp3",
    "contractual_mp4_prsdq_pr",
    "contractual_mp4",
    "contractual_zapt2",
    "contractual_cisd",
    "contractual_qcisd",
    "contractual_qcisd_prt_pr",
    "contractual_fci",
    "contractual_remp2",
    "contractual_lccd",
    "contractual_lccsd",
    "contractual_cepa_pr1_pr",
    "contractual_cepa_pr3_pr",
    "contractual_acpf",
    "contractual_aqcc",
    "contractual_ccd",
    "contractual_bccd",
    "contractual_cc2",
    "contractual_ccsd",
    "contractual_ccsdpt_prccsd_pr",
    "contractual_ccsd_prt_pr",
    "contractual_accsd_prt_pr",
    "contractual_bccd_prt_pr",
    "contractual_cc3",
    "contractual_ccsdt",
    "contractual_ccsdt1a",
    "contractual_ccsdt1b",
    "contractual_ccsdt2",
    "contractual_ccsdt3",
    "contractual_ccsdt_prq_pr",
    "contractual_ccsdtq",
    "contractual_omp2",
    "contractual_omp2p5",
    "contractual_omp3",
    "contractual_oremp2",
    "contractual_olccd",
    "contractual_occd",
    "contractual_occd_prt_pr",
    "contractual_aoccd_prt_pr",
    "contractual_dft_current",
    "contractual_dhdft_current",
    "query_qcvar",
    "query_has_qcvar",
]


from typing import Any, Tuple

from qcengine.programs.qcvar_identities_resources import qcvars_to_atomicproperties


def query_qcvar(obj: Any, pv: str) -> Any:
    """Uniform interface to value of variable `pv` or its QCSchema alias in QCVariable store `obj`."""

    try:
        # psi4.core (module -- P::e.globals)
        # psi4.Wavefunction
        # qcdb (module)
        vval = obj.variable(pv)
    except AttributeError:
        try:
            # qcdb jobrec["qcvars"] qcel.Datum
            vval = obj[pv].data
        except (AttributeError, KeyError):
            # qcel.AtomicResult.extras["qcvars"]
            vval = obj.get(pv)
            if vval is None:
                # qcel.AtomicResult.properties dict
                vval = obj.get(qcvars_to_atomicproperties[pv])
        except TypeError:
            # qcel.AtomicResult.properties object
            vval = getattr(obj, qcvars_to_atomicproperties[pv])

    return vval


def query_has_qcvar(obj: Any, pv: str) -> bool:
    """Uniform interface to whether variable `pv` or its QCSchema alias in QCVariable store `obj`."""

    try:
        bval = obj.has_variable(pv)
    except AttributeError:
        bval = pv in obj

        if not bval:
            try:
                bval = qcvars_to_atomicproperties[pv] in obj
            except KeyError:
                bval = False

    return bval


_contractual_docstring = """
    Parameters
    ----------
    qc_module
        The program or subprogram running the job (e.g., "cfour" or "cfour-ecc").
    driver
        {"energy", "gradient", "hessian"}
        The derivative level that should be expected.
    reference
        {"rhf", "uhf", "rohf"}
        The SCF reference since programs often output differently based on it.
    method
        The target AtomicInput.model.method since "free" methods may not always be
        output (e.g., MP2 available when target is MP2 but not when target is CCSD).
    corl_type
        {"conv", "df", "cd"}
        The algorithm for the target method since programs often output differently
        based on it.
    fcae
        {"ae", "fc"}
        The all-electron vs. frozen-orbital aspect.

    Returns
    -------
    (rpv, pv, expected)
        Of all the QCVariables `pv` that should be available, returns tuple of
        whether `expected` and what key `rpv` in the reference `pv` should match.

"""


def contractual_current(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Given the target method, returns the CURRENT QCVariables that should be produced.

    {_contractual_docstring}
    """
    contractual_qcvars = [
        ("HF TOTAL ENERGY", "SCF TOTAL ENERGY"),
        ("HF TOTAL ENERGY", "CURRENT REFERENCE ENERGY"),
        (f"{method.upper()} CORRELATION ENERGY", "CURRENT CORRELATION ENERGY"),
        (f"{method.upper()} TOTAL ENERGY", "CURRENT ENERGY"),
    ]
    if driver == "gradient":
        contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
    elif driver == "hessian":
        # contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
        contractual_qcvars.append((f"{method.upper()} TOTAL HESSIAN", "CURRENT HESSIAN"))

    for rpv, pv in contractual_qcvars:
        expected = True
        if method == "hf" and rpv == f"{method.upper()} CORRELATION ENERGY":
            expected = False

        yield (rpv, pv, expected)


def contractual_hf(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal HF should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """

    contractual_qcvars = [
        ("HF TOTAL ENERGY", "HF TOTAL ENERGY"),
        ("HF TOTAL ENERGY", "SCF TOTAL ENERGY"),
    ]
    if driver == "gradient" and method == "hf":
        contractual_qcvars.append(("HF TOTAL GRADIENT", "HF TOTAL GRADIENT"))
        # contractual_qcvars.append(("HF TOTAL GRADIENT", "SCF TOTAL GRADIENT"))
    elif driver == "hessian" and method == "hf":
        # contractual_qcvars.append(("HF TOTAL GRADIENT", "HF TOTAL GRADIENT"))
        contractual_qcvars.append(("HF TOTAL HESSIAN", "HF TOTAL HESSIAN"))
        # contractual_qcvars.append(("HF TOTAL GRADIENT", "SCF TOTAL GRADIENT"))

    for rpv, pv in contractual_qcvars:
        yield (rpv, pv, True)


def contractual_mp2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP2 CORRELATION ENERGY",
        "MP2 TOTAL ENERGY",
        "MP2 SAME-SPIN CORRELATION ENERGY",
        "MP2 SINGLES ENERGY",
        "MP2 DOUBLES ENERGY",
        "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "mp2":
        contractual_qcvars.append("MP2 TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp2":
        # contractual_qcvars.append("MP2 TOTAL GRADIENT")
        contractual_qcvars.append("MP2 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module == "cfour" and reference == "rohf" and method == "mp2" and driver == "hessian")
                    or (
                        qc_module in ["gamess-serial", "gamess-ddi"]
                        and reference in ["uhf", "rohf"]
                        and method == "mp2"
                    )
                    or (
                        qc_module in ["gamess-serial", "gamess-ims"]
                        and reference == "rhf"
                        and method == "mp2"
                        and driver in ["gradient", "hessian"]
                    )
                    or (
                        qc_module == "gamess"
                        and reference in ["rhf"]
                        and method in ["lccd", "ccd", "ccsd", "ccsd+t(ccsd)", "ccsd(t)"]
                    )
                    or (qc_module == "nwchem-tce" and method in ["mp2", "mp3", "mp4"])
                    or (
                        qc_module == "nwchem-cc"
                        and reference in ["rhf"]
                        and method in ["ccsd", "ccsd+t(ccsd)", "ccsd(t)"]
                    )
                    or (qc_module == "nwchem-directmp2" and reference == "rhf" and method == "mp2")
                    or (
                        qc_module == "psi4-occ"
                        and reference == "rhf"
                        and corl_type in ["df", "cd"]
                        and method
                        in [
                            "mp2",
                            "mp2.5",
                            "mp3",
                            "remp2",
                            "lccd",
                            "ccd",
                            "ccsd",
                            "ccsd(t)",
                            "a-ccsd(t)",
                            "omp2",
                            "omp2.5",
                            "omp3",
                            "oremp2",
                            "olccd",
                            "occd",
                            "occd(t)",
                            "a-occd(t)",
                        ]
                    )
                    or (
                        qc_module == "psi4-mrcc"
                        and reference in ["rhf", "uhf"]
                        and method
                        in ["ccsd", "ccsd(t)", "a-ccsd(t)", "cc3", "ccsdt-1a", "ccsdt-1b", "ccsdt-3", "ccsdt"]
                    )
                )
                and pv in ["MP2 SAME-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (
                    (qc_module == "psi4-detci" and method in ["mp2", "mp3", "mp4"])
                    or (
                        qc_module == "qchem" and method == "mp2"
                    )  # for structured -- can probably get these from parsing
                )
                and pv
                in [
                    "MP2 SAME-SPIN CORRELATION ENERGY",
                    "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
                    "MP2 SINGLES ENERGY",
                    "MP2 DOUBLES ENERGY",
                ]
            )
            or (
                (
                    (
                        qc_module == "psi4-occ"
                        and reference == "rohf"
                        and corl_type == "conv"
                        and method in ["omp2", "omp2.5", "omp3", "oremp2", "olccd"]
                    )
                )
                and pv
                in [
                    "MP2 CORRELATION ENERGY",
                    "MP2 TOTAL ENERGY",
                    "MP2 SINGLES ENERGY",
                ]
            )
            or (
                (
                    (qc_module == "psi4-ccenergy" and reference == "rohf" and method == "ccsd" and sdsc == "sd")
                    # next two are an evasion as possibly collectable
                    or (qc_module == "cfour-vcc" and reference in ["rhf", "uhf"] and method == "bccd")
                    or (qc_module == "cfour-vcc" and reference in ["rhf"] and method == "bccd(t)")
                    or (
                        qc_module == "nwchem-tce"
                        and method
                        in ["qcisd", "lccd", "lccsd", "ccd", "cc2", "ccsd", "ccsd+t(ccsd)", "ccsd(t)", "ccsdt"]
                    )
                    or (qc_module == "gamess" and reference == "rohf" and method == "ccsd")
                    or (
                        qc_module.startswith("cfour")
                        and reference == "rohf"
                        and method in ["lccsd", "ccd", "ccsd", "ccsd(t)", "ccsdt"]
                        and sdsc == "sd"
                    )  # this is a cop out as c4 perfectly able to produce good rohf mp2 but not with same orbitals as ref definition on ccsd
                    or (
                        qc_module == "psi4-mrcc"
                        and reference in ["rohf"]
                        and method in ["ccsd", "ccsd(t)", "a-ccsd(t)", "ccsdt"]
                    )
                )
                and pv
                in [
                    "MP2 CORRELATION ENERGY",
                    "MP2 TOTAL ENERGY",
                    "MP2 SAME-SPIN CORRELATION ENERGY",
                    "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
                    "MP2 SINGLES ENERGY",
                    "MP2 DOUBLES ENERGY",
                ]
            )
        ):
            expected = False

        yield (pv, pv, expected)


#        # TODO check CUSTOM SCS-MP2 _absent_


def contractual_mp2p5(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal MP2.5 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP2.5 CORRELATION ENERGY",
        "MP2.5 TOTAL ENERGY",
        "MP2.5 SAME-SPIN CORRELATION ENERGY",
        "MP2.5 SINGLES ENERGY",
        "MP2.5 DOUBLES ENERGY",
        "MP2.5 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "mp2.5":
        contractual_qcvars.append("MP2.5 TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp2.5":
        # contractual_qcvars.append("MP2.5 TOTAL GRADIENT")
        contractual_qcvars.append("MP2.5 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (qc_module.startswith("cfour") and method in ["mp3", "mp4(sdq)", "mp4"])
                or (
                    qc_module == "psi4-occ"
                    and reference == "rhf"
                    and corl_type in ["df", "cd"]
                    and method in ["mp2.5", "mp3", "omp2.5"]
                )
                or (qc_module.startswith("nwchem") and method in ["mp3", "mp4"])
            )
            and pv in ["MP2.5 SAME-SPIN CORRELATION ENERGY", "MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"]
        ) or (
            (
                (qc_module == "psi4-detci" and method in ["mp3", "mp4"])
                or (
                    qc_module == "psi4-occ"
                    and reference == "rohf"
                    and corl_type in ["conv", "df", "cd"]
                    and method in ["omp2.5", "omp3"]
                )
            )
            # Note SS/OS might be obtainable but no reference to verify
            and pv
            in [
                "MP2.5 CORRELATION ENERGY",
                "MP2.5 TOTAL ENERGY",
                "MP2.5 SAME-SPIN CORRELATION ENERGY",
                "MP2.5 OPPOSITE-SPIN CORRELATION ENERGY",
                "MP2.5 SINGLES ENERGY",
                "MP2.5 DOUBLES ENERGY",
            ]
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_mp3(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP3 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP3 CORRELATION ENERGY",
        "MP3 TOTAL ENERGY",
        "MP3 SAME-SPIN CORRELATION ENERGY",
        "MP3 SINGLES ENERGY",
        "MP3 DOUBLES ENERGY",
        "MP3 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "mp3":
        contractual_qcvars.append("MP3 TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp3":
        # contractual_qcvars.append("MP3 TOTAL GRADIENT")
        contractual_qcvars.append("MP3 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module.startswith("cfour") and method in ["mp3", "mp4(sdq)", "mp4"])
                    or (
                        qc_module == "psi4-occ"
                        and reference == "rhf"
                        and corl_type in ["df", "cd"]
                        and method in ["mp2.5", "mp3", "omp2.5", "omp3"]
                    )
                    or (qc_module.startswith("nwchem") and method in ["mp3", "mp4"])
                )
                and pv in ["MP3 SAME-SPIN CORRELATION ENERGY", "MP3 OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                ((qc_module == "psi4-detci" and method in ["mp3", "mp4"]))
                and pv
                in [
                    "MP3 SAME-SPIN CORRELATION ENERGY",
                    "MP3 OPPOSITE-SPIN CORRELATION ENERGY",
                    "MP3 SINGLES ENERGY",
                    "MP3 DOUBLES ENERGY",
                ]
            )
            or (
                (
                    (
                        qc_module == "psi4-occ"
                        and reference == "rohf"
                        and corl_type in ["conv", "df", "cd"]
                        and method in ["omp2.5", "omp3"]
                    )
                )
                # Note SS/OS might be obtainable but no reference to verify
                and pv
                in [
                    "MP3 CORRELATION ENERGY",
                    "MP3 TOTAL ENERGY",
                    "MP3 SAME-SPIN CORRELATION ENERGY",
                    "MP3 SINGLES ENERGY",
                    "MP3 DOUBLES ENERGY",
                    "MP3 OPPOSITE-SPIN CORRELATION ENERGY",
                ]
            )
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_mp4_prsdq_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP4(SDQ) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP4(SDQ) CORRELATION ENERGY",
        "MP4(SDQ) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "mp4(sdq)":
        contractual_qcvars.append("MP4(SDQ) TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp4(sdq)":
        # contractual_qcvars.append("MP4(SDQ) TOTAL GRADIENT")
        contractual_qcvars.append("MP4(SDQ) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (qc_module.startswith("psi4-detci") and method == "mp4")
            or (qc_module.startswith("nwchem") and method == "mp4")
        ) and pv in ["MP4(SDQ) TOTAL ENERGY", "MP4(SDQ) CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_mp4(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP4 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "MP4(T) CORRECTION ENERGY",
        "MP4 CORRELATION ENERGY",
        "MP4 CORRECTION ENERGY",
        "MP4 TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "mp4":
        contractual_qcvars.append("MP4 TOTAL GRADIENT")
    elif driver == "hessian" and method == "mp4":
        # contractual_qcvars.append("MP4 TOTAL GRADIENT")
        contractual_qcvars.append("MP4 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (qc_module.startswith("psi4-detci") and method == "mp4")
            or (qc_module.startswith("nwchem") and method == "mp4")
        ) and pv in ["MP4(T) CORRECTION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_zapt2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal ZAPT2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "ZAPT2 CORRELATION ENERGY",
        "ZAPT2 TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "zapt2":
        contractual_qcvars.append("ZAPT2 TOTAL GRADIENT")
    elif driver == "hessian" and method == "zapt2":
        # contractual_qcvars.append("ZAPT2 TOTAL GRADIENT")
        contractual_qcvars.append("ZAPT2 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_cisd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CISD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CISD CORRELATION ENERGY",
        "CISD TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "cisd":
        contractual_qcvars.append("CISD TOTAL GRADIENT")
    elif driver == "hessian" and method == "cisd":
        # contractual_qcvars.append("CISD TOTAL GRADIENT")
        contractual_qcvars.append("CISD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_qcisd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal QCISD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "QCISD CORRELATION ENERGY",
        "QCISD TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "qcisd":
        contractual_qcvars.append("QCISD TOTAL GRADIENT")
    elif driver == "hessian" and method == "qcisd":
        # contractual_qcvars.append("QCISD TOTAL GRADIENT")
        contractual_qcvars.append("QCISD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_qcisd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal QCISD(T) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "QCISD(T) CORRECTION ENERGY",
        "QCISD(T) CORRELATION ENERGY",
        "QCISD(T) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "qcisd(t)":
        contractual_qcvars.append("QCISD(T) TOTAL GRADIENT")
    elif driver == "hessian" and method == "qcisd(t)":
        # contractual_qcvars.append("QCISD(T) TOTAL GRADIENT")
        contractual_qcvars.append("QCISD(T) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_fci(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal FCI should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "FCI CORRELATION ENERGY",
        "FCI TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "fci":
        contractual_qcvars.append("FCI TOTAL GRADIENT")
    elif driver == "hessian" and method == "fci":
        # contractual_qcvars.append("FCI TOTAL GRADIENT")
        contractual_qcvars.append("FCI TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_remp2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal REMP2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "REMP2 CORRELATION ENERGY",
        "REMP2 TOTAL ENERGY",
        "REMP2 SAME-SPIN CORRELATION ENERGY",
        "REMP2 SINGLES ENERGY",
        "REMP2 DOUBLES ENERGY",
        "REMP2 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "remp2":
        contractual_qcvars.append("REMP2 TOTAL GRADIENT")
    elif driver == "hessian" and method == "remp2":
        # contractual_qcvars.append("REMP2 TOTAL GRADIENT")
        contractual_qcvars.append("REMP2 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            ((qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"] and method == "remp2"))
            and pv in ["REMP2 SAME-SPIN CORRELATION ENERGY", "REMP2 OPPOSITE-SPIN CORRELATION ENERGY"]
        ) or (
            (
                qc_module == "psi4-occ"
                and reference in ["rhf", "uhf", "rohf"]
                and corl_type in ["conv", "df", "cd"]
                and method == "oremp2"
            )
            and pv
            in [
                "REMP2 CORRELATION ENERGY",
                "REMP2 TOTAL ENERGY",
                "REMP2 SAME-SPIN CORRELATION ENERGY",
                "REMP2 SINGLES ENERGY",
                "REMP2 DOUBLES ENERGY",
                "REMP2 OPPOSITE-SPIN CORRELATION ENERGY",
            ]
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_lccd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal LCCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "LCCD CORRELATION ENERGY",
        "LCCD TOTAL ENERGY",
        "LCCD SAME-SPIN CORRELATION ENERGY",
        "LCCD SINGLES ENERGY",
        "LCCD DOUBLES ENERGY",
        "LCCD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "lccd":
        contractual_qcvars.append("LCCD TOTAL GRADIENT")
    elif driver == "hessian" and method == "lccd":
        # contractual_qcvars.append("LCCD TOTAL GRADIENT")
        contractual_qcvars.append("LCCD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"] and method == "lccd")
                    or (qc_module == "cfour-ncc" and reference in ["rhf"] and method == "lccd")
                    or (qc_module == "nwchem-tce" and reference in ["rhf", "uhf"] and method == "lccd")
                    or (qc_module == "gamess" and reference in ["rhf"] and method == "lccd")
                )
                and pv in ["LCCD SAME-SPIN CORRELATION ENERGY", "LCCD OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (qc_module == "nwchem-tce" and reference in ["rohf"] and method in ["lccd"])
                and pv
                in [
                    "LCCD SAME-SPIN CORRELATION ENERGY",
                    "LCCD OPPOSITE-SPIN CORRELATION ENERGY",
                    "LCCD SINGLES ENERGY",
                    "LCCD DOUBLES ENERGY",
                ]
            )
            or (
                (
                    qc_module == "psi4-occ"
                    and reference in ["rhf", "uhf", "rohf"]
                    and corl_type in ["conv", "df", "cd"]
                    and method == "olccd"
                )
                and pv
                in [
                    "LCCD CORRELATION ENERGY",
                    "LCCD TOTAL ENERGY",
                    "LCCD SAME-SPIN CORRELATION ENERGY",
                    "LCCD SINGLES ENERGY",
                    "LCCD DOUBLES ENERGY",
                    "LCCD OPPOSITE-SPIN CORRELATION ENERGY",
                ]
            )
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_lccsd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal LCCSD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "LCCSD CORRELATION ENERGY",
        "LCCSD TOTAL ENERGY",
        "LCCSD SAME-SPIN CORRELATION ENERGY",
        "LCCSD SINGLES ENERGY",
        "LCCSD DOUBLES ENERGY",
        "LCCSD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "lccsd":
        contractual_qcvars.append("LCCSD TOTAL GRADIENT")
    elif driver == "hessian" and method == "lccsd":
        # contractual_qcvars.append("LCCSD TOTAL GRADIENT")
        contractual_qcvars.append("LCCSD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module == "cfour-ncc" and reference in ["rhf"] and method == "lccsd")
                    or (qc_module == "nwchem-tce" and reference in ["rhf", "uhf"] and method == "lccsd")
                )
                and pv in ["LCCSD SAME-SPIN CORRELATION ENERGY", "LCCSD OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (qc_module == "cfour-vcc" and reference in ["rohf"] and method in ["lccsd"])
                and pv
                in [
                    "LCCSD SAME-SPIN CORRELATION ENERGY",
                    "LCCSD SINGLES ENERGY",
                    "LCCSD DOUBLES ENERGY",
                ]
            )
            or (
                (qc_module == "nwchem-tce" and reference in ["rohf"] and method in ["lccsd"])
                and pv
                in [
                    "LCCSD SAME-SPIN CORRELATION ENERGY",
                    "LCCSD OPPOSITE-SPIN CORRELATION ENERGY",
                    "LCCSD SINGLES ENERGY",
                    "LCCSD DOUBLES ENERGY",
                ]
            )
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_cepa_pr1_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CEPA(1) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CEPA(1) CORRELATION ENERGY",
        "CEPA(1) TOTAL ENERGY",
        "CEPA(1) SAME-SPIN CORRELATION ENERGY",
        "CEPA(1) SINGLES ENERGY",
        "CEPA(1) DOUBLES ENERGY",
        "CEPA(1) OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "cepa(1)":
        contractual_qcvars.append("CEPA(1) TOTAL GRADIENT")
    elif driver == "hessian" and method == "cepa(1)":
        # contractual_qcvars.append("CEPA(1) TOTAL GRADIENT")
        contractual_qcvars.append("CEPA(1) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (False) and pv in ["CEPA(1) SAME-SPIN CORRELATION ENERGY", "CEPA(1) OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_cepa_pr3_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CEPA(3) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CEPA(3) CORRELATION ENERGY",
        "CEPA(3) TOTAL ENERGY",
        "CEPA(3) SAME-SPIN CORRELATION ENERGY",
        "CEPA(3) SINGLES ENERGY",
        "CEPA(3) DOUBLES ENERGY",
        "CEPA(3) OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "cepa(3)":
        contractual_qcvars.append("CEPA(3) TOTAL GRADIENT")
    elif driver == "hessian" and method == "cepa(3)":
        # contractual_qcvars.append("CEPA(3) TOTAL GRADIENT")
        contractual_qcvars.append("CEPA(3) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (False) and pv in ["CEPA(3) SAME-SPIN CORRELATION ENERGY", "CEPA(3) OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_acpf(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal ACPF should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "ACPF CORRELATION ENERGY",
        "ACPF TOTAL ENERGY",
        "ACPF SAME-SPIN CORRELATION ENERGY",
        "ACPF SINGLES ENERGY",
        "ACPF DOUBLES ENERGY",
        "ACPF OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "acpf":
        contractual_qcvars.append("ACPF TOTAL GRADIENT")
    elif driver == "hessian" and method == "acpf":
        # contractual_qcvars.append("ACPF TOTAL GRADIENT")
        contractual_qcvars.append("ACPF TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (False) and pv in ["ACPF SAME-SPIN CORRELATION ENERGY", "ACPF OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_aqcc(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal AQCC should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "AQCC CORRELATION ENERGY",
        "AQCC TOTAL ENERGY",
        "AQCC SAME-SPIN CORRELATION ENERGY",
        "AQCC SINGLES ENERGY",
        "AQCC DOUBLES ENERGY",
        "AQCC OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "aqcc":
        contractual_qcvars.append("AQCC TOTAL GRADIENT")
    elif driver == "hessian" and method == "aqcc":
        # contractual_qcvars.append("AQCC TOTAL GRADIENT")
        contractual_qcvars.append("AQCC TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (False) and pv in ["AQCC SAME-SPIN CORRELATION ENERGY", "AQCC OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_ccd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCD CORRELATION ENERGY",
        "CCD TOTAL ENERGY",
        "CCD SAME-SPIN CORRELATION ENERGY",
        "CCD SINGLES ENERGY",
        "CCD DOUBLES ENERGY",
        "CCD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "ccd":
        contractual_qcvars.append("CCD TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccd":
        # contractual_qcvars.append("CCD TOTAL GRADIENT")
        contractual_qcvars.append("CCD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (qc_module == "cfour-ecc" and reference in ["rhf"] and method == "ccd")
                or (qc_module == "cfour-ncc" and reference in ["rhf"] and method == "ccd")
                or (qc_module == "nwchem-tce" and reference in ["rhf", "uhf"] and method == "ccd")
                or (qc_module == "gamess" and reference in ["rhf"] and method == "ccd")
                or (qc_module == "psi4-occ" and reference in ["rhf", "uhf"] and method == "ccd")
            )
            and pv in ["CCD SAME-SPIN CORRELATION ENERGY", "CCD OPPOSITE-SPIN CORRELATION ENERGY"]
        ) or (
            (
                (qc_module == "cfour-vcc" and reference in ["rohf"] and method in ["ccd"])
                or (qc_module == "nwchem-tce" and reference in ["rohf"] and method in ["ccd"])
            )
            and pv
            in [
                "CCD SAME-SPIN CORRELATION ENERGY",
                "CCD SINGLES ENERGY",
                "CCD DOUBLES ENERGY",
                "CCD OPPOSITE-SPIN CORRELATION ENERGY",
            ]
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_bccd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal BCCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "BCCD CORRELATION ENERGY",
        "BCCD TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "bccd":
        contractual_qcvars.append("BCCD TOTAL GRADIENT")
    elif driver == "hessian" and method == "bccd":
        # contractual_qcvars.append("BCCD TOTAL GRADIENT")
        contractual_qcvars.append("BCCD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_cc2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CC2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CC2 CORRELATION ENERGY",
        "CC2 TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "cc2":
        contractual_qcvars.append("CC2 TOTAL GRADIENT")
    elif driver == "hessian" and method == "cc2":
        # contractual_qcvars.append("CC2 TOTAL GRADIENT")
        contractual_qcvars.append("CC2 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSD CORRELATION ENERGY",
        "CCSD TOTAL ENERGY",
        "CCSD SAME-SPIN CORRELATION ENERGY",
        "CCSD SINGLES ENERGY",
        "CCSD DOUBLES ENERGY",
        "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "ccsd":
        contractual_qcvars.append("CCSD TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsd":
        # contractual_qcvars.append("CCSD TOTAL GRADIENT")
        contractual_qcvars.append("CCSD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module == "gamess" and reference == "rhf" and method in ["ccsd", "ccsd+t(ccsd)", "ccsd(t)"])
                    or (
                        qc_module == "nwchem-tce"
                        and reference in ["rhf", "uhf"]
                        and method in ["ccsd", "ccsd+t(ccsd)", "ccsd(t)"]
                    )
                    or (
                        qc_module in ["cfour-ncc", "cfour-ecc"]
                        and reference in ["rhf"]
                        and method in ["ccsd", "ccsd+t(ccsd)", "ccsd(t)", "a-ccsd(t)"]
                    )
                    or (
                        qc_module == "psi4-occ"
                        and reference in ["rhf", "uhf"]
                        and corl_type in ["df", "cd"]
                        and method in ["ccsd", "ccsd(t)", "a-ccsd(t)"]
                    )
                    or (qc_module in ["cfour-vcc"] and reference in ["rhf", "uhf"] and method in ["ccsd+t(ccsd)"])
                    or (
                        qc_module == "psi4-mrcc"
                        and reference in ["rhf", "uhf"]
                        and method in ["ccsd", "ccsd(t)", "a-ccsd(t)"]
                    )
                )
                and pv in ["CCSD SAME-SPIN CORRELATION ENERGY", "CCSD OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (qc_module == "cfour-vcc" and reference in ["rohf"] and method in ["ccsd", "ccsd(t)"])
                and pv
                in [
                    "CCSD SAME-SPIN CORRELATION ENERGY",
                    "CCSD SINGLES ENERGY",
                    "CCSD DOUBLES ENERGY",
                ]
            )
            or (
                (qc_module == "cfour-ecc" and reference in ["rohf"] and method in ["ccsd", "ccsd(t)"])
                and pv
                in [
                    "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
                    "CCSD SINGLES ENERGY",
                    "CCSD DOUBLES ENERGY",
                ]
            )
            or (
                (
                    (qc_module == "gamess" and reference in ["rohf"] and method == "ccsd")
                    or (qc_module == "nwchem-tce" and reference in ["rohf"] and method in ["ccsd", "ccsd(t)"])
                    or (
                        qc_module == "psi4-mrcc"
                        and reference in ["rohf"]
                        and method in ["ccsd", "ccsd(t)", "a-ccsd(t)"]
                    )
                )
                and pv
                in [
                    "CCSD SAME-SPIN CORRELATION ENERGY",
                    "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
                    "CCSD SINGLES ENERGY",
                    "CCSD DOUBLES ENERGY",
                ]
            )
            or (
                ((qc_module == "cfour-vcc" and method in ["bccd", "bccd(t)"]))
                and pv
                in [
                    "CCSD CORRELATION ENERGY",
                    "CCSD TOTAL ENERGY",
                    "CCSD SAME-SPIN CORRELATION ENERGY",
                    "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
                    "CCSD SINGLES ENERGY",
                    "CCSD DOUBLES ENERGY",
                ]
            )
        ):
            expected = False

        yield (pv, pv, expected)

    # TODO check CUSTOM SCS-CCSD _absent_


def contractual_ccsdpt_prccsd_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSD+T(CCSD) (aka CCSD[T]) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "T(CCSD) CORRECTION ENERGY",
        "CCSD+T(CCSD) CORRELATION ENERGY",
        "CCSD+T(CCSD) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsd+t(ccsd)":
        contractual_qcvars.append("CCSD+T(CCSD) TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsd+t(ccsd)":
        # contractual_qcvars.append("CCSD+T(CCSD) TOTAL GRADIENT")
        contractual_qcvars.append("CCSD+T(CCSD) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_ccsd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSD(T) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "(T) CORRECTION ENERGY",
        "CCSD(T) CORRELATION ENERGY",
        "CCSD(T) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsd(t)":
        contractual_qcvars.append("CCSD(T) TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsd(t)":
        # contractual_qcvars.append("CCSD(T) TOTAL GRADIENT")
        contractual_qcvars.append("CCSD(T) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        # print("WW", qc_module, driver, reference, method, corl_type, fcae, pv)
        expected = True
        if (
            (qc_module == "psi4-ccenergy" and method == "bccd(t)") or (qc_module == "cfour-vcc" and method == "bccd(t)")
        ) and pv in [
            "(T) CORRECTION ENERGY",
            "CCSD(T) CORRELATION ENERGY",
            "CCSD(T) TOTAL ENERGY",
        ]:
            expected = False

        yield (pv, pv, expected)


def contractual_accsd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal A-CCSD(T) (aka Lambda-CCSD(T), aka CCSD(aT) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "A-(T) CORRECTION ENERGY",
        "A-CCSD(T) CORRELATION ENERGY",
        "A-CCSD(T) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "a-ccsd(t)":
        contractual_qcvars.append("A-CCSD(T) TOTAL GRADIENT")
    elif driver == "hessian" and method == "a-ccsd(t)":
        # contractual_qcvars.append("A-CCSD(T) TOTAL GRADIENT")
        contractual_qcvars.append("A-CCSD(T) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_bccd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal BCCD(T) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "B(T) CORRECTION ENERGY",
        "BCCD(T) CORRELATION ENERGY",
        "BCCD(T) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "bccd(t)":
        contractual_qcvars.append("BCCD(T) TOTAL GRADIENT")
    elif driver == "hessian" and method == "bccd(t)":
        # contractual_qcvars.append("BCCD(T) TOTAL GRADIENT")
        contractual_qcvars.append("BCCD(T) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_cc3(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CC3 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CC3 CORRELATION ENERGY",
        "CC3 TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "cc3":
        contractual_qcvars.append("CC3 TOTAL GRADIENT")
    elif driver == "hessian" and method == "cc3":
        # contractual_qcvars.append("CC3 TOTAL GRADIENT")
        contractual_qcvars.append("CC3 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT CORRELATION ENERGY",
        "CCSDT TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt":
        contractual_qcvars.append("CCSDT TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt":
        # contractual_qcvars.append("CCSDT TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt1a(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT-1A should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT-1A CORRELATION ENERGY",
        "CCSDT-1A TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt-1a":
        contractual_qcvars.append("CCSDT-1A TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt-1a":
        # contractual_qcvars.append("CCSDT-1A TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT-1A TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt1b(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT-1B should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT-1B CORRELATION ENERGY",
        "CCSDT-1B TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt-1b":
        contractual_qcvars.append("CCSDT-1B TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt-1b":
        # contractual_qcvars.append("CCSDT-1B TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT-1B TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT-2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT-2 CORRELATION ENERGY",
        "CCSDT-2 TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt-2":
        contractual_qcvars.append("CCSDT-2 TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt-2":
        # contractual_qcvars.append("CCSDT-2 TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT-2 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt3(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    f"""Of the list of QCVariables an ideal CCSDT-3 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDT-3 CORRELATION ENERGY",
        "CCSDT-3 TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt-3":
        contractual_qcvars.append("CCSDT-3 TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt-3":
        # contractual_qcvars.append("CCSDT-3 TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT-3 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdt_prq_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSDT(Q) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "(Q) CORRECTION ENERGY",
        "CCSDT(Q) CORRELATION ENERGY",
        "CCSDT(Q) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdt(q)":
        contractual_qcvars.append("CCSDT(Q) TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdt(q)":
        # contractual_qcvars.append("CCSDT(Q) TOTAL GRADIENT")
        contractual_qcvars.append("CCSDT(Q) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_ccsdtq(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSDTQ should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSDTQ CORRELATION ENERGY",
        "CCSDTQ TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "ccsdtq":
        contractual_qcvars.append("CCSDTQ TOTAL GRADIENT")
    elif driver == "hessian" and method == "ccsdtq":
        # contractual_qcvars.append("CCSDTQ TOTAL GRADIENT")
        contractual_qcvars.append("CCSDTQ TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)


def contractual_omp2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal OMP2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "OMP2 CORRELATION ENERGY",
        "OMP2 TOTAL ENERGY",
        "OMP2 REFERENCE CORRECTION ENERGY",
        "OMP2 SAME-SPIN CORRELATION ENERGY",
        "OMP2 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "omp2":
        contractual_qcvars.append("OMP2 TOTAL GRADIENT")
    elif driver == "hessian" and method == "omp2":
        # contractual_qcvars.append("OMP2 TOTAL GRADIENT")
        contractual_qcvars.append("OMP2 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"] and method == "omp2")
        ) and pv in ["OMP2 SAME-SPIN CORRELATION ENERGY", "OMP2 OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_omp2p5(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal OMP2.5 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "OMP2.5 CORRELATION ENERGY",
        "OMP2.5 TOTAL ENERGY",
        "OMP2.5 REFERENCE CORRECTION ENERGY",
        "OMP2.5 SAME-SPIN CORRELATION ENERGY",
        "OMP2.5 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "omp2.5":
        contractual_qcvars.append("OMP2.5 TOTAL GRADIENT")
    elif driver == "hessian" and method == "omp2.5":
        # contractual_qcvars.append("OMP2.5 TOTAL GRADIENT")
        contractual_qcvars.append("OMP2.5 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"] and method == "omp2.5")
        ) and pv in ["OMP2.5 SAME-SPIN CORRELATION ENERGY", "OMP2.5 OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_omp3(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal OMP3 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "OMP3 CORRELATION ENERGY",
        "OMP3 TOTAL ENERGY",
        "OMP3 REFERENCE CORRECTION ENERGY",
        "OMP3 SAME-SPIN CORRELATION ENERGY",
        "OMP3 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "omp3":
        contractual_qcvars.append("OMP3 TOTAL GRADIENT")
    elif driver == "hessian" and method == "omp3":
        # contractual_qcvars.append("OMP3 TOTAL GRADIENT")
        contractual_qcvars.append("OMP3 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"] and method == "omp3")
        ) and pv in ["OMP3 SAME-SPIN CORRELATION ENERGY", "OMP3 OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_oremp2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal OREMP2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "OREMP2 CORRELATION ENERGY",
        "OREMP2 TOTAL ENERGY",
        "OREMP2 REFERENCE CORRECTION ENERGY",
        "OREMP2 SAME-SPIN CORRELATION ENERGY",
        "OREMP2 OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "oremp2":
        contractual_qcvars.append("OREMP2 TOTAL GRADIENT")
    elif driver == "hessian" and method == "oremp2":
        # contractual_qcvars.append("OREMP2 TOTAL GRADIENT")
        contractual_qcvars.append("OREMP2 TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            # usual for oo (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"] and method == "oremp2")
            (
                qc_module == "psi4-occ"
                and reference in ["rhf", "uhf", "rohf"]
                and corl_type in ["df", "cd"]
                and method == "oremp2"
            )
        ) and pv in ["OREMP2 SAME-SPIN CORRELATION ENERGY", "OREMP2 OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_olccd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal OLCCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "OLCCD CORRELATION ENERGY",
        "OLCCD TOTAL ENERGY",
        "OLCCD REFERENCE CORRECTION ENERGY",
        "OLCCD SAME-SPIN CORRELATION ENERGY",
        "OLCCD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "olccd":
        contractual_qcvars.append("OLCCD TOTAL GRADIENT")
    elif driver == "hessian" and method == "olccd":
        # contractual_qcvars.append("OLCCD TOTAL GRADIENT")
        contractual_qcvars.append("OLCCD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"] and method == "olccd")
        ) and pv in ["OLCCD SAME-SPIN CORRELATION ENERGY", "OLCCD OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_occd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal OCCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "OCCD CORRELATION ENERGY",
        "OCCD TOTAL ENERGY",
        "OCCD REFERENCE CORRECTION ENERGY",
        "OCCD SAME-SPIN CORRELATION ENERGY",
        "OCCD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient" and method == "occd":
        contractual_qcvars.append("OCCD TOTAL GRADIENT")
    elif driver == "hessian" and method == "occd":
        # contractual_qcvars.append("OCCD TOTAL GRADIENT")
        contractual_qcvars.append("OCCD TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                qc_module == "psi4-occ"
                and reference in ["rhf", "uhf", "rohf"]
                and corl_type in ["df", "cd"]
                and method in ["occd", "occd(t)", "a-occd(t)"]
            )
        ) and pv in ["OCCD SAME-SPIN CORRELATION ENERGY", "OCCD OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_occd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal OCCD(T) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "O(T) CORRECTION ENERGY",
        "OCCD(T) CORRELATION ENERGY",
        "OCCD(T) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "occd(t)":
        contractual_qcvars.append("OCCD(T) TOTAL GRADIENT")
    elif driver == "hessian" and method == "occd(t)":
        # contractual_qcvars.append("OCCD(T) TOTAL GRADIENT")
        contractual_qcvars.append("OCCD(T) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_aoccd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal A-OCCD(T) (aka Lambda-OCCD(T), aka OCCD(aT) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    {_contractual_docstring}
    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "A-O(T) CORRECTION ENERGY",
        "A-OCCD(T) CORRELATION ENERGY",
        "A-OCCD(T) TOTAL ENERGY",
    ]
    if driver == "gradient" and method == "a-occd(t)":
        contractual_qcvars.append("A-OCCD(T) TOTAL GRADIENT")
    elif driver == "hessian" and method == "a-occd(t)":
        # contractual_qcvars.append("A-OCCD(T) TOTAL GRADIENT")
        contractual_qcvars.append("A-OCCD(T) TOTAL HESSIAN")

    for pv in contractual_qcvars:
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_dft_current(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Given the target DFT method, returns the CURRENT QCVariables that should be produced.

    {_contractual_docstring}
    """
    contractual_qcvars = [
        (f"{method.upper()} TOTAL ENERGY", "DFT TOTAL ENERGY"),
        (f"{method.upper()} TOTAL ENERGY", "CURRENT REFERENCE ENERGY"),
        (f"{method.upper()} TOTAL ENERGY", "CURRENT ENERGY"),
        # assert missing
        ("HF TOTAL ENERGY", "HF TOTAL ENERGY"),
    ]
    if driver == "gradient":
        contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
    elif driver == "hessian":
        # contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
        contractual_qcvars.append((f"{method.upper()} TOTAL HESSIAN", "CURRENT HESSIAN"))

    for rpv, pv in contractual_qcvars:
        expected = True
        if True and rpv in [  # assert canonical wfn quantities missing in DFT
            "HF TOTAL ENERGY",
        ]:
            expected = False

        yield (rpv, pv, expected)


def contractual_dhdft_current(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str, sdsc: str
) -> Tuple[str, str, bool]:
    """Given the target DFT method, returns the CURRENT QCVariables that should be produced.

    {_contractual_docstring}
    """
    contractual_qcvars = [
        (f"{method.upper()} TOTAL ENERGY", "DFT TOTAL ENERGY"),
        (f"{method.upper()} FUNCTIONAL TOTAL ENERGY", "CURRENT REFERENCE ENERGY"),
        (f"{method.upper()} TOTAL ENERGY", "CURRENT ENERGY"),
        # assert missing
        ("HF TOTAL ENERGY", "HF TOTAL ENERGY"),
        ("MP2 CORRELATION ENERGY", "MP2 CORRELATION ENERGY"),
        ("MP2 TOTAL ENERGY", "MP2 TOTAL ENERGY"),
        ("MP2 SAME-SPIN CORRELATION ENERGY", "MP2 SAME-SPIN CORRELATION ENERGY"),
        ("MP2 SINGLES ENERGY", "MP2 SINGLES ENERGY"),
        ("MP2 DOUBLES ENERGY", "MP2 DOUBLES ENERGY"),
        ("MP2 OPPOSITE-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY"),
    ]
    if driver == "gradient":
        contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
    elif driver == "hessian":
        # contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))
        contractual_qcvars.append((f"{method.upper()} TOTAL HESSIAN", "CURRENT HESSIAN"))

    for rpv, pv in contractual_qcvars:
        expected = True
        if True and rpv in [  # assert canonical wfn quantities missing in DFT
            "HF TOTAL ENERGY",
            "MP2 CORRELATION ENERGY",
            "MP2 TOTAL ENERGY",
            "MP2 SAME-SPIN CORRELATION ENERGY",
            "MP2 SINGLES ENERGY",
            "MP2 DOUBLES ENERGY",
            "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
        ]:
            expected = False

        yield (rpv, pv, expected)
