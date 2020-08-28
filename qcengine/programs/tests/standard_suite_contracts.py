# This file consolidates the expectations of what properties/QCVariables
#  a method should produce and what in practice a QC program does produce.

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
                # qcel.AtomicResult.properties
                vval = obj.get(qcvars_to_atomicproperties[pv])

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


def contractual_current(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Given the target method, returns the CURRENT QCVariables that should be produced.

    Parameters
    ----------
    qc_module : str
        The program or subprogram running the job (e.g., "cfour" or "cfour-ecc").
    driver : {"energy", "gradient", "hessian"}
        The derivative level that should be expected.
    reference: {"rhf", "uhf", "rohf"}
        The SCF reference since programs often output differently based on it.
    method: str
        The target AtomicInput.model.method since "free" methods may not always be
        output (e.g., MP2 available when target is MP2 but not when target is CCSD).
    corl_type: {"conv", "df", "cd"}
        The algorithm for the target method since programs often output differently
        based on it.
    fcae: {"ae", "fc"}
        The all-electron vs. frozen-orbital aspect.

    Returns
    -------
    (rpv, pv, expected)
        Of all the QCVariables `pv` that should be available, returns tuple of
        whether `expected` and what key `rpv` in the reference `pv` should match.

    """
    contractual_qcvars = [
        ("HF TOTAL ENERGY", "SCF TOTAL ENERGY"),
        ("HF TOTAL ENERGY", "CURRENT REFERENCE ENERGY"),
        (f"{method.upper()} CORRELATION ENERGY", "CURRENT CORRELATION ENERGY"),
        (f"{method.upper()} TOTAL ENERGY", "CURRENT ENERGY"),
    ]
    if driver == "gradient":
        contractual_qcvars.append((f"{method.upper()} TOTAL GRADIENT", "CURRENT GRADIENT"))

    for rpv, pv in contractual_qcvars:
        expected = True
        if method == "hf" and rpv == f"{method.upper()} CORRELATION ENERGY":
            expected = False

        yield (rpv, pv, expected)


def contractual_hf(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal HF should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    """

    contractual_qcvars = [
        ("HF TOTAL ENERGY", "HF TOTAL ENERGY"),
        ("HF TOTAL ENERGY", "SCF TOTAL ENERGY"),
    ]
    if driver == "gradient" and method == "hf":
        contractual_qcvars.append(("HF TOTAL GRADIENT", "HF TOTAL GRADIENT"))
        # contractual_qcvars.append(("HF TOTAL GRADIENT", "SCF TOTAL GRADIENT"))

    for rpv, pv in contractual_qcvars:
        yield (rpv, pv, True)


def contractual_mp2(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP2 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    Parameters
    ----------
    qc_module : str
        The program or subprogram running the job (e.g., "cfour" or "cfour-ecc").
    driver : {"energy", "gradient", "hessian"}
        The derivative level that should be expected.
    reference: {"rhf", "uhf", "rohf"}
        The SCF reference since programs often output differently based on it.
    method: str
        The target AtomicInput.model.method since "free" methods may not always be
        output (e.g., MP2 available when target is MP2 but not when target is CCSD).
    corl_type: {"conv", "df", "cd"}
        The algorithm for the target method since programs often output differently
        based on it.
    fcae: {"ae", "fc"}
        The all-electron vs. frozen-orbital aspect.

    Returns
    -------
    (rpv, pv, expected)
        Of all the QCVariables `pv` that should be available, returns tuple of
        whether `expected` and what key `rpv` in the reference `pv` should match.

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

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module == "gamess" and reference in ["uhf", "rohf"] and method == "mp2")
                    or (qc_module == "gamess" and reference in ["rhf"] and method in ["ccsd", "ccsd(t)"])
                    or (qc_module == "nwchem-tce" and method == "mp2")
                    or (qc_module == "nwchem" and reference in ["rhf"] and method in ["ccsd", "ccsd(t)"])
                    or (
                        qc_module == "psi4-occ"
                        and reference == "rhf"
                        and corl_type in ["df", "cd"]
                        and method in ["mp2", "mp2.5", "mp3", "lccd", "ccsd", "ccsd(t)"]
                    )
                )
                and pv in ["MP2 SAME-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (
                    (qc_module == "psi4-detci" and method in ["mp2", "mp3"])
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
                ((qc_module == "psi4-occ" and reference == "rohf" and method in ["olccd"]))
                and pv in ["MP2 CORRELATION ENERGY", "MP2 TOTAL ENERGY", "MP2 SINGLES ENERGY",]
            )
            or (
                (
                    (qc_module == "psi4-ccenergy" and reference == "rohf" and method == "ccsd")
                    or (qc_module == "nwchem-tce" and method in ["ccsd", "ccsd(t)"])
                    or (qc_module == "gamess" and reference == "rohf" and method == "ccsd")
                    or (
                        qc_module.startswith("cfour")
                        and reference == "rohf"
                        and fcae == "fc"
                        and method in ["ccsd", "ccsd(t)"]
                    )  # this is a cop out as c4 perfectly able to produce good rohf mp2 but not with same orbitals as ref definition on ccsd
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
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP2.5 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

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

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                qc_module == "psi4-occ"
                and reference == "rhf"
                and corl_type in ["df", "cd"]
                and method in ["mp2.5", "mp3"]
            )
            and pv in ["MP2.5 SAME-SPIN CORRELATION ENERGY", "MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"]
        ) or (
            (qc_module == "psi4-detci" and method in ["mp3"])
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
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal MP3 should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

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

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (qc_module.startswith("cfour") and method == "mp3")
                or (
                    qc_module == "psi4-occ"
                    and reference == "rhf"
                    and corl_type in ["df", "cd"]
                    and method in ["mp2.5", "mp3"]
                )
            )
            and pv in ["MP3 SAME-SPIN CORRELATION ENERGY", "MP3 OPPOSITE-SPIN CORRELATION ENERGY"]
        ) or (
            ((qc_module == "psi4-detci" and method == "mp3"))
            and pv
            in [
                "MP3 SAME-SPIN CORRELATION ENERGY",
                "MP3 OPPOSITE-SPIN CORRELATION ENERGY",
                "MP3 SINGLES ENERGY",
                "MP3 DOUBLES ENERGY",
            ]
        ):
            expected = False

        yield (pv, pv, expected)


def contractual_lccd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal LCCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    Parameters
    ----------
    qc_module : str
        The program or subprogram running the job (e.g., "cfour" or "cfour-ecc").
    driver : {"energy", "gradient", "hessian"}
        The derivative level that should be expected.
    reference: {"rhf", "uhf", "rohf"}
        The SCF reference since programs often output differently based on it.
    method: str
        The target AtomicInput.model.method since "free" methods may not always be
        output (e.g., MP2 available when target is MP2 but not when target is CCSD).
    corl_type: {"conv", "df", "cd"}
        The algorithm for the target method since programs often output differently
        based on it.
    fcae: {"ae", "fc"}
        The all-electron vs. frozen-orbital aspect.

    Returns
    -------
    (rpv, pv, expected)
        Of all the QCVariables `pv` that should be available, returns tuple of
        whether `expected` and what key `rpv` in the reference `pv` should match.

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

    for pv in contractual_qcvars:
        expected = True
        if (
            (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"] and method == "lccd")
        ) and pv in ["LCCD SAME-SPIN CORRELATION ENERGY", "LCCD OPPOSITE-SPIN CORRELATION ENERGY"]:
            expected = False

        yield (pv, pv, expected)


def contractual_lccsd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal LCCSD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    Parameters
    ----------
    qc_module : str
        The program or subprogram running the job (e.g., "cfour" or "cfour-ecc").
    driver : {"energy", "gradient", "hessian"}
        The derivative level that should be expected.
    reference: {"rhf", "uhf", "rohf"}
        The SCF reference since programs often output differently based on it.
    method: str
        The target AtomicInput.model.method since "free" methods may not always be
        output (e.g., MP2 available when target is MP2 but not when target is CCSD).
    corl_type: {"conv", "df", "cd"}
        The algorithm for the target method since programs often output differently
        based on it.
    fcae: {"ae", "fc"}
        The all-electron vs. frozen-orbital aspect.

    Returns
    -------
    (rpv, pv, expected)
        Of all the QCVariables `pv` that should be available, returns tuple of
        whether `expected` and what key `rpv` in the reference `pv` should match.

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
    if driver == "gradient" and method == "lccd":
        contractual_qcvars.append("LCCSD TOTAL GRADIENT")

    for pv in contractual_qcvars:
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_ccsd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    Parameters
    ----------
    qc_module : str
        The program or subprogram running the job (e.g., "cfour" or "cfour-ecc").
    driver : {"energy", "gradient", "hessian"}
        The derivative level that should be expected.
    reference: {"rhf", "uhf", "rohf"}
        The SCF reference since programs often output differently based on it.
    method: str
        The target AtomicInput.model.method since "free" methods may not always be
        output (e.g., MP2 available when target is MP2 but not when target is CCSD).
    corl_type: {"conv", "df", "cd"}
        The algorithm for the target method since programs often output differently
        based on it.
    fcae: {"ae", "fc"}
        The all-electron vs. frozen-orbital aspect.

    Returns
    -------
    (rpv, pv, expected)
        Of all the QCVariables `pv` that should be available, returns tuple of
        whether `expected` and what key `rpv` in the reference `pv` should match.

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

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module == "gamess" and reference == "rhf" and method in ["ccsd", "ccsd(t)"])
                    or (qc_module == "nwchem-tce" and reference in ["rhf", "uhf"] and method in ["ccsd", "ccsd(t)"])
                    or (
                        qc_module in ["cfour-ncc", "cfour-ecc"]
                        and reference in ["rhf"]
                        and method in ["ccsd", "ccsd(t)"]
                    )
                    or (
                        qc_module == "psi4-occ"
                        and reference == "rhf"
                        and corl_type in ["df", "cd"]
                        and method in ["ccsd", "ccsd(t)"]
                    )
                )
                and pv in ["CCSD SAME-SPIN CORRELATION ENERGY", "CCSD OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (qc_module == "cfour-vcc" and reference in ["rohf"] and method in ["ccsd", "ccsd(t)"])
                and pv in ["CCSD SAME-SPIN CORRELATION ENERGY", "CCSD SINGLES ENERGY", "CCSD DOUBLES ENERGY",]
            )
            or (
                (qc_module == "cfour-ecc" and reference in ["rohf"] and method in ["ccsd", "ccsd(t)"])
                and pv in ["CCSD OPPOSITE-SPIN CORRELATION ENERGY", "CCSD SINGLES ENERGY", "CCSD DOUBLES ENERGY",]
            )
            or (
                (
                    (qc_module == "gamess" and reference in ["rohf"] and method == "ccsd")
                    or (qc_module == "nwchem-tce" and reference in ["rohf"] and method in ["ccsd", "ccsd(t)"])
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
                (False)
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


def contractual_ccsd_prt_pr(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal CCSD(T) should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

    Parameters
    ----------
    qc_module : str
        The program or subprogram running the job (e.g., "cfour" or "cfour-ecc").
    driver : {"energy", "gradient", "hessian"}
        The derivative level that should be expected.
    reference: {"rhf", "uhf", "rohf"}
        The SCF reference since programs often output differently based on it.
    method: str
        The target AtomicInput.model.method since "free" methods may not always be
        output (e.g., MP2 available when target is MP2 but not when target is CCSD).
    corl_type: {"conv", "df", "cd"}
        The algorithm for the target method since programs often output differently
        based on it.
    fcae: {"ae", "fc"}
        The all-electron vs. frozen-orbital aspect.

    Returns
    -------
    (rpv, pv, expected)
        Of all the QCVariables `pv` that should be available, returns tuple of
        whether `expected` and what key `rpv` in the reference `pv` should match.

    """
    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "(T) CORRECTION ENERGY",
        "CCSD(T) CORRELATION ENERGY",
        "CCSD(T) TOTAL ENERGY",
    ]
    if driver == "gradient":
        contractual_qcvars.append("CCSD(T) TOTAL GRADIENT")

    for pv in contractual_qcvars:
        # print("WW", qc_module, driver, reference, method, corl_type, fcae, pv)
        expected = True
        if False:
            expected = False

        yield (pv, pv, expected)


def contractual_olccd(
    qc_module: str, driver: str, reference: str, method: str, corl_type: str, fcae: str
) -> Tuple[str, str, bool]:
    """Of the list of QCVariables an ideal OLCCD should produce, returns whether or
    not each is expected, given the calculation circumstances (like QC program).

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

    for pv in contractual_qcvars:
        expected = True

        yield (pv, pv, expected)
