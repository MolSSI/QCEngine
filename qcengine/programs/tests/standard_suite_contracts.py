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
            bval = qcvars_to_atomicproperties[pv] in obj

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
                    (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"])
                    or (qc_module == "gamess" and reference in ["uhf", "rohf"])
                    or (qc_module == "nwchem-tce")
                )
                and method == "mp2"
                and pv in ["MP2 SAME-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                ((qc_module == "psi4-detci"))
                and method == "mp2"
                and pv
                in [
                    "MP2 SAME-SPIN CORRELATION ENERGY",
                    "MP2 OPPOSITE-SPIN CORRELATION ENERGY",
                    "MP2 SINGLES ENERGY",
                    "MP2 DOUBLES ENERGY",
                ]
            )
            or (
                ((qc_module == "gamess" and reference in ["rhf"]) or (qc_module == "nwchem" and reference in ["rhf"]))
                and method == "ccsd"
                and pv in ["MP2 SAME-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY",]
            )
            or (
                ((qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"]))
                and method == "ccsd"
                and pv in ["MP2 SAME-SPIN CORRELATION ENERGY", "MP2 OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (
                    (qc_module == "psi4-ccenergy" and reference == "rohf")
                    or (qc_module == "nwchem-tce")
                    or (qc_module == "gamess" and reference == "rohf")
                    or (
                        qc_module.startswith("cfour") and reference == "rohf" and fcae == "fc"
                    )  # this is a cop out as c4 perfectly able to produce good rohf mp2 but not with same orbitals as ref definition on ccsd
                )
                and method == "ccsd"
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
    """For CCSD, check that the expected QCVars are present in P::e.globals and wfn and match expected ref_block."""

    contractual_qcvars = [
        "HF TOTAL ENERGY",
        "CCSD CORRELATION ENERGY",
        "CCSD TOTAL ENERGY",
        "CCSD SAME-SPIN CORRELATION ENERGY",
        "CCSD SINGLES ENERGY",
        "CCSD DOUBLES ENERGY",
        "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
    ]
    if driver == "gradient":
        contractual_qcvars.append("CCSD TOTAL GRADIENT")

    for pv in contractual_qcvars:
        expected = True
        if (
            (
                (
                    (qc_module == "gamess" and reference == "rhf")
                    or (qc_module == "nwchem-tce" and reference in ["rhf", "uhf"])
                    or (qc_module in ["cfour-ncc", "cfour-ecc"] and reference in ["rhf"])
                    or (qc_module == "psi4-occ" and reference == "rhf" and corl_type in ["df", "cd"])
                )
                and method == "ccsd"
                and pv in ["CCSD SAME-SPIN CORRELATION ENERGY", "CCSD OPPOSITE-SPIN CORRELATION ENERGY"]
            )
            or (
                (qc_module == "cfour-vcc" and reference in ["rohf"])
                and method == "ccsd"
                and pv in ["CCSD SAME-SPIN CORRELATION ENERGY", "CCSD SINGLES ENERGY", "CCSD DOUBLES ENERGY",]
            )
            or (
                (qc_module == "cfour-ecc" and reference in ["rohf"])
                and method == "ccsd"
                and pv in ["CCSD OPPOSITE-SPIN CORRELATION ENERGY", "CCSD SINGLES ENERGY", "CCSD DOUBLES ENERGY",]
            )
            or (
                (
                    (qc_module == "gamess" and reference in ["rohf"])
                    or (qc_module == "nwchem-tce" and reference in ["rohf"])
                )
                and method == "ccsd"
                and pv
                in [
                    "CCSD SAME-SPIN CORRELATION ENERGY",
                    "CCSD OPPOSITE-SPIN CORRELATION ENERGY",
                    "CCSD SINGLES ENERGY",
                    "CCSD DOUBLES ENERGY",
                ]
            )
            or (
                ()
                and method == "ccsd(t)"
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
