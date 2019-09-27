from subprocess import Popen, PIPE, TimeoutExpired
from typing import Any, Dict, Optional

from .methods import ricc2_methods


def execute_define(stdin: str, cwd: Optional["Path"] = None) -> str:
    # TODO: replace this with a call to the default execute provided by QCEngine
    # if possible. May be difficult though, as we have to pipe in stdin and
    # be careful with the encoding.

    # We cant use univeral_newlines=True or text=True in Popen as some of the
    # data that define returns isn't proper UTF-8, so the decoding will crash.
    # We will decode it later on manually.
    proc = Popen("define",
                 stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=cwd,
                 # stdin=PIPE, stderr=PIPE, cwd=cwd,
    )
    # TODO: add timeout? Unless the disk hangs this should never take long...
    # TODO: how to get the stdout when the process hangs? Maybe write it to a file?
    stdout, _ = proc.communicate(str.encode(stdin))
    proc.terminate()
    try:
        stdout = stdout.decode("utf-8")
    except UnicodeDecodeError:
        # Some of the basis files (cbas, I'm looking at you ...) are saved
        # in ISO-8859-15 but most of them are in UTF-8. Decoding will
        # crash in the former cases so here we try the correct decoding.
        stdout = stdout.decode("latin-1")

    # try:
        # stdout, _ = proc.communicate(str.encode(stdin), timeout=5)
        # proc.terminate()
    # except TimeoutExpired as err:
        # print(err)
        # import pdb; pdb.set_trace()

    return stdout


def prepare_stdin(method: str, basis: str, keywords: Dict[str, Any],
                  charge: int, mult: int, geoopt: Optional[str] = "")  -> str:
    """Prepares a str that can be sent to define to produce the desired
    input for Turbomole."""

    # Load data from keywords
    unrestricted = keywords.get("unrestricted", False)
    grid = keywords.get("grid", "m3")

    def occ_num_mo_data(charge: int, mult: int,
                        unrestricted: Optional[bool] = False) -> str:
        """Handles the 'Occupation Number & Molecular Orbital' section
        of define. Sets appropriate charge and multiplicity in the
        system and decided between restricted and unrestricted calculation.

        RHF and UHF are supported. ROHF could be implemented later on
        by using the 's' command to list the available MOs and then
        close the appropriate number of MOs to doubly occupied MOs
        by 'c' by comparing the number of total MOs and the desired
        multiplicity."""

        # Do unrestricted calculation if explicitly requested or mandatory
        unrestricted = unrestricted or (mult != 1)
        unpaired = mult - 1
        charge = int(charge)

        occ_num_mo_data_stdin = f"""eht
        y
        {charge}
        y
        """
        if unrestricted:
            # Somehow Turbomole/define asks us if we want to write
            # natural orbitals... we don't want to.
            occ_num_mo_data_stdin = f"""eht
            y
            {charge}
            n
            u {unpaired}
            *
            n
            """
        return occ_num_mo_data_stdin

    def set_method(method, grid):
        if method == "hf":
            method_stdin = ""
        elif method in ricc2_methods:
            # Setting geoopt in $ricc2 will make the ricc2 module to produce
            # a gradient.
            # Drop the 'ri'-prefix of the method string.
            geoopt_stdin = f"geoopt {method[2:]} ({geoopt})" if geoopt else ""
            method_stdin = f"""cc
                               freeze
                               *
                               cbas
                               *
                               ricc2
                               {method}
                               list models

                               {geoopt_stdin}
                               list geoopt

                               *
                               *
                            """
        # Fallback: assume method corresponds to a DFT functional
        #
        # TODO: handle xcfuncs that aren't defined in define, e.g.
        # new functionals introduced in 7.4 from libxc. ...
        # Maybe the best idea would be to not set the functional here
        # but just turn on DFT and add it to the control file later on.
        else:
            method_stdin = f"""dft
                               on
                               func
                               {method}
                               grid
                               {grid}

                            """
        return method_stdin

    def set_ri():
        ri_stdin = ""
        return ri_stdin

    kwargs = {
        "init_guess": occ_num_mo_data(charge, mult, unrestricted),
        "set_method": set_method(method, grid),
        "title": "QCEngine Turbomole",
        "scf_conv": 8,
        "scf_iters": 150,
        "basis": basis,
    }

    stdin = """
    {title}
    a coord
    *
    no
    b
    all {basis}
    *
    {init_guess}
    {set_method}
    scf
    conv
    {scf_conv}
    iter
    {scf_iters}

    *
    """.format(**kwargs)

    return stdin
