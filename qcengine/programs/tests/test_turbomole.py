import json
import re
from subprocess import Popen, PIPE

import qcelemental
import qcengine
import pytest


def test_turbomole():
    mol = qcelemental.models.Molecule.from_data("""
            O 0.000000000000     0.000000000000    -0.068516245955
            H 0.000000000000    -0.790689888800     0.543701278274
            H 0.000000000000     0.790689888800     0.543701278274
    """)

    task = {
            "schema_name": "qcschema_input",
            "schema_version": 1,
            "molecule": mol,
            "driver": "energy",
            "model" : {"method": "SCF", "basis": "sto-3g"},
            "keywords": {},
    }

    # from qcengine.programs.turbomole import TurbomoleHarness
    # tm = TurbomoleHarness()
    # v = tm.get_version()
    # print(v)
    res = qcengine.compute(task, "turbomole")
    print(res)


def test_harvest():
    from qcengine.programs.turbomole.harvester import harvest
    with open("dscf.out") as handle:
        text = handle.read()
    res = harvest(None, text)
    print(res)


def wrap_define():
    charge = 0
    mult = 1
    unrestricted = True

    def occ_num_mo_data(charge, mult, unrestricted=False):
        """Handles the 'Occupation Number & Molecular Orbital' section
        of define."""
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

    kwargs = {
        "charge": 0,
        "init_guess": occ_num_mo_data_stdin(charge, mult, unrestricted),
        "title": "QCEngine Turbomole",
        "xcfunc": "pbe0",
        "grid": "m4",
        "scf_conv": 8,
        "scf_iters": 150,
    }

    stdin = """
    {title}
    a coord
    *
    no
    b
    all def2-SVP
    *
    {init_guess}
    dft
    on
    func
    {xcfunc}
    grid
    {grid}

    scf
    conv
    {scf_conv}
    iter
    {scf_iters}

    *
    """.format(**kwargs)

    proc = Popen("define", universal_newlines=True,
                 stdin=PIPE, stdout=PIPE, stderr=PIPE
    )
    stdout, _ = proc.communicate(stdin)#, timeout=3)
    proc.terminate()

# @pytest.mark.parametrize(
    # 'program,basis,keywords',
    # [
        # pytest.param('turbomole', 'aug-pvdz', {}),
    # ])  # yapf: disable
# def test_sp_hf_rhf(program, basis, keywords):
    # mol = qcelemental.models.Molecule.from_data("""
            # O 0.000000000000     0.000000000000    -0.068516245955
            # H 0.000000000000    -0.790689888800     0.543701278274
            # H 0.000000000000     0.790689888800     0.543701278274
    # """)

    # resi = {
        # "molecule": mol,
        # "driver": "energy",
        # "model": {
            # "method": "hf",
            # "basis": basis
        # },
        # "keywords": keywords,
    # }

    # res = qcengine.compute(resi, program, raise_error=True, return_dict=True)

    # assert res["driver"] == "energy"
    # assert "provenance" in res
    # assert res["success"] is True

    # # aug-cc-pvdz
    # scf_tot = -76.0413815332

    # atol = 1.e-6
    # assert compare_values(scf_tot, res["return_result"], atol=atol)


if __name__ == "__main__":
    test_turbomole()
    # test_harvest()
    # wrap_define()
