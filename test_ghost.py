import pprint
import re

import numpy as np
import qcelemental as qcel
import psi4
from qcelemental.testing import compare, compare_values
import qcengine as qcng





bimol_ref = {}
dmm = ["dimer", "mA", "mB", "mAgB", "gAmB"]
bimol_ref["eneyne"] = {}
bimol_ref["eneyne"]["natom"] = dict(zip(dmm, [10, 6, 4, 10, 10]))
bimol_ref["eneyne"]["nreal"] = dict(zip(dmm, [10, 6, 4, 6, 4]))
bimol_ref["eneyne"]["nre"] = dict(zip(dmm, [85.1890645313, 33.3580722134, 24.6979461998, 33.3580722134, 24.6979461998]))
bimol_ref["eneyne"]["pg"] = dict(zip(dmm, ["C2v", ["D2h", "C2v"], ["D4h", "C4v", "C2v"], "C2v", "C2v"]))
# 6-31G*
bimol_ref["eneyne"]["nbasis"] = dict(zip(dmm, [72, 38, 34, 72, 72]))
bimol_ref["eneyne"]["nmo"] = dict(zip(dmm, [72, 38, 34, 72, 72]))
bimol_ref["eneyne"]["mp2"] = dict(
    zip(dmm, [-155.3738614073, -78.2942587466, -77.0760547257, -78.2957592762, -77.0762583352])
)
# fmt: off
bimol_ref["eneyne"]["mp2_gradient"] = dict(zip(dmm, [
    np.array([
              0.000000000000,    -0.000593266109,    -0.000461056876,
             -0.000000000000,     0.000593266109,    -0.000461056876,
             -0.001061612039,     0.001517372852,     0.000162497702,
              0.001061612039,     0.001517372852,     0.000162497702,
              0.001061612039,    -0.001517372852,     0.000162497702,
             -0.001061612039,    -0.001517372852,     0.000162497702,
              0.000000000000,     0.000000000000,    -0.017632146748,
              0.000000000000,     0.000000000000,     0.018013978833,
             -0.000000000000,     0.000000000000,     0.001709075893,
              0.000000000000,     0.000000000000,    -0.001818785037]).reshape(-1, 3),
    np.array([
              0.000000000000,    -0.002441448034,     0.000180001772,
             -0.000000000000,     0.002441448034,     0.000180001772,
             -0.000901391196,     0.001441168636,    -0.000090000886,
              0.000901391196,     0.001441168636,    -0.000090000886,
              0.000901391196,    -0.001441168636,    -0.000090000886,
             -0.000901391196,    -0.001441168636,    -0.000090000886]).reshape(-1, 3),
    np.array([
              0.000000000000,     0.000000000000,    -0.016477969697,
              0.000000000000,     0.000000000000,     0.018475631244,
              0.000000000000,     0.000000000000,    -0.000114435361,
              0.000000000000,     0.000000000000,    -0.001883226187]).reshape(-1, 3),
    np.array([
              0.000000000000,    -0.001485687659,    -0.000311350386,
             -0.000000000000,     0.001485687659,    -0.000311350386,
             -0.000964480529,     0.001442449761,    -0.000001369355,
              0.000964480529,     0.001442449761,    -0.000001369355,
              0.000964480529,    -0.001442449761,    -0.000001369355,
             -0.000964480529,    -0.001442449761,    -0.000001369355,
             -0.000000000000,     0.000000000000,     0.000120621150,
             -0.000000000000,     0.000000000000,     0.000279898955,
              0.000000000000,     0.000000000000,     0.000207604252,
             -0.000000000000,     0.000000000000,     0.000020053834]).reshape(-1, 3),
    np.array([
              0.000000000000,    -0.000012504382,    -0.000103233192,
             -0.000000000000,     0.000012504382,    -0.000103233192,
             -0.000000589727,    -0.000001580375,    -0.000009047576,
              0.000000589727,    -0.000001580375,    -0.000009047576,
              0.000000589727,     0.000001580375,    -0.000009047576,
             -0.000000589727,     0.000001580375,    -0.000009047576,
              0.000000000000,     0.000000000000,    -0.016334419378,
             -0.000000000000,     0.000000000000,     0.018457229773,
             -0.000000000000,     0.000000000000,     0.000025489498,
              0.000000000000,     0.000000000000,    -0.001905643204]).reshape(-1, 3),
]))
# fmt: on


if True:
    # def test_tricky_ghost():
    driver = "energy"
    subject = "dimer"
    qcprog = "psi4"
    basis = "6-31g*"
    keywords = {"scf_type": "pk", "mp2_type": "conv"}

    # dmol = eneyne_ne_qcschemamols()["eneyne"][subject]

    dmol = {'fragments': [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9]],
            'geometry': [[ 0.        , -1.26153959, -4.01502362],
                         [ 0.        ,  1.26153959, -4.01502362],
                         [ 1.74539073, -2.32862069, -4.01790734],
                         [-1.74539073, -2.32862069, -4.01790734],
                         [-1.74539073,  2.32862069, -4.01790734],
                         [ 1.74539073,  2.32862069, -4.01790734],
                         [ 0.        ,  0.        ,  5.4811563 ],
                         [ 0.        ,  0.        ,  3.19975986],
                         [ 0.        ,  0.        ,  1.18552346],
                         [ 0.        ,  0.        ,  7.49074019]],
            'symbols': ["C", "C", "H", "H", "H", "H", "C", "C", "H", "H"],
            }

    # Freeze the input orientation so that output arrays are aligned to input
    #   and all programs match gradient.
    dmol["fix_com"] = True
    dmol["fix_orientation"] = True
    kmol = qcel.models.Molecule(**dmol)
    pmol = psi4.core.Molecule.from_schema(kmol.dict())
    ref = bimol_ref["eneyne"]
    retres = {"energy": ref["mp2"][subject], "gradient": ref["mp2_gradient"][subject]}[driver]

    assert len(kmol.symbols) == ref["natom"][subject]
    assert sum([int(at) for at in kmol.real]) == ref["nreal"][subject]

    print("CALC")
    psi4.set_options(keywords)
    ene = psi4.energy(f"mp2/{basis}", molecule=pmol)
    print("CALC FINISHED", ene)
    assert compare_values(
        retres, ene, atol=3.0e-6, label="return"
    ), f"return: {ene} != {retres}"

    print("QCENGINE CALC")
    atin = qcel.models.AtomicInput(
        **{"molecule": kmol, "model": {"method": "mp2", "basis": basis}, "driver": driver, "keywords": keywords}
    )
    atres = qcng.compute(atin, qcprog, raise_error=True)
    pprint.pprint(atres.dict(), width=200)
    print("QCENGINE CALC FINISHED")

    assert compare_values(
        ref["nre"][subject], atres.properties.nuclear_repulsion_energy, atol=1.0e-4, label="nre"
    ), f'nre: {atres.properties.nuclear_repulsion_energy} != {ref["nre"][subject]}'
    assert compare(
        ref["nbasis"][subject], atres.properties.calcinfo_nbasis, label="nbasis"
    ), f'nbasis: {atres.properties.calcinfo_nbasis} != {ref["nbasis"][subject]}'
    assert compare(
        ref["nmo"][subject], atres.properties.calcinfo_nmo, label="nmo"
    ), f'nmo: {atres.properties.calcinfo_nmo} != {ref["nmo"][subject]}'
    assert compare_values(
        ref["mp2"][subject], atres.properties.return_energy, atol=3.0e-6, label="ene"
    ), f'ene: {atres.properties.return_energy} != {ref["mp2"][subject]}'
    assert compare_values(
        retres, atres.return_result, atol=3.0e-6, label="return"
    ), f"return: {atres.return_result} != {retres}"

    pgline = {
        "cfour": r"Computational point group: (?P<pg>\w+)",
        "gamess": r"THE POINT GROUP IS (?P<pg>[\w\s,=]+)",
        "nwchem": r"Group name\s+(?P<pg>\w+)",
        "psi4": r"Running in (?P<pg>\w+) symmetry.",
    }
    mobj = re.search(pgline[qcprog], atres.stdout)
    if mobj:
        pg = mobj.group("pg").strip()
        if pg == "CNV, NAXIS= 2, ORDER= 4":
            pg = "C2v"
        elif pg == "C1 , NAXIS= 0, ORDER= 1":
            pg = "C1"
        pg = pg.capitalize()

    assert pg in ref["pg"][subject], f'pg: {pg} != {ref["pg"][subject]}'

