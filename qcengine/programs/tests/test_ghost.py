import pprint
import re

import pytest
import numpy as np
import qcelemental as qcel
from qcelemental.testing import compare, compare_values

import qcengine as qcng
from qcengine.programs.tests.test_dftd3_mp2d import eneyne_ne_qcschemamols
from qcengine.testing import using


@pytest.fixture
def hene():
    smol = """
 0 1
 He 0 0 0
 @Ne 2.5 0 0
 nocom
 noreorient
"""
    return qcel.models.Molecule.from_data(smol)


@pytest.mark.parametrize("driver", ["energy", "gradient"])
@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param("cfour", "aug-pvdz", {}, marks=using("cfour")),
        pytest.param("gamess", "ccd", {"contrl__ispher": 1}, marks=using("gamess")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "scf__thresh": 1.0e-6}, marks=using("nwchem")),
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"basis__spherical": True, "scf__thresh": 1.0e-6, "qc_module": "tce"},
            marks=using("nwchem"),
        ),
        pytest.param("psi4", "aug-cc-pvdz", {"scf_type": "direct"}, marks=using("psi4")),
    ],
)
def test_simple_ghost(driver, program, basis, keywords, hene):
    resi = {"molecule": hene, "driver": driver, "model": {"method": "hf", "basis": basis}, "keywords": keywords}

    if program == "gamess":
        with pytest.raises(qcng.exceptions.InputError) as e:
            qcng.compute(resi, program, raise_error=True, return_dict=True)
        pytest.xfail("no ghosts with gamess")

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == driver
    assert "provenance" in res
    assert res["success"] is True

    pprint.pprint(res, width=200)

    atol = 1.0e-6
    assert compare_values(0.0, res["properties"]["nuclear_repulsion_energy"], atol=atol, label="nre")
    assert compare(32, res["properties"]["calcinfo_nbasis"], label="nbas")
    assert compare(32, res["properties"]["calcinfo_nmo"], label="nmo")

    ene = -2.8557143339397539
    grad = np.array(
        [
            0.000004317771,
            0.000000000000,
            0.000000000000,
            -0.000004317771,
            0.000000000000,
            0.000000000000,
        ]
    )
    retres = {"energy": ene, "gradient": grad}[driver]
    assert compare_values(ene, res["properties"]["return_energy"], atol=atol, label="ene")
    assert compare_values(retres, res["return_result"], atol=atol, label="return")


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


@pytest.mark.parametrize("driver", ["energy", "gradient"])
@pytest.mark.parametrize("subject", dmm)
@pytest.mark.parametrize(
    "qcprog, basis, keywords",
    [
        pytest.param("cfour", "6-31g*", {"spherical": 0, "scf_conv": 12}, id="cfour", marks=using("cfour")),
        pytest.param(
            "gamess",
            "n31",
            {"basis__ngauss": 6, "basis__ndfunc": 1, "mp2__nacore": 0},
            id="gamess",
            marks=using("gamess"),
        ),
        pytest.param("nwchem", "6-31g*", {"scf__thresh": 1.0e-8}, id="nwchem", marks=using("nwchem")),
        pytest.param(
            "psi4",
            "6-31g*",
            {
                "scf_type": "pk",
                "mp2_type": "conv",
            },
            id="psi4",
            marks=using("psi4"),
        ),
    ],
)
def test_tricky_ghost(driver, qcprog, subject, basis, keywords):
    dmol = eneyne_ne_qcschemamols()["eneyne"][subject]
    # Freeze the input orientation so that output arrays are aligned to input
    #   and all programs match gradient.
    dmol["fix_com"] = True
    dmol["fix_orientation"] = True
    if qcprog == "nwchem" and subject in ["mA", "mB"]:
        # NWChem symmetry detection uses a loose tolerance (0.01) that catches
        #   the monomers and symmetrizes them, belying the atoms_map=True arg
        #   to the in_mol vs. calc_mol aligner. Gradients don't change much
        #   w/symm geometry but enough that the symmetrizer must be defeated.
        keywords["geometry__autosym"] = "1d-4"
    kmol = qcel.models.Molecule(**dmol)
    ref = bimol_ref["eneyne"]

    assert len(kmol.symbols) == ref["natom"][subject]
    assert sum([int(at) for at in kmol.real]) == ref["nreal"][subject]

    atin = qcel.models.AtomicInput(
        **{"molecule": kmol, "model": {"method": "mp2", "basis": basis}, "driver": driver, "keywords": keywords}
    )

    if qcprog == "gamess" and subject in ["mAgB", "gAmB"]:
        with pytest.raises(qcng.exceptions.InputError) as e:
            res = qcng.compute(atin, qcprog, raise_error=True)
        pytest.xfail("no ghosts with gamess")

    atres = qcng.compute(atin, qcprog)
    pprint.pprint(atres.dict(), width=200)

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
    retres = {"energy": ref["mp2"][subject], "gradient": ref["mp2_gradient"][subject]}[driver]
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

    if qcprog == "gamess":  # and subject in ["mAgB", "gAmB"]:
        # don't know how to get master frame w/ghosts in gamess, so C1 forced
        assert pg == "C1", f"pg: {pg} != C1"
    else:
        assert pg in ref["pg"][subject], f'pg: {pg} != {ref["pg"][subject]}'


@pytest.mark.parametrize(
    "qcprog, basis, keywords",
    [
        pytest.param("cfour", "aug-pvdz", {"scf_conv": 12}, id="cfour", marks=using("cfour")),
        pytest.param(
            "gamess",
            "accd",
            {"mp2__nacore": 0, "contrl__ispher": 1},
            id="gamess",
            marks=using("gamess"),
        ),
        pytest.param(
            "nwchem", "aug-cc-pvdz", {"scf__nr": 1.0, "scf__thresh": 1.0e-8}, id="nwchem", marks=using("nwchem")
        ),
        pytest.param(
            "psi4",
            "aug-cc-pvdz",
            {
                "scf_type": "pk",
                "mp2_type": "conv",
            },
            id="psi4",
            marks=using("psi4"),
        ),
    ],
)
def test_atom_labels(qcprog, basis, keywords):
    kmol = qcel.models.Molecule.from_data(
        """
      H       0 0 0
      H5      5 0 0
      H_other 0 5 0
      H_4sq   5 5 0
      units au
    """
    )

    assert compare(["H", "H", "H", "H"], kmol.symbols, "elem")
    assert compare(["", "5", "_other", "_4sq"], kmol.atom_labels, "elbl")

    atin = qcel.models.AtomicInput(
        **{"molecule": kmol, "model": {"method": "mp2", "basis": basis}, "driver": "energy", "keywords": keywords}
    )

    atres = qcng.compute(atin, qcprog)
    pprint.pprint(atres.dict(), width=200)

    nre = 1.0828427
    assert compare_values(
        nre, atres.properties.nuclear_repulsion_energy, atol=1.0e-4, label="nre"
    ), f"nre: {atres.properties.nuclear_repulsion_energy} != {nre}"

    nmo = 36
    assert compare(nmo, atres.properties.calcinfo_nmo, label="nmo"), f"nmo: {atres.properties.calcinfo_nmo} != {nmo}"

    scf = -1.656138508
    scf_alt = -1.705577613  # some versions of NWChem land on this value
    mp2 = -1.7926264513
    mp2_alt = -1.870251459939
    try:
        assert compare_values(
            scf, atres.properties.scf_total_energy, atol=3.0e-6, label="scf ene"
        ), f"scf ene: {atres.properties.scf_total_energy} != {scf}"
    except AssertionError as exc:
        if qcprog == "nwchem":
            assert compare_values(
                scf_alt, atres.properties.scf_total_energy, atol=3.0e-6, label="scf ene"
            ), f"scf ene: {atres.properties.scf_total_energy} != {scf_alt}"
            assert compare_values(
                mp2_alt, atres.return_result, atol=3.0e-6, label="ene"
            ), f"ene: {atres.return_result} != {mp2_alt}"
        else:
            raise AssertionError from exc
    else:
        assert compare_values(
            mp2, atres.return_result, atol=3.0e-6, label="ene"
        ), f"ene: {atres.return_result} != {mp2}"
