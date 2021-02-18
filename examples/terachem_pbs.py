import os
import qcengine as qcng
import qcelemental as qcel

os.environ["TERACHEM_PBS_HOST"] = "127.0.0.1"
os.environ["TERACHEM_PBS_PORT"] = "11111"

prog = qcng.get_program("terachem_pbs")


mol = qcel.models.Molecule.from_data(
    """
    O  0.0  0.000  -0.129
    H  0.0 -1.494  1.027
    H  0.0  1.494  1.027
"""
)

inp = qcel.models.AtomicInput(
    molecule=mol,
    driver="energy",
    model={"method": "pbe0", "basis": "6-31g"},
)
ret = prog.compute(inp)
print(ret)