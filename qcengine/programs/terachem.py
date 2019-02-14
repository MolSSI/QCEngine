"""
Calls the TeraChem executable.
"""

from qcelemental.models import ComputeError, FailedOperation, Provenance, Result

from qcengine.units import ureg

def _format_input(input_data, config):
    print("")
    xyz_file = []

    s = str(len(input_data.molecule.symbols))
    xyz_file.append(s+"\n")
    for sym, geom in zip(input_data.molecule.symbols, input_data.molecule.geometry):
        s = "{:<4s} {:>{width}.{prec}f} {:>{width}.{prec}f} {:>{width}.{prec}f}".format(sym, *geom, width=14, prec=10)
        xyz_file.append(s)

    xyz_file = "\n".join(xyz_file)
    print(xyz_file)
    print("")

    input_file = []
    input_file.append("# molecule definition")
    input_file.append( "charge " + str(int(input_data.molecule.molecular_charge)))
    input_file.append( "spinmult " + str(input_data.molecule.molecular_multiplicity))
    input_file.append( "coordinates geometry.xyz")

    input_file.append("\n# model")
    input_file.append("basis " + str(input_data.model.basis))

    
    input_file.append("\n# driver")
    input_file.append("run " +input_data.driver)

    input_file.append("\n# keywords")
    for k, v in input_data.keywords.items():
        input_file.append("{} {}".format(k, v))
    print("\n".join(input_file))

def terachem(input_data, config):
    _format_input(input_data, config)

