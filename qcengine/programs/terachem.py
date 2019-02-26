"""
Calls the TeraChem executable.
"""

from qcelemental.models import ComputeError, FailedOperation, Provenance, Result

from qcengine.units import ureg

import os

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

def _output_parse(output_file):
    output_data=open(output_file).readlines()
    energy = 0
    natom = 0
    gradients = []
    spin_S2 = -1 # calculated S(S+1)
    for idx,line in enumerate(output_data):
        if "FINAL ENERGY" in line:
            energy = float(line.strip('\n').split()[2])
        elif "Total atoms" in line:
            natom = int(line.strip('\n').split()[2])
        elif "Gradient units are Hartree/Bohr" in line:
            #Gradient is stored as (dE/dx1,dE/dy1,dE/dz1,dE/dx2,dE/dy2,...)
            for i in range(idx+3,idx+3+natom):
               grad = output_data[i].strip('\n').split() 
               for x in grad:
                   gradients.append( float(x) )
        elif "SPIN S-SQUARED" in line:
            spin_S2 = float(line.strip('\n').split()[2])

    print("Energy "+ str(energy))
    print("Natom " + str(natom))
    print("Gradient " + str(gradients))
    print("S^2: " + str(spin_S2))

def _scr_parse(scrdir):
    atomic_charge = []
    atomic_charge_file = os.path.join(scrdir,"charge_mull.xls")
    if os.path.exists(atomic_charge_file):
        atomic_charge=_atomic_charge_parse(atomic_charge_file)

    print("\n")
    print(atomic_charge)

def _atomic_charge_parse(atomic_charge_file):
    atomic_charge = []
    data = open(atomic_charge_file).readlines() 
    for line in data:
        atomic_charge.append(line.strip('\n').split()[-1])
    return atomic_charge
