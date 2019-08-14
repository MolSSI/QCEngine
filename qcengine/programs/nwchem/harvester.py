import re
from decimal import Decimal

import qcelemental as qcel
from qcelemental.models import Molecule

from ..util import PreservingDict

def harvest_output(outtext):
    """Function to separate portions of a NWChem output file *outtext*,
    divided by " Line search:".

    """
    pass_psivar = []
    pass_coord = []
    pass_grad = []
    print(pass_coord)
    for outpass in re.split(r' Line search:', outtext, re.MULTILINE):
        psivar, nwcoord, nwgrad, version, error = harvest_outfile_pass(outpass)
        pass_psivar.append(psivar)
        pass_coord.append(nwcoord)
        pass_grad.append(nwgrad)


#       print ('\n\nXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n')
        #print (outpass)
        #print (psivar, nwcoord, nwgrad)
        #print (psivar, nwgrad)
#       print ('\n\nxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n')

    retindx = -1 if pass_coord[-1] else -2

    #    print ('    <<<  NW PSIVAR  >>>')
    #    for item in pass_psivar[retindx]:
    #       print('       %30s %16.8f' % (item, pass_psivar[retindx][item]))
    #    print ('    <<<  NW COORD   >>>')
    #    for item in pass_coord[retindx]:
    #    print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
    #    print ('    <<<   NW GRAD   >>>')
    #    for item in pass_grad[retindx]:
    #       print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

    return pass_psivar[retindx], pass_coord[retindx], pass_grad[retindx], version, error

def harvest_outfile_pass(outtext):
    """Function to read NWChem output file *outtext* and parse important 
    quantum chemical information from it in

    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None
    version = ''
    error = ''

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    # Process version
    mobj = re.search(
        r'^\s+' + r'Northwest Computational Chemistry Package (NWChem)' + r'\s+' + r'(?:<version>\d+.\d+)' + r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched version')
        version = mobj.group('version')

    #Process SCF
    #1)Fail to converge
    mobj = re.search(r'^\s+' + r'(?:Calculation failed to converge)' + r'\s*$', outtext, re.MULTILINE)
    if mobj:
        print('failed to converge')

    #2)Calculation converged
    else:
        mobj = re.search(r'^\s+' + r'(?:Total SCF energy)' + r'\s+=\s*' + NUMBER + r's*$', outtext, re.MULTILINE)
        if mobj:
            print('matched HF')
            psivar['HF TOTAL ENERGY'] = mobj.group(1)

    #Process Effective nuclear repulsion energy (a.u.)
        mobj = re.search(r'^\s+' + r'Effective nuclear repulsion energy \(a\.u\.\)' + r'\s+' + NUMBER + r'\s*$',
                         outtext, re.MULTILINE)
        if mobj:
            print('matched NRE')
            #print (mobj.group(1))
            psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)

        #Process DFT (RDFT, RODFT,UDFT, SODFT [SODFT for nwchem versions before nwchem 6.8])
        mobj = re.search(
                r'^\s+' + r'(?:Total DFT energy)' + r'\s+=\s*' + NUMBER + r'\s*$',
                outtext, re.MULTILINE)
        if mobj:
            print('matched DFT')
            print (mobj.group(1))
            psivar['DFT TOTAL ENERGY'] = mobj.group(1)

        #SODFT [for nwchem 6.8+]
        mobj = re.search(
            r'^\s+' + r'Total SO-DFT energy' + r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'Nuclear repulsion energy' +
            r'\s+' + NUMBER + r'\s*$', outtext, re.MULTILINE)
        if mobj:
            print('matched DFT')
            #print (mobj.group(1))
            psivar['DFT TOTAL ENERGY'] = mobj.group(1)
            psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(2)

        #MCSCF
        sym = []
        mobj = re.findall(
                r'^\s+' + r'Total SCF energy' + r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'One-electron energy' +
                r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'Two-electron energy' + r'\s+' + NUMBER + r'\s*' +
                r'^\s+' + r'Total MCSCF energy' + r'\s+' + NUMBER + r'\s*$', 
                outtext, re.MULTILINE | re.DOTALL)

        #for mobj_list in mobj:

        if mobj:# Need to change to accommodate find all instances
            print('matched mcscf')  #MCSCF energy calculation
            psivar['HF TOTAL ENERGY'] = mobj.group(1)
            psivar['ONE-ELECTRON ENERGY'] = mobj.group(2)
            psivar['TWO-ELECTRON ENERGY'] = mobj.group(3)
            psivar['MCSCF TOTAL ENERGY'] = mobj.group(4)
	#for mobj_list in mobj:
	    #for i in mobj_list:
		#count += 0
		#print('matched mcscf iteration %i', count)
            	#psivar['HF TOTAL ENERGY'] = mobj.group(1)
            	#psivar['ONE-ELECTRON ENERGY'] = mobj.group(2)
            	#psivar['TWO-ELECTRON ENERGY'] = mobj.group(3)
            	#psivar['MCSCF TOTAL ENERGY'] = mobj.group(4)

        #Process MP2 (Restricted, Unrestricted(RO n/a))
        #1)SCF-MP2
        mobj = re.search(
            r'^\s+' + r'SCF energy' + r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'Correlation energy' + r'\s+' + NUMBER +
            r'\s*' + r'^\s+' + r'Singlet pairs' + r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'Triplet pairs' + r'\s+' +
            NUMBER + r'\s*' + r'^\s+' + r'Total MP2 energy' + r'\s+' + NUMBER + r'\s*$', outtext, re.MULTILINE)  #MP2
        if mobj:
            print('matched scf-mp2')
            psivar['HF TOTAL ENERGY'] = mobj.group(1)
            psivar['MP2 CORRELATION ENERGY'] = mobj.group(2)
            psivar['MP2 TOTAL ENERGY'] = mobj.group(5)
        #SCS-MP2
        mobj = re.search(
            r'^\s+' + r'Same spin pairs' + r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'Same spin scaling factor' + r'\s+' +
            NUMBER + r'\s*' + r'^\s+' + r'Opposite spin pairs' + r'\s+' + NUMBER + r'\s*' + r'^\s+' +
            r'Opposite spin scaling fact.' + r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'SCS-MP2 correlation energy' +
            r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'Total SCS-MP2 energy' + r'\s+' + NUMBER + r'\s*$', outtext,
            re.MULTILINE)
        if mobj:
            print('matched scs-mp2')
            psivar['MP2 SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) * Decimal(mobj.group(2))
            psivar['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(3)) * Decimal(mobj.group(4))
            
            print(mobj.group(1))  #ess
            print(mobj.group(2))  #fss
            print(mobj.group(3))  #eos
            print(mobj.group(4))  #fos
            print(mobj.group(5))  #scs corl
            print(mobj.group(6))  #scs-mp2

        #2) DFT-MP2
        mobj = re.search(
            r'^\s+' + r'DFT energy' + r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'Unscaled MP2 energy' + r'\s+' + NUMBER +
            r'\s*' + r'^\s+' + r'Total DFT+MP2 energy' + r'\s+' + NUMBER + r'\s*$', outtext, re.MULTILINE)
        if mobj:
            print('matched dft-mp2')
            psivar['DFT TOTAL ENERGY'] = mobj.group(1)
            psivar['MP2 CORRELATION ENERGY'] = mobj.group(2)
            psivar['MP2 TOTAL ENERGY'] = mobj.group(3)

        #3) MP2 with CCSD or CCSD(T) calculation (through CCSD(T) directive)
        mobj = re.search(
            r'^\s+' + r'MP2 Energy \(coupled cluster initial guess\)' + r'\s*' + r'^\s+' +
            r'------------------------------------------' + r'\s*' + r'^\s+' + r'Reference energy:' + r'\s+' + NUMBER +
            r'\s*' + r'^\s+' + r'MP2 Corr\. energy:' + r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'Total MP2 energy:' +
            r'\s+' + NUMBER + r'\s*$', outtext, re.MULTILINE)

        if mobj:
            print('matched coupled cluster-mp2')
            psivar['MP2 CORRELATION ENERGY'] = mobj.group(2)
            psivar['MP2 TOTAL ENERGY'] = mobj.group(3)

        #4) Direct MP2

        #5) RI-MP2

        #Process calculation through tce [dertype] command
        cc_name = ''
        for cc_name in [r'MBPT\(2\)', r'MBPT\(3\)', r'MBPT\(4\)']:
            mobj = re.search(
                r'^\s+' + cc_name + r'\s+' + r'correlation energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*' +
                r'^\s+' + cc_name + r'\s+' + r'total energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*$', outtext, re.MULTILINE)
        
            if mobj:
                mbpt_plain = cc_name.replace('\\', '').replace('MBPT', 'MP').replace('(', '').replace(')', '')
                print(f'matched tce mbpt {mbpt_plain}', mobj.groups())

                if mbpt_plain == 'MP2':
                    psivar[f'{mbpt_plain} CORRELATION ENERGY'] = mobj.group(1)
                else:
                    psivar[f'{mbpt_plain} CORRECTION ENERGY'] = mobj.group(1)
                psivar[f'{mbpt_plain} TOTAL ENERGY'] = mobj.group(2)

        # Process CC '()' correction part through tce [dertype] command
        for cc_name in [r'CCSD\(T\)', r'CCSD\[T\]']:
            mobj = re.search(
                r'^\s+' + cc_name + r'\s+' + r'correction energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*' +
                r'^\s+' + cc_name + r'\s+' + r'correlation energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*' +
                r'^\s+' + cc_name + r'\s+' + r'total energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*$', outtext, re.MULTILINE)

            if mobj:
                cc_plain = cc_name.replace('\\', '')
                cc_corr = cc_plain.replace('CCSD', '')
                print(f'matched tce cc {cc_plain}')

                psivar[f'{cc_corr} CORRECTION ENERGY'] = mobj.group(1)
                psivar[f'{cc_plain} CORRELATION ENERGY'] = mobj.group(2)
                psivar[f'{cc_plain} TOTAL ENERGY'] = mobj.group(3)
        #Process other TCE cases
        for cc_name in [r'CISD', r'CISDT', r'CISDTQ', r'CCD', r'CCSD', r'CCSDT', r'CCSDTQ', r'LCCSD', r'LCCD']:
            mobj = re.search(
                r'^\s+' + r'Iterations converged' + r'\s*' +
                r'^\s+' + cc_name + r'\s+' + r'correlation energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*'+
                r'^\s+' + cc_name + r'\s+' + r'total energy / hartree' + r'\s+=\s*' + NUMBER + r'\s*$', 
                outtext, re.MULTILINE)

            if mobj:
                print(f'matched {cc_name}')
                print(mobj)
                psivar[f'{cc_name} CORRELATION ENERGY'] = mobj.group(1)
                psivar[f'{cc_name} TOTAL ENERGY'] = mobj.group(2)

        #Process CCSD/CCSD(T) using nwchem CCSD/CCSD(T) [dertype] command
        mobj = re.search(
            r'^\s+' + r'-----------' + r'\s*' + r'^\s+' + r'CCSD Energy' + r'\s*' + r'^\s+' + r'-----------' + r'\s*' +
            r'^\s+' + r'Reference energy:' + r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'CCSD corr\. energy:' + r'\s+' +
            NUMBER + r'\s*' + r'^\s+' + r'Total CCSD energy:' + r'\s+' + NUMBER + r'\s*$', outtext,
            re.MULTILINE | re.DOTALL)

        if mobj:
            print('matched ccsd')
            psivar['CCSD CORRELATION ENERGY'] = mobj.group(2)
            psivar['CCSD TOTAL ENERGY'] = mobj.group(3)

        mobj = re.search(
            r'^\s+' + r'--------------' + r'\s*' + r'^\s+' + r'CCSD\(T\) Energy' + r'\s*' + r'^\s+' + r'--------------'
            + r'\s*' + r'(?:.*?)' + r'^\s+' + r'\(T\) corr\. energy:' + r'\s+' + NUMBER + r'\s*' + r'^\s+' +
            r'Total CCSD\(T\) energy:' + r'\s+' + NUMBER + r'\s*$', outtext, re.MULTILINE | re.DOTALL)

        if mobj:
            print('matched ccsd(t)')
            psivar['(T) CORRECTION ENERGY'] = mobj.group(1)
            psivar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(2)) - psivar['HF TOTAL ENERGY']
            psivar['CCSD(T) TOTAL ENERGY'] = mobj.group(2)

        mobj = re.search(
            r'^\s+' + r'Spin Component Scaled (SCS) CCSD' + r'\s*' + r'(?:.*?)' + r'^\s+' + r'Same spin contribution' +
            r'\s+' + NUMBER + r'\s*' + r'^\s+' + r'Same spin scaling factor' + r'\s+' + NUMBER + r'\s*'
            r'^\s+' + r'Opposite spin contribution' + r'\s+' + NUMBER + r'\s*' + r'^\s+' +
            r'Opposite spin scaling factor' + r'\s+' + NUMBER + r'\s*'
            r'\^s+' + r'SCS-CCSD correlation energy' + r'\s+' + NUMBER + r'\s*' + r'\^s+' + r'Total SCS-CCSD energy' +
            r'\s+' + NUMBER + r'\s*$', outtext, re.MULTILINE | re.DOTALL)
        #SCS-CCSD included
        if mobj:
            print('matched scs-ccsd')
            psivar['CCSD SAME-SPIN CORRELATION ENERGY'] = psivar['SCS-CCSD SAME-SPIN CORRELATION ENERGY'] = (
                Decimal(mobj.group(1)) * Decimal(mobj.group(2)))
            psivar['CCSD OPPOSITE-SPIN CORRELATION ENERGY'] = psivar['SCS-CCSD OPPOSITE-SPIN CORRELATION ENERGY'] = (
                Decimal(mobj.group(4)) * Decimal(mobj.group(3)))
            psivar['SCS-CCSD CORRELATION ENERGY'] = mobj.group(5)
            psivar['SCS-CCSD TOTAL ENERGY'] = mobj.group(6)
            psivar['CUSTOM SCS-CCSD CORRELATION ENERGY'] = 0.5 * (float(
                psivar['CCSD SAME-SPIN CORRELATION ENERGY']) + float(psivar['CCSD OPPOSITE-SPIN CORRELATION ENERGY']))
            #psivar['CUSTOM SCS-CCSD TOTAL ENERGY'] = float(mobj.group(6)) + float(
            #    psivar['CUSTOM SCS-CCSD CORRERLATION ENERGY'])

        #Process EOM-[cc_name] #nwchem_tce_dipole = false
        # Parsed information: each symmetry, root excitation energy in eV and total energy in hartree
        # psivar name might need to be fixed
        # each root excitation energy is extracted from the last iteration of right hand side
        mobj = re.findall(
            #r'^\s+' + r'Excited-state calculation' + r'\s+' + r'\(\s+' + r'(.*?)' + r'\s+' + r'symmetry\)' + r'\s*' +
            #r'^\s+' + r'=========================================' + r'\s*' + r'(?:.*?)'
            #+  #Skip to line before right-hand side iterations
            #r'^\s+' + r'EOM-' + cc_name + r' right-hand side iterations' + r'\s*' + r'(?:.*?)' + r'^\s+' + r'Iteration'
            #+ r'\s+\d+\s+' + r'using' + r'\s+\d+\s+' + r'trial vectors' + r'\s*' +
            #r'((?:\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+(?:(\s+\d+\.\d+\s+\d+\.\d+)?)\s*\n)+)' + r'^\s+' +
            #r'--------------------------------------------------------------' + r'\s*' + r'^\s+' +
            #r'Iterations converged' + r'\s*$',
            #outtext,
            #re.MULTILINE | re.DOTALL)
            r'^\s+' + r'Excited-state calculation' + r'\s+' + r'\((.*?)\)' + r'\s*' +
            r'^\s+' + r'EOM-' + cc_name + r'\s+right-hand side iterations' + r'\s*'+  #capture cc_name
            r'^\s+' + r'Iterations converged' + r'\s*' + #skip to excitation
            r'^\s+' + r'Excited state root\s+\d{1,2}\s*'+
            r'^\s+' + r'Excitation energy / hartree' + r'\s+=\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'/ eV\s+' + r'\s+=\s+' + NUMBER + r'\s*$',
            outtext, re.MULTILINE | re.DOTALL)
        
        if mobj:
            ext_energy = {}  #dic

            ext_energy_list = []
            for mobj_list in mobj:
                print('matched EOM-%s - %s symmetry' % (cc_name, mobj_list[0]))  #cc_name, symmetry
                print(mobj_list)
                count = 0
                for line in mobj_list[1].splitlines():
                    lline = line.split()
                    print(lline[1])  #in hartree
                    print(lline[2])  #in eV
                    count += 1

                    print('matched excitation energy #%d - %s symmetry' % (count, mobj_list[0]))

                    ext_energy_list.append(lline[1])  #Collect all excitation energies

                    sym = str(mobj_list[0])
                    ext_energy.setdefault(sym, [])
                    ext_energy[sym].append(lline[1])  #Dictionary: symmetries(key), energies(value)

            ext_energy_list.sort(key=float)

            for nroot in range(len(ext_energy_list)):
                for k, e_val in ext_energy.items():
                    if ext_energy_list[nroot] in e_val:
                        symm = k
                        psivar ['EOM-%s ROOT 0 -> ROOT %d EXCITATION ENERGY - %s SYMMETRY' %(cc_name, nroot+1, symm)] = \
                            ext_energy_list[nroot] #in hartree
                        psivar ['EOM-%s ROOT 0 -> ROOT %d TOTAL ENERGY - %s SYMMETRY' %(cc_name, nroot+1, symm)] = \
                            psivar['%s TOTAL ENERGY' %(cc_name)] + Decimal(ext_energy_list[nroot]) #in hartree
        gssym = ''
        gs = re.search(
            r'^\s+' + r'Ground-state symmetry is' + gssym + r'\s*$', outtext, re.MULTILINE)
        
        if gs:
            print('matched ground-state symmetry')
            psivar['GROUND-STATE SYMMETRY'] = gssym.group(1)

#Process TDDFT
#       1) Spin allowed
        mobj = re.findall(
            r'^\s+' + r'----------------------------------------------------------------------------' + r'\s*' +
            r'^\s+' + r'Root' + r'\s+' + r'(\d+)' + r'\s+' + r'(\w+)' + r'\s+' + r'(.*?)' + r'\s+' + NUMBER +
            r'\s+a\.u\.\s+' + NUMBER + r'\s+eV\s*' + r'^\s+' +
            r'----------------------------------------------------------------------------' + r'\s*' + r'^\s+' +
            r'Transition Moments' + r'\s+X\s+' + NUMBER + r'\s+Y\s+' + NUMBER + r'\s+Z\s+' + NUMBER + r'\s*' + r'^\s+'
            + r'Transition Moments' + r'\s+XX\s+' + NUMBER + r'\s+XY\s+' + NUMBER + r'\s+XZ\s+' + NUMBER + r'\s*' +
            r'^\s+' + r'Transition Moments' + r'\s+YY\s+' + NUMBER + r'\s+YZ\s+' + NUMBER + r'\s+ZZ\s+' + NUMBER +
            r'\s*$', outtext, re.MULTILINE)

        if mobj:
            print('matched TDDFT with transition moments')
            for mobj_list in mobj:
                #print (mobj_list)
                #print (mobj_list[0]) #root number
                #print (mobj_list[1]) #singlet or triplet?
                #print (mobj_list[2]) #symmetry

                #### temporary psivars ####
                psivar['TDDFT ROOT %s %s %s EXCITATION ENERGY' % (mobj_list[0], mobj_list[1],
                                                                  mobj_list[2])] = mobj_list[3]  # in a.u.
                psivar ['TDDFT ROOT %s %s %s EXCITED STATE ENERGY' %(mobj_list[0],mobj_list[1],mobj_list[2])] = \
                    psivar ['DFT TOTAL ENERGY'] + Decimal(mobj_list[3])
                psivar['TDDFT ROOT %s DIPOLE X' % (mobj_list[0])] = mobj_list[5]
                psivar['TDDFT ROOT %s DIPOLE Y' % (mobj_list[0])] = mobj_list[6]
                psivar['TDDFT ROOT %s DIPOLE Z' % (mobj_list[0])] = mobj_list[7]
                psivar['TDDFT ROOT %s QUADRUPOLE XX' % (mobj_list[0])] = mobj_list[8]
                psivar['TDDFT ROOT %s QUADRUPOLE XY' % (mobj_list[0])] = mobj_list[9]
                psivar['TDDFT ROOT %s QUADRUPOLE XZ' % (mobj_list[0])] = mobj_list[10]
                psivar['TDDFT ROOT %s QUADRUPOLE YY' % (mobj_list[0])] = mobj_list[11]
                psivar['TDDFT ROOT %s QUADRUPOLE YZ' % (mobj_list[0])] = mobj_list[12]
                psivar['TDDFT ROOT %s QUADRUPOLE ZZ' % (mobj_list[0])] = mobj_list[13]


#       2) Spin forbidden
        mobj = re.findall(
            r'^\s+' + r'----------------------------------------------------------------------------' + r'\s*' +
            r'^\s+' + r'Root' + r'\s+' + r'(\d+)' + r'\s+' + r'(\w+)' + r'\s+' + r'(.*?)' + r'\s+' + NUMBER +
            r'\s+a\.u\.\s+' + NUMBER + r'\s+eV\s*' + r'^\s+' +
            r'----------------------------------------------------------------------------' + r'\s*' + r'^\s+' +
            r'Transition Moments' + r'\s+Spin forbidden\s*$', outtext, re.MULTILINE)

        if mobj:
            print('matched TDDFT - spin forbidden')
            for mobj_list in mobj:
                #### temporary psivars ####
                psivar['TDDFT ROOT %s %s %s EXCITATION ENERGY' % (mobj_list[0], mobj_list[1],
                                                                  mobj_list[2])] = mobj_list[3]  # in a.u.
                psivar ['TDDFT ROOT %s %s %s EXCITED STATE ENERGY' %(mobj_list[0],mobj_list[1],mobj_list[2])] = \
                    psivar ['DFT TOTAL ENERGY'] + Decimal(mobj_list[3])
            if mobj:
                print('Non-variation initial energy')  #prints out energy, 5 counts

        #Process geometry
        # 1) CHARGE
        # Read charge from SCF module
        mobj = re.search(r'^\s+' + r'charge          =' + r'\s+' + NUMBER + r'\s*$', outtext,
                         re.MULTILINE | re.IGNORECASE)

        if mobj:
            print('matched charge')
            out_charge = int(float(mobj.group(1)))

        # Read charge from General information (not scf module)
        mobj = re.search(r'^\s+' + r'Charge           :' + r'\s+' + r'(-?\d+)' + r'\s*$', outtext,
                         re.MULTILINE | re.IGNORECASE)

        if mobj:
            print('matched charge')
            out_charge = int(float(mobj.group(1)))

        # 2) MULTIPLICITY
        # Read multiplicity from SCF module
        mobj = re.search(r'^\s+' + r'open shells     =' + r'\s+' + r'(\d+)' + r'\s*$', outtext,
                         re.MULTILINE | re.IGNORECASE)

        if mobj:
            print('matched multiplicity')
            out_mult = int(mobj.group(1)) + 1

        # Read multiplicity from SCF module through alpha, beta electrons
        mobj = re.search(
            r'^\s+' + r'alpha electrons =' + r'\s+' + r'(\d+)' + r'\s*' + r'^\s+' + r'beta  electrons =' + r'\s+' +
            r'(\d+)' + r'\s*$', outtext, re.MULTILINE | re.IGNORECASE)

        if mobj:
            print('matched multiplicity')
            out_mult = int(mobj.group(1)) - int(mobj.group(2)) + 1  #nopen + 1

        # Read multiplicity from General information (not scf module)
        mobj = re.search(r'^\s+' + r'Spin multiplicity:' + r'\s+' + r'(\d+)' + r'\s*$', outtext,
                         re.MULTILINE | re.IGNORECASE)

        if mobj:
            print('matched multiplicity')
            out_mult = int(mobj.group(1))

        #3) Initial geometry
        mobj = re.search(
            r'^\s+' + r'Geometry' + r'.*' + r'\s*' + r'^\s+' + r'(?:-+)\s*' + r'\s+' + r'\n' + r'^\s' +
            r'Output coordinates in ' + r'(.*?)' + r'\s' + r'\(scale by' + r'.*' + r'\s' + r'to convert to a\.u\.\)' +
            r'\s+' + r'\n' + r'^\s+' + r'No\.\       Tag          Charge          X              Y              Z' +
            r'\s*' + r'^\s+' + r'---- ---------------- ---------- -------------- -------------- --------------' +
            r'\s*' +
            r'((?:\s+([1-9][0-9]*)+\s+([A-Z][a-z]*)+\s+\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)'
            + r'\s*$', outtext, re.MULTILINE | re.IGNORECASE)

        if mobj:
            print('matched geom')

            #dinky molecule w/ charge and multiplicity
            if mobj.group(1) == 'angstroms':
                molxyz = '%d \n%d %d tag\n' % (len(mobj.group(2).splitlines()), out_charge, out_mult)  # unit = angstrom
                for line in mobj.group(2).splitlines():
                    lline = line.split()
                    molxyz += '%s %16s %16s %16s\n' % (lline[-5], lline[-3], lline[-2], lline[-1])
                    # Jiyoung was collecting charge (-4)? see if this is ok for ghosts
                    # Tag    ,    X,        Y,        Z
                psivar_coord = Molecule(validate=False, **qcel.molparse.to_schema(qcel.molparse.from_string(molxyz, dtype='xyz+', fix_com=True, fix_orientation=True)["qm"], dtype=2))

            else:  # unit = a.u.
                molxyz = '%d au\n%d %d tag\n' % (len(mobj.group(2).splitlines()), out_charge, out_mult)
                for line in mobj.group(2).splitlines():
                    lline = line.split()
                    molxyz += '%s %16s %16s %16s\n' % (int(float(lline[-4])), lline[-3], lline[-2], lline[-1])
                    # Tag    ,    X,        Y,        Z
                psivar_coord = Molecule(validate=False, **qcel.molparse.to_schema(qcel.molparse.from_string(molxyz, dtype='xyz+', fix_com=True, fix_orientation=True)["qm"], dtype=2))

        #Process gradient
        mobj = re.search(
            r'^\s+' + r'.*' + r'ENERGY GRADIENTS' + r'\s*' + r'\s+' + r'\n' + r'^\s+' +
            r'atom               coordinates                        gradient' + r'\s*' + r'^\s+' +
            r'x          y          z           x          y          z' + r'\s*' +
            r'((?:\s+([1-9][0-9]*)+\s+([A-Z][a-x]*)+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)'
            + r'\s*$', outtext, re.MULTILINE)

        if mobj:
            print('matched molgrad')
            atoms = []
            psivar_grad = []
            for line in mobj.group(1).splitlines():
                lline = line.split()  # Num, Tag, coord x, coord y, coord z, grad x, grad y, grad z
                #print (lline)
                if lline == []:
                    pass
                else:
                    atoms.append(lline[1])  #Tag
                    psivar_grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])

        #Process dipole
        mobj = re.search(
            r'^\s+' + r'Dipole moment' + r'\s+' + NUMBER + r'\s' + r'A\.U\.' + r'\s*' + r'^\s+' + r'DMX' + r'\s+' +
            NUMBER + r'.*' + r'\s*' + r'^\s+' + r'DMY' + r'\s+' + NUMBER + r'.*' + r'\s*' + r'^\s+' + r'DMZ' + r'\s+' +
            NUMBER + r'.' + r'\s*' + r'^\s+' + r'.*' + r'\s*' + r'^\s+' + r'Total dipole' + r'\s+' + NUMBER + r'\s' +
            r'A\.U\.' + r'\s*' + r'^\s+' + r'Dipole moment' + r'\s+' + NUMBER + r'\s' + r'Debye\(s\)' + r'\s*' +
            r'^\s+' + r'DMX' + r'\s+' + NUMBER + r'.*' + r'\s*' + r'^\s+' + r'DMY' + r'\s+' + NUMBER + r'.*' + r'\s*' +
            r'^\s+' + r'DMZ' + r'\s+' + NUMBER + r'.*' + r'\s*' + r'^\s+' + r'.*' + r'\s*' + r'^\s+' + r'Total dipole'
            + r'\s+' + NUMBER + r'\s' + r'DEBYE\(S\)' + r'\s*$', outtext, re.MULTILINE)

        if mobj:
            print('matched total dipole')
            #print (mobj.group(5), 'a.u.')
            #print (mobj.group(10),'debye(s)')

            #UNIT = DEBYE(S)
            psivar['CURRENT DIPOLE X'] = mobj.group(7)
            psivar['CURRENT DIPOLE Y'] = mobj.group(8)
            psivar['CURRENT DIPOLE Z'] = mobj.group(9)
            # total?

            #Process error code
            mobj = re.search(
                r'^\s+' + r'current input line \:' + r'\s*' + r'^\s+' + r'([1-9][0-9]*)' + r'\:' + r'\s+' + r'(.*)' +
                r'\s*' + r'^\s+'
                r'------------------------------------------------------------------------' + r'\s*' + r'^\s+'
                r'------------------------------------------------------------------------' + r'\s*' + r'^\s+' +
                r'There is an error in the input file' + r'\s*$', outtext, re.MULTILINE)
            if mobj:
                print('matched error')
            #print (mobj.group(1)) #error line number
        #print (mobj.group(2)) #error reason
            psivar['NWCHEM ERROR CODE'] = mobj.group(1)
            # TODO process errors into error var

    # Process CURRENT energies (TODO: needs better way)
    if 'HF TOTAL ENERGY' in psivar:
        psivar['SCF TOTAL ENERGY'] = psivar['HF TOTAL ENERGY']
        psivar['CURRENT REFERENCE ENERGY'] = psivar['HF TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['HF TOTAL ENERGY']

    if 'MP2 TOTAL ENERGY' in psivar and 'MP2 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP2 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP2 TOTAL ENERGY']
    if 'MP3 TOTAL ENERGY' in psivar and 'MP3 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP3 TOTAL ENERGY']
    if 'MP4 TOTAL ENERGY' in psivar and 'MP4 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP4 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP4 TOTAL ENERGY']
    
    if 'DFT TOTAL ENERGY' in psivar:
        psivar['CURRENT REFERENCE ENERGY'] = psivar['DFT TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['DFT TOTAL ENERGY']

    # Process TCE CURRENT energies
    # Need to be fixed
    # HOW TO KNOW options['NWCHEM']['NWCHEM_TCE']['value']?
    # TODO: CURRENT ENERGY = TCE ENERGY
    if ('%s TOTAL ENERGY' % (cc_name) in psivar and \
       ('%s CORRELATION ENERGY' % (cc_name) in psivar)):
        psivar['CURRENT CORRELATION ENERGY'] = psivar['%s CORRELATION ENERGY' % (cc_name)]
        psivar['CURRENT ENERGY'] = psivar['%s TOTAL ENERGY' % (cc_name)]

    if 'CCSD TOTAL ENERGY' in psivar and 'CCSD CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD TOTAL ENERGY']

    if 'CCSD(T) TOTAL ENERGY' in psivar and 'CCSD(T) CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD(T) CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD(T) TOTAL ENERGY']

    if ('EOM-%s TOTAL ENERGY' % (cc_name) in psivar) and \
       ('%s EXCITATION ENERGY' %(cc_name) in psivar):
        psivar['CURRENT ENERGY'] = psivar['EOM-%s TOTAL ENERGY' % (cc_name)]
        psivar['CURRENT EXCITATION ENERGY'] = psivar['%s EXCITATION ENERGY' % (cc_name)]

    return psivar, psivar_coord, psivar_grad, version, error


def harvest_hessian(hess):
    """Parses the contents *hessian* of the NWChem hess file into a hessian array.
    Hess file name has to be "nwchem.hess". (default)

    """
    raise NotImplementedError()


def harvest(p4Mol, nwout, **largs):  #check orientation and scratch files
    """Parses all the pieces of output from NWChem: the stdout in
    *nwout* Scratch files are not yet considered at this moment. 

    """

    outPsivar, outMol, outGrad, version, error = harvest_output(nwout)
    # print ('outMol', outMol)
    #print ('outGrad', outGrad)
    # print ('P4Mol', p4Mol)
    #Coordinate info check

    if outMol:
        if p4Mol:
            if abs(outMol.nuclear_repulsion_energy() - p4Mol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValidationError("""NWChem outfile (NRE: %f) inconsistent with Psi4 input (NRE: %f).""" % \
                            (outMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy()))
            ## TEST##########
            #else:
            #   print ( """NWChem outfile (NRE: %f) consistent with Psi4 input (NRE: %f).""" % \
            #             (outMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy()))
            ################
    else:
        raise ValidationError("""No coordinate information extracted from NWChem output.""")

    return outPsivar, None, outGrad, outMol, version, error


def nwchem_psivar_list():
    """Return a dict with keys of most NWChem methods and values of dicts
    with the PSI Variables returned by those methods. Used by cbs()
    wrapper to avoid unnecessary computations in compound methods.
    Result is appended to ``VARH``.
        
    """
    VARH = {}
    VARH['nwc-scf'] = {'nwc-scf': 'HF TOTAL ENERGY'}
    VARH['nwc-hf'] = {'nwc-hf': 'HF TOTAL ENERGY'}
    VARH['nwc-mp2'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY'}
    VARH['nwc-mp3'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY', 'nwc-mp3':'MP3 TOTAL ENERGY'}
    VARH['nwc-mp4'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY', 'nwc-mp3':'MP3 TOTAL ENERGY', 'nwc-mp4':'MP4 TOTAL ENERGY'}
    VARH['nwc-lccd'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-lccd': 'LCCD TOTAL ENERGY'}
    VARH['nwc-cisd'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY', 'nwc-cisd': 'CISD TOTAL ENERGY'}
    VARH['nwc-cisdt'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY', 'nwc-cisd': 'CISD TOTAL ENERGY', 'nwc-cisdt': 'CISDT TOTAL ENERGY'}
    VARH['nwc-cisdtq'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY','nwc-cisd': 'CISD TOTAL ENERGY', 'nwc-cisdt': 'CISDT TOTAL ENERGY', 'nwc-cisdtq': 'CISDTQ TOTAL ENERGY'}
    VARH['nwc-ccsd'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-mp2': 'MP2 TOTAL ENERGY', 'nwc-ccsd': 'CCSD TOTAL ENERGY'}

    VARH['nwc-ccsdt'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-ccsdt': 'CCSDT TOTAL ENERGY'}

    VARH['nwc-ccsdtq'] = {'nwc-hf': 'HF TOTAL ENERGY', 'nwc-ccsdtq': 'CCSDTQ TOTAL ENERGY'}

    VARH['nwc-ccsd(t)'] = {
        'nwc-hf': 'HF TOTAL ENERGY',
        'nwc-mp2': 'MP2 TOTAL ENERGY',
        'nwc-ccsd': 'CCSD TOTAL ENERGY',
        'nwc-ccsd(t)': 'CCSD(T) TOTAL ENERGY'
    }

    #VARH['nwc-dft'] = {
    #                        'nwc-dfttot' : 'DFT TOTAL ENERGY'}

    VARH['nwc-eom-ccsd'] = {
        'nwc-hf': 'HF TOTAL ENERGY',
        #'nwc-dfttot' : 'DFT TOTAL ENERGY',
        'nwc-ccsd': 'CCSD TOTAL ENERGY'
    }
    return VARH
