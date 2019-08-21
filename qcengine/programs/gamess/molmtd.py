
def muster_modelchem(mtd, dertype):
    ropts = {}

    runtyp = {0: 'energy',
              1: 'gradient',
              2: 'hessian',
              #'properties': 'prop',
             }[dertype]

    ropts['contrl__runtyp'] = runtyp

    if mtd == 'gms-gamess':
        pass

    #elif mtd in ['scf', 'hf']:
    #    ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
    #    ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
    #    ropts.require('GAMESS', 'contrl__cctyp', 'none', accession=accession, verbose=verbose)

    elif mtd == 'mp2':
        ropts['contrl__mplevl'] = 2
        #ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        #ropts.require('GAMESS', 'contrl__cctyp', 'none', accession=accession, verbose=verbose)

    elif mtd == 'ccsd':
        #ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
        #ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
        ropts['contrl__cctyp'] = 'ccsd'

    #elif mtd == 'gms-ccsd(t)':
    #    ropts.require('GAMESS', 'contrl__mplevl', 0, accession=accession, verbose=verbose)
    #    ropts.require('GAMESS', 'contrl__cityp', 'none', accession=accession, verbose=verbose)
    #    ropts.require('GAMESS', 'contrl__cctyp', 'ccsd(t)', accession=accession, verbose=verbose)

    return ropts
