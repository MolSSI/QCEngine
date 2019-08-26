def muster_modelchem(method, derint):
    method = method.lower()
    opts = {}

    runtyp = {
        0: 'energy',
        1: 'gradient',
        2: 'hessian',
        #'properties': 'prop',
    }[derint]

    opts['contrl__runtyp'] = runtyp

    if method == 'gamess':
        pass

    elif method in ['scf', 'hf']:
        pass

        # opts['contrl__mplevl'] = 0
        # opts['contrl__cityp'] = 'none'
        # opts['contrl__cctyp'] = 'none'

    elif method == 'mp2':
        opts['contrl__mplevl'] = 2

    elif method == 'ccsd':
        opts['contrl__cctyp'] = 'ccsd'

    elif method == 'ccsd(t)':
        opts['contrl__cctyp'] = 'ccsd(t)'

    return opts
