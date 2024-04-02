#!/usr/bin/env python
# Author: Nike Dattani, nike@hpqc.org

import pyscf
from pyscf import cc, lib, tools
from pyscf.gto.basis import parse_gaussian

name = 'out'
mol = pyscf.M(
    atom = '''       
        O
    ''',
    unit = 'angstrom',
    basis = {
            'O' : parse_gaussian.load('/home/jaafar1/projects/def-nike-ab/jaafar1/PySCF/basis/aVDZ-EMSL.gbs', 'O')
    },
    charge = 0,
    spin = 2,
    verbose = 9,
    symmetry = True,
    output = name +'.txt',
    symmetry_subgroup = 'Dooh',
)
mol.max_memory =4000
mf = mol.ROHF().set(conv_tol=1e-10,max_cycle=999,direct_scf_tol=1e-14,chkfile=name+'.chk',init_guess='atom',irrep_nelec={'A1g': (2,2), 'E1ux': (1,0), 'E1uy': (1,0),'A1u': (1,1)})
mf.kernel()
pyscf.tools.fcidump.from_chkfile('out.fcidump', name+'.chk', tol=1e-18, float_format=' %.16g', molpro_orbsym=False, orbsym=None)
