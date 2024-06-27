import textwrap
import sys
import glob
import os
import subprocess
import pyscf
import numpy
from functools import reduce
from pyscf import cc, lib, tools, scf, symm, ao2mo
from pyscf.tools.fcidump import from_mo
from pyscf.tools.fcidump import from_integrals
from pyscf.gto.basis import parse_gaussian
from io import StringIO

name = 'out'
mol = pyscf.M(
    atom = '''
        O
    ''',
    unit = 'angstrom',
    basis = {
            'O' : parse_gaussian.load('STO-3G.gbs', 'O')
    },
    charge = 2,
    spin = 2,
    verbose = 9,
    symmetry = True,
    output = name +'.txt',
    symmetry_subgroup = 'D2h',
    max_memory = 4000,
)
mol.max_memory =4000
mf = mol.UHF().set(conv_tol=1e-10,max_cycle=999,direct_scf_tol=1e-14,chkfile=name+'.chk',init_guess='atom',irrep_nelec={'Ag': 4, 'B3u':1 , 'B2u':1 ,'B1u':0 })
mf.kernel()
#pyscf.tools.fcidump.from_chkfile('fcidump', name+'.chk', tol=1e-18, float_format=' %.16g')
mol, scf_rec = scf.chkfile.load_scf(name+'.chk')
mo_coeff = numpy.array(scf_rec['mo_coeff'])
orbsym_alpha = symm.label_orb_symm(mol, mol.irrep_id,mol.symm_orb, mo_coeff[0], check=False)
orbsym_beta = symm.label_orb_symm(mol, mol.irrep_id,mol.symm_orb, mo_coeff[1], check=False)
h1ao = scf.hf.get_hcore(mol)
h1e_alpha = reduce(numpy.dot, (mo_coeff[0].T, h1ao, mo_coeff[0]))
h1e_beta = reduce(numpy.dot, (mo_coeff[1].T, h1ao, mo_coeff[1]))
eri_alpha = ao2mo.full(mol, mo_coeff[0], verbose=0)
eri_beta = ao2mo.full(mol, mo_coeff[1], verbose=0)
nuc = mol.energy_nuc()
ms=mol.spin
from_integrals('fcidump_alpha', h1e_alpha, eri_alpha, h1e_alpha.shape[0], mol.nelec, nuc, ms, orbsym_alpha, tol=1e-18, float_format=' %.16g')
from_integrals('fcidump_beta', h1e_beta, eri_beta, h1e_beta.shape[0], mol.nelec, nuc, ms, orbsym_beta, tol=1e-18, float_format=' %.16g')

