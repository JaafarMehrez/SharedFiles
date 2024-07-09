import os, sys, glob, subprocess, textwrap, numpy
import pyscf
from functools import reduce
from pyscf import cc, lib, tools, scf, symm, ao2mo
from pyscf.tools.fcidump import from_mo
from pyscf.tools.fcidump import from_integrals
from pyscf.gto.basis import parse_gaussian
import pyscf.symm.param as param
import pyscf.lib.logger as logger
from subprocess import call
from io import StringIO
name = 'out'
mol = pyscf.M(
    atom = '''
        HE
    ''',
    unit = 'angstrom',
    basis = {
            'HE' : parse_gaussian.load('STO-3G.gbs', 'HE')
    },
    charge = 0,
    spin = 0,
    verbose = 9,
    symmetry = True,
    output = name +'.txt',
    symmetry_subgroup = 'D2h',
    max_memory = 4000,
)
mol.max_memory =4000
mf = mol.UHF().set(conv_tol=1e-10,max_cycle=999,direct_scf_tol=1e-14,chkfile=name+'.chk',init_guess='atom',irrep_nelec={'Ag': 2})
mf.kernel()
orbs = mf.mo_coeff
''' nmo is number of MO orbitals per spin channel
        note that ordering is abababa...   '''
nmo = orbs[0].shape[1]
eri_aaaa = pyscf.ao2mo.restore(8,pyscf.ao2mo.incore.general(mf._eri, (orbs[0],orbs[0],orbs[0],orbs[0]), compact=False),nmo)
eri_bbbb = pyscf.ao2mo.restore(8,pyscf.ao2mo.incore.general(mf._eri, (orbs[1],orbs[1],orbs[1],orbs[1]), compact=False),nmo)
eri_aabb = pyscf.ao2mo.restore(8,pyscf.ao2mo.incore.general(mf._eri, (orbs[0],orbs[0],orbs[1],orbs[1]), compact=False),nmo)
eri_bbaa = pyscf.ao2mo.restore(8,pyscf.ao2mo.incore.general(mf._eri, (orbs[1],orbs[1],orbs[0],orbs[0]), compact=False),nmo)
h_core = mf.get_hcore(mol)
h_aa = reduce(numpy.dot, (orbs[0].T, h_core, orbs[0]))
h_bb = reduce(numpy.dot, (orbs[1].T, h_core, orbs[1]))
nuc = mol.energy_nuc()
float_format = ' %.16g'
if mol.symmetry:
        groupname = mol.groupname
        if groupname in ('SO3', 'Dooh'):
            logger.info(mol, 'Lower symmetry from %s to D2h', groupname)
            raise RuntimeError('Lower symmetry from %s to D2h' % groupname)
        elif groupname == 'Coov':
            logger.info(mol, 'Lower symmetry from Coov to C2v')
            raise RuntimeError('''Lower symmetry from Coov to C2v''')
mol.orbsym = pyscf.symm.label_orb_symm(mol,mol.irrep_name,mol.symm_orb,orbs[0])
tmp_orblist = mol.orbsym.tolist()
tmp_orblist += pyscf.symm.label_orb_symm(mol,mol.irrep_name,mol.symm_orb,orbs[1]).tolist()
mol.orbsym = numpy.array(tmp_orblist)
orbsym = [param.IRREP_ID_TABLE[groupname][i]+1 for i in mol.orbsym]
#NECI wants its orbitals as a,b,a,b,a,b rather than aaaabbbb
assert(len(orbsym) % 2 == 0)
orbsym_reorder = [i for tup in zip(orbsym[:int(len(orbsym)/2)], orbsym[int(len(orbsym)/2):]) for i in tup] 
a_inds = [i*2+1 for i in range(orbs[0].shape[1])]
b_inds = [i*2+2 for i in range(orbs[1].shape[1])]
nelec = mol.nelec
tol=1e-15
with open('FCIDUMP', 'w') as fout:
        if not isinstance(nelec, (int, numpy.number)):
            ms = abs(nelec[0] - nelec[1])
            nelec = nelec[0] + nelec[1]
        else: ms=0
        fout.write(' &FCI NORB=%4d,NELEC=%2d,MS2=%d,\n' % (nmo*2, nelec, ms))
        if orbsym is not None and len(orbsym_reorder) > 0:
            fout.write('  ORBSYM=%s\n' % ','.join([str(x) for x in orbsym_reorder]))
        else:
            fout.write('  ORBSYM=%s\n' % ('1,' * 2*nmo))
        fout.write('  ISYM=1, UHF=TRUE\n')
        fout.write(' &END\n')
        # Assume 8-fold symmetry
        npair = nmo*(nmo+1)//2
        output_format = float_format + ' %4d %4d %4d %4d\n'
        ij = 0
        ijkl = 0
        for i in range(nmo):
            for j in range(0, i+1):
                kl = 0
                for k in range(0, i+1):
                    for l in range(0, k+1):
                        if ij >= kl:
                            if abs(eri_aaaa[ijkl]) > tol:
                                fout.write(output_format % (eri_aaaa[ijkl], a_inds[i], a_inds[j], a_inds[k], a_inds[l]))
                            if abs(eri_bbbb[ijkl]) > tol:
                                fout.write(output_format % (eri_bbbb[ijkl], b_inds[i], b_inds[j], b_inds[k], b_inds[l]))
                            if abs(eri_aabb[ijkl]) > tol:
                                fout.write(output_format % (eri_aabb[ijkl], a_inds[i], a_inds[j], b_inds[k], b_inds[l]))
                            if abs(eri_bbaa[ijkl]) > tol:
                                fout.write(output_format % (eri_bbaa[ijkl], b_inds[i], b_inds[j], a_inds[k], a_inds[l]))
                            ijkl += 1
                        kl += 1
                ij += 1
        h_aa = h_aa.reshape(nmo,nmo)
        h_bb = h_bb.reshape(nmo,nmo)
        output_format = float_format + ' %4d %4d  0  0\n'
        for i in range(nmo):
            for j in range(0, i+1):
                if abs(h_aa[i,j]) > tol:
                    fout.write(output_format % (h_aa[i,j], a_inds[i], a_inds[j]))
                if abs(h_bb[i,j]) > tol:
                    fout.write(output_format % (h_bb[i,j], b_inds[i], b_inds[j]))
        output_format = float_format + '  0  0  0  0\n'
        fout.write(output_format % nuc)
