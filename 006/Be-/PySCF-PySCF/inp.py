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
        BE
    ''',
    unit = 'angstrom',
    basis = {
            'BE' : parse_gaussian.load('BE-STO-XG.gbs', 'BE')
    },
    charge = -1,
    spin = 1,
    symmetry = True,
    verbose = 9,
    symmetry_subgroup = 'D2h',
    output = name +'.txt',
    max_memory = 4000,
)
mf = mol.UHF().set(conv_tol=1e-10,max_cycle=999,ddm_tol=1e-14,direct_scf_tol=1e-13,chkfile=name+'.chk',init_guess='atom',irrep_nelec={'Ag':4,'B3u':1})
mf.kernel()
mycc = cc.UCCSD(mf).set(conv_tol=1e-10).run()
et = mycc.ccsd_t()
print('CCSD total energy', mycc.e_tot)
print('CCSD(T) total energy', mycc.e_tot + et)

orbs = mf.mo_coeff
nmo = orbs[0].shape[0]
eri_aaaa = pyscf.ao2mo.restore(4,pyscf.ao2mo.incore.general(mf._eri, (orbs[0],orbs[0],orbs[0],orbs[0]), compact=False),nmo)
eri_bbbb = pyscf.ao2mo.restore(4,pyscf.ao2mo.incore.general(mf._eri, (orbs[1],orbs[1],orbs[1],orbs[1]), compact=False),nmo)
eri_aabb = pyscf.ao2mo.restore(4,pyscf.ao2mo.incore.general(mf._eri, (orbs[0],orbs[0],orbs[1],orbs[1]), compact=False),nmo)
h_core = mf.get_hcore(mol)
h_aa = reduce(numpy.dot, (orbs[0].T, h_core, orbs[0]))
h_bb = reduce(numpy.dot, (orbs[1].T, h_core, orbs[1]))
nuc = mol.energy_nuc()
float_format = '% 0.20E'
if mol.symmetry:
        groupname = mol.groupname
        if groupname in ('SO3', 'Dooh'):
            logger.info(mol, 'Lower symmetry from %s to D2h', groupname)
            raise RuntimeError('Lower symmetry from %s to D2h' % groupname)
        elif groupname == 'Coov':
            logger.info(mol, 'Lower symmetry from Coov to C2v')
            raise RuntimeError('''Lower symmetry from Coov to C2v''')
orbsym = pyscf.symm.label_orb_symm(mol,mol.irrep_name,mol.symm_orb,orbs[0])
orbsym = numpy.array(orbsym)
orbsym = [param.IRREP_ID_TABLE[groupname][i]+1 for i in orbsym]
a_inds = [i+1 for i in range(orbs[0].shape[0])]
b_inds = [i+1 for i in range(orbs[1].shape[1])]
nelec = mol.nelec
tol=1e-18
with open('fort.55', 'w') as fout:
        if not isinstance(nelec, (int, numpy.number)):
            ms    = abs(nelec[0] - nelec[1])
            nelec =     nelec[0]  + nelec[1]
        else: ms=0
        fout.write(f"{nmo:1d} {nelec:1d}\n")
        if orbsym is not None and len(orbsym) > 0:
            fout.write(f"{' '.join([str(x) for x in orbsym])}\n")
        else:
            fout.write(f"{' 1' * nmo}\n")
        fout.write(' 150000\n')
        output_format = float_format + ' %5d %5d %5d %5d\n'
        #4-fold symmetry
        kl = 0
        for l in range(nmo):
            for k in range(0, l+1):
                ij = 0
                for i in range(0, nmo):
                    for j in range(0, i+1):
                        if i >= k:
                            if abs(eri_aaaa[ij,kl]) > tol:
                                fout.write(output_format % (eri_aaaa[ij,kl], a_inds[i], a_inds[j], a_inds[k], a_inds[l]))
                        ij += 1
                kl += 1
        fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')
        kl = 0
        for l in range(nmo):
            for k in range(0, l+1):
                ij = 0
                for i in range(0, nmo):
                    for j in range(0, i+1):
                        if i >= k:
                            if abs(eri_bbbb[ij,kl]) > tol:
                                fout.write(output_format % (eri_bbbb[ij,kl], b_inds[i], b_inds[j], b_inds[k], b_inds[l]))
                        ij += 1
                kl += 1
        fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')
        ij = 0
        for j in range(nmo):
            for i in range(0, j+1):
                kl = 0
                for k in range(nmo):
                    for l in range(0, k+1):
                        if abs(eri_aabb[ij,kl]) > tol:
                            fout.write(output_format % (eri_aabb[ij,kl], a_inds[i], a_inds[j], b_inds[k], b_inds[l]))
                        kl += 1
                ij +=1

        fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')
        h_aa = h_aa.reshape(nmo,nmo)
        h_bb = h_bb.reshape(nmo,nmo)
        output_format = float_format + ' %5d %5d     0     0\n'
        for i in range(nmo):
            for j in range(nmo):
                fout.write(output_format % (h_aa[i,j], a_inds[i], a_inds[j]))
        fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')
        for i in range(nmo):
            for j in range(nmo):
                fout.write(output_format % (h_bb[i,j], b_inds[i], b_inds[j]))
        fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')
        output_format = float_format + '     0     0     0     0\n'
        fout.write(output_format % nuc)
