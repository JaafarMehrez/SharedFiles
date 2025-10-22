import pyscf
import numpy as np
import warnings
from functools import reduce
from pyscf import gto
from pyscf.gto.basis import parse_gaussian
from pyscf.mp.mp2 import get_frozen_mask
from pyscf import lib, ao2mo
from pyscf.scf import atom_hf
from pyscf import __config__

np.set_printoptions(threshold=np.inf)

def _normalize_spin_masks(mask, nmo):
    if isinstance(mask, (tuple, list)) and len(mask) == 2:
        ma = np.asarray(mask[0], dtype=bool)
        mb = np.asarray(mask[1], dtype=bool)
    else:
        m = np.asarray(mask, dtype=bool)
        if m.size == 2 * nmo:
            ma = m[:nmo].copy()
            mb = m[nmo:].copy()
        elif m.size == nmo:
            ma = mb = m.copy()
        else:
            raise ValueError(f"Mask length {m.size} not compatible with nmo={nmo}.")
    return ma, mb

def freezeCore(oneBody_a, oneBody_b, twoBody_a, twoBody_b, twoBody_ab, frozen_core, active=None):
    nmo_a = np.asarray(oneBody_a).shape[0]
    nmo_b = np.asarray(oneBody_b).shape[0]
    if nmo_a != nmo_b:
        raise ValueError(f"Different nmo in oneBody_a ({nmo_a}) and oneBody_b ({nmo_b}).")
    nmo = nmo_a
    core_a, core_b = _normalize_spin_masks(frozen_core, nmo)
    if active is None:
        valence_a = ~core_a
        valence_b = ~core_b
    else:
        valence_a, valence_b = _normalize_spin_masks(active, nmo)
    assert core_a.shape == (nmo,) and core_b.shape == (nmo,)
    assert valence_a.shape == (nmo,) and valence_b.shape == (nmo,)
    core_idx_a = np.where(core_a)[0]
    core_idx_b = np.where(core_b)[0]
    val_idx_a  = np.where(valence_a)[0]
    val_idx_b  = np.where(valence_b)[0]
    constant   = np.einsum('ii->', oneBody_a[np.ix_(core_idx_a, core_idx_a)])
    constant  += np.einsum('ii->', oneBody_b[np.ix_(core_idx_b, core_idx_b)])
    core_aa    = twoBody_a[np.ix_(core_idx_a, core_idx_a, core_idx_a, core_idx_a)]
    core_bb    = twoBody_b[np.ix_(core_idx_b, core_idx_b, core_idx_b, core_idx_b)]
    core_ab    = twoBody_ab[np.ix_(core_idx_a, core_idx_a, core_idx_b, core_idx_b)]
    constant  += 0.5 * ( np.einsum('iijj->', core_aa) - np.einsum('ijji->', core_aa)
              + np.einsum('iijj->', core_bb) - np.einsum('ijji->', core_bb)
              + 2.0 * np.einsum('iijj->', core_ab) )
    
    h_active_a = oneBody_a[np.ix_(val_idx_a, val_idx_a)].copy()
    coul_aa = twoBody_a[np.ix_(val_idx_a, val_idx_a, core_idx_a, core_idx_a)]
    coul_ab = twoBody_ab[np.ix_(val_idx_a, val_idx_a, core_idx_b, core_idx_b)]
    h_active_a += np.einsum('pqkk->pq', coul_aa)
    h_active_a += np.einsum('pqkk->pq', coul_ab)
    exch_aa = twoBody_a[np.ix_(val_idx_a, core_idx_a, core_idx_a, val_idx_a)]
    h_active_a -= np.einsum('pkkq->pq', exch_aa)
    
    h_active_b = oneBody_b[np.ix_(val_idx_b, val_idx_b)].copy()
    coul_bb = twoBody_b[np.ix_(val_idx_b, val_idx_b, core_idx_b, core_idx_b)]
    h_active_b += np.einsum('pqkk->pq', coul_bb)
    twoBody_ab_T = twoBody_ab.transpose(2,3,0,1)
    coul_ba = twoBody_ab_T[np.ix_(val_idx_b, val_idx_b, core_idx_a, core_idx_a)]
    h_active_b += np.einsum('pqkk->pq', coul_ba)
    exch_bb = twoBody_b[np.ix_(val_idx_b, core_idx_b, core_idx_b, val_idx_b)]
    h_active_b -= np.einsum('pkkq->pq', exch_bb)
    
    twoBody_active_a  = twoBody_a[np.ix_(val_idx_a, val_idx_a, val_idx_a, val_idx_a)].copy()
    twoBody_active_b  = twoBody_b[np.ix_(val_idx_b, val_idx_b, val_idx_b, val_idx_b)].copy()
    twoBody_active_ab = twoBody_ab[np.ix_(val_idx_a, val_idx_a, val_idx_b, val_idx_b)].copy()
    
    return h_active_a, h_active_b, twoBody_active_a, twoBody_active_b, twoBody_active_ab, constant

def main():
    from pyscf import gto, scf, ao2mo
    name = 'out'
    mol = pyscf.M(
        atom='O',
        unit='angstrom',
        basis={'O': parse_gaussian.load('O-STO-3G-EMSL.gbs', 'O')},
        charge=0,
        spin=2,
        verbose=9,
        symmetry=True,
        output=name + '.txt',
        symmetry_subgroup='D2h',
        max_memory=4000,
    )

    original_AtomSphAverageRHF = atom_hf.AtomSphAverageRHF

    class CustomAtomSphAverageRHF(original_AtomSphAverageRHF):
        def __init__(self, mol):
            super().__init__(mol)
            self.max_cycle = 9999
            self.direct_scf = False

    atom_hf.AtomSphAverageRHF = CustomAtomSphAverageRHF
    mymf = mol.UHF().set(
        conv_tol=1e-14,
        max_cycle=9999,
        ddm_tol=1e-15,
        direct_scf=False,
        chkfile=name + '.chk',
        init_guess='atom',
        irrep_nelec={'Ag': 4, 'B3u': 2, 'B2u': 1, 'B1u': 1}
    )
    mymf.kernel()
    atom_hf.AtomSphAverageRHF = original_AtomSphAverageRHF
    
    def compute_mo_irreps(mol, mo_coeff):
        symm_orbs = mol.symm_orb
        irrep_labels = mol.irrep_name
        mo_irreps = []
        for mo in mo_coeff.T:
            projections = [np.linalg.norm(symm_orbs[i].T @ mo) for i in range(len(symm_orbs))]
            irrep_idx = np.argmax(projections)
            mo_irreps.append(irrep_labels[irrep_idx])
        return mo_irreps

    def align_beta_orbitals_to_alpha(mol, mo_coeff):
        alpha_orbs, beta_orbs = mo_coeff[0], mo_coeff[1]
        alpha_irreps = compute_mo_irreps(mol, mo_coeff[0])
        beta_irreps = compute_mo_irreps(mol, mo_coeff[1])
        beta_orbs_sorted = []
        used_indices = set()
        for target_irrep in alpha_irreps:
            for idx, beta_irrep in enumerate(beta_irreps):
                if beta_irrep == target_irrep and idx not in used_indices:
                    beta_orbs_sorted.append(beta_orbs[:, idx])
                    used_indices.add(idx)
                    break
            else:
                raise ValueError(f"No matching beta orbital found for alpha irrep: {target_irrep}")
        beta_orbs_sorted = np.column_stack(beta_orbs_sorted)
        return alpha_orbs, beta_orbs_sorted
    
    
    mo_coeff = mymf.mo_coeff
    mol = mymf.mol
    
    alpha_irreps = compute_mo_irreps(mol, mo_coeff[0])
    beta_irreps = compute_mo_irreps(mol, mo_coeff[1])
    print(alpha_irreps)
    print(beta_irreps)
    
    assert mo_coeff[0].dtype == np.double and mo_coeff[1].dtype == np.double
    
    mo_coeff_a, mo_coeff_b = align_beta_orbitals_to_alpha(mol, mo_coeff)
    
    mymf.mo_coeff = (mo_coeff_a, mo_coeff_b)
    
    from pyscf import cc
    mycc = cc.UCCSD(mymf, frozen=1)
        
    orbsym_full = getattr(mo_coeff_a, 'orbsym', None)
    nuc = mymf.energy_nuc()
    
    nmo = mo_coeff_a.shape[0]
    
    h1e_a  = reduce(np.dot, (mo_coeff_a.T, mymf.get_hcore(), mo_coeff_a))
    h1e_b  = reduce(np.dot, (mo_coeff_b.T, mymf.get_hcore(), mo_coeff_b))
    
    eri_a  = ao2mo.restore(1,ao2mo.incore.general(mymf._eri,(mo_coeff_a, mo_coeff_a, mo_coeff_a, mo_coeff_a),compact=False),nmo)
    eri_b  = ao2mo.restore(1,ao2mo.incore.general(mymf._eri,(mo_coeff_b, mo_coeff_b, mo_coeff_b, mo_coeff_b),compact=False),nmo)
    eri_ab = ao2mo.restore(1,ao2mo.incore.general(mymf._eri,(mo_coeff_a, mo_coeff_a, mo_coeff_b, mo_coeff_b),compact=False),nmo)

    active = get_frozen_mask(mycc)
    active_in = active

    mo_occ_obj = getattr(getattr(mycc, '_scf', None) or getattr(mycc, 'mf', None) or mycc, 'mo_occ', None)
    if isinstance(mo_occ_obj, np.ndarray) and mo_occ_obj.ndim == 2 and mo_occ_obj.shape[0] == 2:
        mo_occ_a = mo_occ_obj[0]
        mo_occ_b = mo_occ_obj[1]
    elif isinstance(mo_occ_obj, (list, tuple)) and len(mo_occ_obj) == 2:
        mo_occ_a = np.asarray(mo_occ_obj[0])
        mo_occ_b = np.asarray(mo_occ_obj[1])
    else:
        mo_occ_a = mo_occ_b = np.asarray(mo_occ_obj)

    nmo = mo_occ_a.size

    if isinstance(active_in, (tuple, list)) and len(active_in) == 2:
        act_a = np.asarray(active_in[0], dtype=bool)
        act_b = np.asarray(active_in[1], dtype=bool)
    else:
        act = np.asarray(active_in, dtype=bool)
        if act.size == 2 * nmo:
            act_a = act[:nmo].copy()
            act_b = act[nmo:].copy()
        elif act.size == nmo:
            act_a = act_b = act.copy()
        else:
            raise ValueError(f"Unexpected active length {act.size}; expected {nmo} or {2*nmo}")
        
    shared = np.asarray(act_a, dtype=bool) & np.asarray(act_b, dtype=bool)
    
    act_a[:] = shared
    act_b[:] = shared
    
    active_full = np.empty(2 * nmo, dtype=bool)
    active_full[:nmo] = act_a
    active_full[nmo:] = act_b
    active = active_full
    frozen_core = np.zeros_like(active, dtype=np.bool_)
    nocc_full = mol.nelectron // 2
    frozen_core[:nocc_full] = ~active[:nocc_full]

    nocc_a = int(np.count_nonzero(np.asarray(mo_occ_a) > 0.5))
    nocc_b = int(np.count_nonzero(np.asarray(mo_occ_b) > 0.5))

    frozen_core_a = np.zeros_like(act_a, dtype=bool)
    frozen_core_b = np.zeros_like(act_b, dtype=bool)

    frozen_core_a[:nocc_a] = ~act_a[:nocc_a]
    frozen_core_b[:nocc_b] = ~act_b[:nocc_b]

    h1e_a, h1e_b, eri_a, eri_b, eri_ab, constant = freezeCore(h1e_a, h1e_b, eri_a, eri_b, eri_ab, frozen_core=(frozen_core_a, frozen_core_b), active=(act_a, act_b))
    
    nmo_active = h1e_a.shape[0]

    if orbsym_full is None:
        orbsym_active = None
    else:
        orbsym_arr = np.asarray(orbsym_full, dtype=int)
        valence_idx = np.where(act_a)[0]
        orbsym_active = [int(x) for x in orbsym_arr[valence_idx]]
        
    if orbsym_active is not None:
        if len(orbsym_active) != nmo_active:
            raise RuntimeError(f"Length mismatch: orbsym_active has length {len(orbsym_active)} "
                               f"but number of active orbitals is {nmo_active}.")

    filename = 'fort.55'
    nelec_full = mol.nelectron
    spin = mol.spin       
    nalpha_full = (nelec_full + spin) // 2
    nbeta_full  = (nelec_full - spin) // 2
    nocc_full = mol.nelectron // 2
    nfrozen_core = int(np.count_nonzero(frozen_core[:nocc_full]))
    nalpha_active = nalpha_full - nfrozen_core
    nbeta_active  = nbeta_full  - nfrozen_core
    nelec_active = (nalpha_active, nbeta_active)
    
    
    DEFAULT_FLOAT_FORMAT = getattr(__config__, 'fcidump_float_format', ' %.16g')
    TOL = getattr(__config__, 'fcidump_write_tol', 1e-15)
    
    def write_hcore_uhf(fout, h1e_a, h1e_b, nmo, tol=TOL, float_format=DEFAULT_FLOAT_FORMAT):
        h1e_a = h1e_a.reshape(nmo,nmo)
        h1e_b = h1e_b.reshape(nmo,nmo)
        indx = [i+1 for i in range(nmo)]
        output_format = float_format + ' %5d %5d     0     0\n'
        for i in range(nmo):
            for j in range(i, nmo):
                if abs(h1e_a[i,j]) > TOL:
                    fout.write(output_format % (h1e_a[i,j], indx[i], indx[j]))
        fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')
        for i in range(nmo):
            for j in range(i, nmo):
                if abs(h1e_b[i,j]) > TOL:
                    fout.write(output_format % (h1e_b[i,j], indx[i], indx[j]))
        fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')
    
    print("DEBUG: eri shapes / dtypes")
    print("eri_a.shape =", np.asarray(eri_a).shape, "dtype=", np.asarray(eri_a).dtype)
    print("eri_b.shape =", np.asarray(eri_b).shape, "dtype=", np.asarray(eri_b).dtype)
    print("eri_ab.shape=", np.asarray(eri_ab).shape, "dtype=", np.asarray(eri_ab).dtype)
    print("sample eri_a[0,0,0,0] =", np.asarray(eri_a).flatten()[0] if np.asarray(eri_a).size>0 else None)
    
    
    def write_eri_uhf(fout, eri_a, eri_b, eri_ab, nmo, tol=TOL, float_format=DEFAULT_FLOAT_FORMAT):
        eri_a = np.asarray(eri_a)
        eri_b = np.asarray(eri_b)
        eri_ab = np.asarray(eri_ab)
        npair = nmo * (nmo + 1) // 2
        output_format = float_format + ' %5d %5d %5d %5d\n'
        indx = [i + 1 for i in range(nmo)]
        
        def pair_index(i, j):
            return i * (i + 1) // 2 + j

        if eri_a.ndim == 2 and eri_b.ndim == 2 and eri_ab.ndim == 2:
            assert eri_a.shape == (npair, npair) and eri_b.shape == (npair, npair) and eri_ab.shape == (npair, npair)
            kl = 0
            for l in range(nmo):
                for k in range(0, l+1):
                    ij = 0
                    for i in range(0, nmo):
                        for j in range(0, i+1):
                            if i >= k:
                                if abs(eri_a[ij, kl]) > tol:
                                    fout.write(output_format % (eri_a[ij, kl], indx[i], indx[j], indx[k], indx[l]))
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
                                if abs(eri_b[ij, kl]) > tol:
                                    fout.write(output_format % (eri_b[ij, kl], indx[i], indx[j], indx[k], indx[l]))
                            ij += 1
                    kl += 1
            fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')

            ij = 0
            for j in range(nmo):
                for i in range(0, j+1):
                    kl = 0
                    for k in range(nmo):
                        for l in range(0, k+1):
                            if abs(eri_ab[ij, kl]) > tol:
                                fout.write(output_format % (eri_ab[ij, kl], indx[i], indx[j], indx[k], indx[l]))
                            kl += 1
                    ij += 1
            fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')
            return

        # CASE B: full 4D arrays
        if eri_a.ndim == 4 and eri_b.ndim == 4 and eri_ab.ndim == 4:
            for i in range(nmo):
                for j in range(0, i + 1):
                    ij_idx = pair_index(i, j)
                    for k in range(nmo):
                        for l in range(0, k + 1):
                            kl_idx = pair_index(k, l)
                            if ij_idx >= kl_idx:
                                val = eri_a[i, j, k, l]
                                if abs(val) > tol:
                                    fout.write(output_format % (val, indx[i], indx[j], indx[k], indx[l]))
            fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')

            for i in range(nmo):
                for j in range(0, i + 1):
                    ij_idx = pair_index(i, j)
                    for k in range(nmo):
                        for l in range(0, k + 1):
                            kl_idx = pair_index(k, l)
                            if ij_idx >= kl_idx:
                                val = eri_b[i, j, k, l]
                                if abs(val) > tol:
                                    fout.write(output_format % (val, indx[i], indx[j], indx[k], indx[l]))
            fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')

            for i in range(nmo):
                for j in range(0, i + 1):
                    for k in range(nmo):
                        for l in range(0, k + 1):
                            val = eri_ab[i, j, k, l]
                            if abs(val) > tol:
                                fout.write(output_format % (val, indx[i], indx[j], indx[k], indx[l]))
            fout.write(' 0.00000000000000000000E+00' + '     0     0     0     0\n')
            return

        raise RuntimeError(f"Unsupported ERI shapes: eri_a {eri_a.shape}, eri_b {eri_b.shape}, eri_ab {eri_ab.shape}")
        
    def write_head(fout, nmo, nelec, ms=0, orbsym=None):
        is_uhf = isinstance(nelec, (list, tuple)) and len(nelec) == 2 and nelec[0] != nelec[1]
        if not isinstance(nelec, (int, np.number)):
            ms = abs(nelec[0] - nelec[1])
            nelec = nelec[0] + nelec[1]
        fout.write(' &FCI NORB=%4d,NELEC=%2d,MS2=%d,\n' % (nmo, nelec, ms))
        if orbsym is not None and len(orbsym) > 0:
            fout.write('  ORBSYM=%s\n' % ','.join([str(x) for x in orbsym]))
        else:
            fout.write('  ORBSYM=%s\n' % ('1,' * nmo))
        fout.write('  ISYM=1,\n')
        if is_uhf:
            fout.write('  IUHF=1,\n')
        fout.write(' &END\n')
    
    
    def write_head55(fout, nmo, nelec, ms=0, orbsym=None):
        if not isinstance(nelec, (int, np.number)):
            ms = abs(nelec[0] - nelec[1])
            nelec = nelec[0] + nelec[1]
        fout.write(f"{nmo:1d} {nelec:1d}\n")
        if orbsym is not None and len(orbsym) > 0:
            orbsym = [x + 1 for x in orbsym]
            fout.write(f"{' '.join([str(x) for x in orbsym])}\n")
        else:
            fout.write(f"{' 1' * nmo}\n")
        fout.write(' 150000\n')
    
    
    def from_integrals_uhf(filename, h1e_a, h1e_b, eri_a, eri_b, eri_ab, nmo, nelec, nuc=0, ms=0, orbsym=None,
                       tol=TOL, float_format=DEFAULT_FLOAT_FORMAT):
        with open(filename, 'w') as fout:
            if filename == 'fort.55':
                write_head55(fout, nmo, nelec, ms, orbsym)
            else:
                write_head(fout, nmo, nelec, ms, orbsym)
            write_eri_uhf(fout, eri_a, eri_b, eri_ab, nmo, tol=tol, float_format=float_format)
            write_hcore_uhf(fout, h1e_a, h1e_b, nmo, tol=tol, float_format=float_format)
            output_format = float_format + '     0     0     0     0\n'
            fout.write(output_format % nuc)
    from_integrals_uhf(filename, h1e_a, h1e_b, eri_a, eri_b, eri_ab, nmo_active, nelec_active, nuc + constant, 0, orbsym_active, tol=1e-18, float_format='% 0.20E')
    
if __name__ == '__main__':
    main()

