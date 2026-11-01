"""Full-source particle-hole identity and certified centered cube sectors."""
from functools import lru_cache
from fractions import Fraction as F
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from ..cube_eta_protection import cube as source
from ..interacting_vacuum import ground, global_ground


@lru_cache(maxsize=1)
def source_symmetry():
    ops=source.operators()
    if ops.scale!=16 or ops.hubbard_H.imag.nnz:
        raise ArithmeticError('expected the exact hidden-frame numerator16H')
    H=ops.hubbard_H.real.tocsr()
    masks=np.arange(1<<16,dtype=np.int64)
    number=np.array([int(m).bit_count() for m in masks],dtype=np.int64)
    complement=((1<<16)-1)^masks
    parity_mask=sum(1<<(2*x+spin) for x,site in enumerate(source.GEOMETRY.sites)
                    for spin in (0,1) if (sum(site)+spin)%2)
    phase=np.array([(-1)**((int(m)&parity_mask).bit_count()) for m in masks],dtype=np.int64)
    # U|m>=phase[m]|complement(m)>; full sixteen-mode U is its own adjoint.
    if not np.array_equal(phase[complement]*phase,np.ones(len(masks),dtype=np.int64)):
        raise ArithmeticError('particle-hole map is not the claimed unitary involution')
    transformed=H[complement][:,complement].multiply(phase[:,None]).multiply(phase[None,:]).tocsr()
    expected=H+sparse.diags(16*(8-number),format='csr',dtype=np.int64)
    if (transformed-expected).nnz:
        raise ArithmeticError('full-source particle-hole Hamiltonian identity failed')
    K=H-sparse.diags(8*number,format='csr',dtype=np.int64)
    transformed_K=K[complement][:,complement].multiply(phase[:,None]).multiply(phase[None,:]).tocsr()
    if (transformed_K-K).nnz:
        raise ArithmeticError('centered Hamiltonian is not particle-hole invariant')
    if not np.array_equal(number[complement],16-number):
        raise ArithmeticError('number did not transform to16-N')
    return {'full_Fock_dimension':len(masks),'source_scale':16,
            'exact_H_transform':'PH(H)=H+8-N',
            'exact_K_transform':'PH(H-N/2)=H-N/2',
            'all_matrix_entries_checked':True,'number_transform':'N ->16-N'}


def sector_bounds():
    # These inequalities use the authenticated existing all-N ground theorem.
    # sqrt3<7/4, sqrt17>4 establish the positive rational floors independently.
    gap_floor=F(1,16)
    improved_floor=F(1,8)
    if not F(103,25)**2<17 or not F(26,15)**2>3:
        raise ArithmeticError('square-root comparison changed')
    if not F(103,50)-F(26,60)-F(3,2)>improved_floor:
        raise ArithmeticError('charged gap floor failed')
    return {'imported_original_H_lower':'-sqrt(3)/4',
            'old_N4_ground_centered_energy':'-2-sqrt(3)/4',
            'N8_Neel_trial_energy':'-4',
            'N8_dimer_product_trial_energy':'-2-sqrt(17)/2',
            'centered_global_ground_number':8,
            'centered_ground_energy_lower':'-4-sqrt(3)/4',
            'centered_ground_energy_upper':'-2-sqrt(17)/2',
            'all_other_N_energy_lower':'-7/2-sqrt(3)/4',
            'charged_sector_gap_lower':'sqrt(17)/2-sqrt(3)/4-3/2',
            'charged_sector_gap_strict_rational_floor':str(improved_floor),
            'Neel_only_gap_lower':'1/2-sqrt(3)/4',
            'Neel_only_gap_strict_rational_floor':str(gap_floor),
            'unique_N8_ground_proved_by_imported_Lieb':True,
            'neutral_gap_from_charge_bounds':False,
            'thermodynamic_ground_uniqueness_proved':False,
            'old_N4_ground_remains_centered_ground':False}


def lieb_premises():
    """Exact hypotheses for the explicitly imported repulsive spin theorem."""
    a,_,even,odd,_,_=ground.band_data()
    if a.dtype.kind not in 'iu' or not np.array_equal(a,a.T) or np.any(np.diag(a)):
        raise ArithmeticError('real symmetric zero-diagonal source hopping required')
    if np.any(a[np.ix_(even,even)]) or np.any(a[np.ix_(odd,odd)]):
        raise ArithmeticError('source graph is not bipartite')
    reached={0}
    while True:
        expanded=reached|{int(y) for x in reached for y in np.flatnonzero(a[x])}
        if expanded==reached:break
        reached=expanded
    if reached!=set(range(8)) or (len(even),len(odd))!=(4,4):
        raise ArithmeticError('connected balanced eight-site graph required')
    if F(source.GEOMETRY.interaction)!=1:
        raise ArithmeticError('strict positive onsite interaction changed')
    return {'theorem':'Lieb, The Hubbard Model: Some Rigorous Results and Open Problems, Theorem1, printed p8 (1993)',
            'primary_review':'https://arxiv.org/abs/cond-mat/9311033',
            'original_result_doi':'https://doi.org/10.1103/PhysRevLett.62.1201',
            'original_erratum_doi':'https://doi.org/10.1103/PhysRevLett.62.1927.5',
            'proof_imported_not_reproved':True,
            'real_symmetric_hopping_checked':True,
            'uniform_hopping_sign_required':False,
            'connected_vertices':sorted(reached),'bipartition_sizes':[4,4],
            'strict_positive_U':'1','particle_number':8,'even_half_filling':True,
            'energy_constant_from_review_convention':'2',
            'ground_hidden_spin':0,'ground_multiplicity':1}


def numerical_diagnostic():
    """Ritz values are diagnostics; sector/ground claims use exact bounds."""
    H=global_ground.source_block(8).astype(float)/16
    H=H-4*sparse.eye(H.shape[0],format='csr')
    v0=np.arange(1,H.shape[0]+1,dtype=float);v0/=np.linalg.norm(v0)
    values,vectors=eigsh(H,k=2,which='SA',tol=1e-11,v0=v0,maxiter=10000)
    order=np.argsort(values);values=values[order];vectors=vectors[:,order]
    residuals=[float(np.linalg.norm(H@vectors[:,j]-values[j]*vectors[:,j])) for j in range(2)]
    if max(residuals)>1e-8:raise ArithmeticError('centered Ritz diagnostic failed residual tolerance')
    return {'balanced_N8_dimension':H.shape[0],
            'lowest_two_Ritz_values':[float(x) for x in values],
            'residual_norms':residuals,'used_for_ground_or_gap_admission':False}


def certificate():
    return {'symmetry':source_symmetry(),'bounds':sector_bounds(),
            'imported_Lieb_theorem':lieb_premises()}
