"""Exact local CAR and frame identities behind the infinite-state theorem."""
from fractions import Fraction as F
from itertools import product

import numpy as np

from ..field_symmetry import symmetry as s


def anti(a, b):
    return a@b+b@a


def contact_identities():
    a=[s._annihilator(2,j) for j in range(2)]
    n=[x.adjoint()@x for x in a]
    b=[n[1]@a[0],n[0]@a[1]]
    expected=((n[1],s._scale(a[1].adjoint()@a[0],-1)),
              (s._scale(a[0].adjoint()@a[1],-1),n[0]))
    checked=0
    for j,k in product(range(2),repeat=2):
        value=anti(b[j],b[k].adjoint())
        if (not (value-anti(a[j],b[k].adjoint())).is_zero() or
                not (value-anti(b[j],a[k].adjoint())).is_zero() or
                not (value-expected[j][k]).is_zero()):
            raise ArithmeticError('full onsite contact CAR identity failed')
        checked+=1
    # Disjoint odd contact polynomials anticommute, including mixed spins.
    aa=[s._annihilator(4,j) for j in range(4)]
    nn=[x.adjoint()@x for x in aa]
    bb=[nn[j^1]@aa[j] for j in range(4)]
    for j,k in product(range(2),range(2,4)):
        if (not anti(bb[j],bb[k].adjoint()).is_zero() or
                not anti(aa[j],bb[k].adjoint()).is_zero()):
            raise ArithmeticError('different-site contact CAR locality failed')
    return {'onsite_Fock_dimension':4,'all_spin_pairs_checked':checked,
            'two_site_Fock_dimension':16,'different_site_mixed_pairs_checked':4,
            'operator_identity':'{b_sigma,b_tau^dagger}={c_sigma,b_tau^dagger}={b_sigma,c_tau^dagger}',
            'contact_operator':'b_sigma=n_opposite c_sigma',
            'one_body_convention':'rho_sigma,tau=omega(c_tau^dagger c_sigma)',
            'matrix_R':'R=tr(rho) I-rho; block diagonal by site',
            'dilated_identity':'V^dagger J=J^dagger V=J^dagger J=R',
            'exact_memory_Gram':'W^dagger W=U^2(R-R^2)',
            'arbitrary_state_bound':'0<=W^dagger W<=U^2 I/4',
            'diagonal_density_only_is_sufficient':False}


def frame_translation_certificate():
    # The source frame is exactly period four in each coordinate. These
    # 4096 cases therefore certify all integer x,t modulo central signs.
    sites=list(product(range(4),repeat=3))
    frames={x:np.asarray(s._frame(x),dtype=np.complex128) for x in sites}
    if any(not np.all(np.isin(v.real,(-1,0,1))) or
           not np.all(np.isin(v.imag,(-1,0,1))) for v in frames.values()):
        raise ArithmeticError('source frame lost Gaussian-integer entries')
    counts={'+':0,'-':0}
    for x,t in product(sites,repeat=2):
        shifted=frames[tuple((a+b)%4 for a,b in zip(x,t,strict=True))]
        candidate=frames[t]@frames[x]
        if np.array_equal(shifted,candidate):counts['+']+=1
        elif np.array_equal(shifted,-candidate):counts['-']+=1
        else:raise ArithmeticError('translation failed to normalize hidden SU2')
    return {'frame':'R_x=(i sigma_x)^x (i sigma_y)^y (i sigma_z)^z',
            'exact_translation_cases':sum(counts.values()),'central_sign_counts':counts,
            'identity':'R_(x+t)=sign(x,t) R_t R_x, sign in{+1,-1}',
            'consequence':'Spatial translations normalize the source global hidden SU2; central signs cancel in onsite conjugation.',
            'ground_join':'Every normalized full-box ground projector and exterior vacuum is hidden-SU2 invariant. Translates, averages and weak-star limits retain this invariance.',
            'limiting_onsite_covariance':'rho_x=r I with r independent of x and 1/8<=r<=1/2',
            'memory_Gram_lower_factor':'7/64','memory_Gram_upper_factor':'1/4',
            'ground_state_functional_faithfulness_inferred':False,
            'hidden_SU2_identified_as_physical_spin':False}


def memory_factor(density_per_spin):
    if type(density_per_spin) not in (int,F) or not F(1,8)<=density_per_spin<=F(1,2):
        raise ValueError('exact selected-ground per-spin density in[1/8,1/2] required')
    r=F(density_per_spin)
    return r*(1-r)


def certificate():
    return {'local_CAR':contact_identities(),
            'selected_ground_symmetry':frame_translation_certificate(),
            'density_examples':[{'r':str(r),'memory_factor':str(memory_factor(r))}
                                for r in (F(1,8),F(1,4),F(1,2))],
            'strictly_nonzero_infinite_ground_response_memory':True,
            'same_CAR_number_conserving_static_quadratic_response_excluded':True,
            'nontrivial_scattering_proved':False,
            'operator_norm_large_frequency_limit_claimed':False,
            'large_frequency_statement':'i y Sigma(i y) converges strongly to W^dagger W as y->infinity; nonzero Gram excludes identically zero memory.'}
