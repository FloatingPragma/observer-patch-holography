"""Exact source spectral residue, broadening bounds and moment controls.

The block comparisons retain the source R and L_perp. Their finite algebra
controls do not realize the interacting infinite ground or select its atoms.
"""
from fractions import Fraction as F

import numpy as np
import sympy as s

from . import composite, hole_dispersion, nodal_memory


def rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('exact rational input required')
    value=F(value)
    if value<0 or (positive and not value):
        raise ValueError('nonnegative input required; scale and frequency must be positive')
    return value


def score_bounds(alpha,eta,memory_lower,memory_upper,residual_upper,
                 hopping_upper=F(0),composite_upper=F(0)):
    """Original Lorentzian score with all kinetic comparison error retained."""
    alpha,eta=(rational(v,positive=True)for v in(alpha,eta))
    lo,hi,r2,h,c=(rational(v)for v in(memory_lower,memory_upper,residual_upper,hopping_upper,composite_upper))
    if lo>hi or hi>r2:raise ValueError('memory score interval exceeds its total Gram bound')
    beta=alpha*alpha+eta*eta;error=max(h,c)/eta
    lower=eta*eta/beta+alpha*alpha*lo/(beta*(beta+r2))
    upper=eta*eta/beta+alpha*alpha*(beta+r2)*hi/beta**3
    return {'original_score_lower':max(F(0),lower-error),
            'original_score_upper':min(F(1),upper+error),
            'nodal_comparison_score_lower':lower,
            'nodal_comparison_score_upper':min(F(1),upper),
            'kinetic_comparison_error':error,
            'memory_score_lower':lo,'memory_score_upper':hi,
            'complete_original_measure_mass':1,
            'shrinking_window_weight_assumed':False}


def window_from_score(score_lower,tail_ratio):
    """Keep original probability when passing from a Lorentzian score to a window."""
    q=rational(score_lower);n=rational(tail_ratio,positive=True)
    if q>1:raise ValueError('a probability lower bound must not exceed one')
    tail=1/(1+n*n)
    return {'outside_kernel_upper':tail,
            'original_window_probability_lower':max(F(0),(q-tail)/(1-tail)),
            'success_conditioning_used':False}


def default_memory_floor(eta):
    """Actual source moment floor; it vanishes at small unit-hopping frequency."""
    eta=rational(eta,positive=True)
    source=nodal_memory.default_budget(F(1),eta)
    r0=source['neighbourhood_Gram_lower'];m=source['Lperp_R_norm_upper']
    lower=eta*eta*r0*r0/(eta*eta*r0+m*m)
    return {'coupling':F(1),'contact':F(8),'S_upper':source['S_upper'],
            'residual_Gram_lower':r0,'D_R_norm_upper':m,
            'hidden_Lorentzian_score_lower':lower,
            'hidden_score_has_proved_nonzero_limit':False,
            'source_state':'proper_cubic_average_of_invariant_centered_full_ground'}


def subcontact_band():
    """Positive source weight below U/2 on a fixed unit-hopping nodal region."""
    source=nodal_memory.default_budget(F(1))
    assert source['neighbourhood_Gram_lower']==F(18,125)
    assert source['Lperp_R_norm_upper']==5175
    # S<=1, |kappa|<=3, R*R<=18 and ||D R||<=5175.
    m2=F(17);m4=F(833)
    m6=F(81**2+188**2+20780**2)
    third=m6+3*m2*m4+4*m2**3
    variance=16*source['neighbourhood_Gram_lower']
    negative_mean=variance**2/(2*third)
    delta=negative_mean/4
    assert 0<delta<source['S_upper']
    return {'unit_hopping_second_moment_upper':m2,
            'unit_hopping_fourth_moment_upper':m4,
            'unit_hopping_sixth_moment_upper':m6,
            'third_absolute_centered_square_energy_moment_upper':third,
            'square_energy_variance_lower':variance,
            'negative_square_energy_deviation_mean_lower':negative_mean,
            'nodal_region_S_upper':delta,
            'absolute_energy_squared_over_c_squared_upper':16-delta,
            'original_combined_charge_probability_lower':negative_mean/34,
            'comparison_atomic_contact_energy_over_c':4,
            'positive_lower_energy_edge_proved':False,
            'physical_ballistic_window_proved':False,
            'original_charge_probabilities_conditioned':False}


def _simple(a):return s.Matrix(a).applyfunc(s.simplify)


def _positive(a):
    """Exact Hermitian positivity through rational/algebraic principal minors."""
    a=_simple(a)
    if _simple(a-a.H)!=s.zeros(a.rows):raise ArithmeticError('non-Hermitian comparison')
    from itertools import combinations
    for n in range(1,a.rows+1):
        for ids in combinations(range(a.rows),n):
            if s.simplify(a.extract(ids,ids).det()).is_nonnegative is not True:
                raise ArithmeticError('uncertified positive-semidefinite matrix')
    return True


def block(alpha,r,d,h=None,c=None):
    alpha=s.Rational(alpha);r=s.Matrix(r);d=s.Matrix(d)
    n=r.cols;m=r.rows
    if alpha<=0 or d.shape!=(m,m) or _simple(d-d.H)!=s.zeros(m):
        raise ValueError('positive alpha and selfadjoint compatible hidden generator required')
    h=s.zeros(n)if h is None else s.Matrix(h)
    c=s.zeros(n)if c is None else s.Matrix(c)
    if h.shape!=(n,n)or c.shape!=(n,n)or _simple(h-h.H)!=s.zeros(n)or _simple(c-c.H)!=s.zeros(n):
        raise ValueError('selfadjoint response blocks required')
    return s.BlockMatrix([[h,alpha*s.eye(n),s.zeros(n,m)],
                          [alpha*s.eye(n),c,r.H],
                          [s.zeros(m,n),r,d]]).as_explicit()


def _kernel_projection(a):
    cols=a.nullspace()
    if not cols:return s.zeros(a.rows)
    w=s.Matrix.hstack(*cols)
    return _simple(w*(w.H*w).inv()*w.H)


def noncommuting_control():
    """Exact complex comparison with noncommuting memory parts."""
    alpha=s.Integer(4);eta=s.Rational(1,4)
    r=s.Matrix([[1,s.I],[2,1],[0,1-s.I]])/3
    d=s.Matrix([[0,1,s.I],[1,2,1-s.I],[-s.I,1+s.I,-1]])
    L=block(alpha,r,d);n=2;eye=s.eye(n)
    sigma=_simple(r.H*(s.I*eta*s.eye(3)-d).inv()*r)
    M=_simple(r.H*eta**2*(d*d+eta**2*s.eye(3)).inv()*r)
    J=_simple(eta*(sigma+sigma.H)/2)
    if _simple(M*J-J*M)==s.zeros(n):raise ArithmeticError('commuting control cannot test matrix ordering')
    g=(s.I*eta*s.eye(L.rows)-L).inv()[:n,:n]
    q=_simple(-eta*(g-g.H)/(2*s.I));beta=alpha**2+eta**2
    lower=_simple((eta**2*eye+M)*(beta*eye+M).inv())
    _positive(eye-r.H*r);_positive(M-J*J)
    _positive(q-lower)
    linear_lower=eta**2/beta*eye+alpha**2/(beta*(beta+1))*M
    linear_upper=eta**2/beta*eye+alpha**2*(beta+1)/beta**3*M
    _positive(q-linear_lower);_positive(linear_upper-q)
    # The source identity removes neither a real memory part nor its commutator.
    if J==s.zeros(n):raise ArithmeticError('real memory part was discarded')
    return {'response_dimension':n,'hidden_dimension':3,
            'exact_noncommuting_real_imaginary_memory':True,
            'matrix_lower_and_linear_upper_PSD':True,
            'R_norm_squared_admitted_upper':F(1),
            'is_actual_many_body_source_realization':False}


def kernel_control():
    r=s.Matrix([[1,s.I],[1,1],[2,-s.I]])/3
    d=s.diag(0,0,2);alpha=s.Integer(4);L=block(alpha,r,d)
    m0=_simple(r.H*_kernel_projection(d)*r)
    predicted=_simple(m0*(alpha**2*s.eye(2)+m0).inv())
    actual=_kernel_projection(L)[:2,:2]
    if _simple(actual-predicted)!=s.zeros(2):raise ArithmeticError('kernel graph residue formula differs')
    return {'response_dimension':2,'hidden_kernel_dimension':2,
            'full_kernel_dimension':len(L.nullspace()),
            'exact_matrix_residue':[[str(v)for v in row]for row in predicted.tolist()],
            'kernel_graph_identity':True,'physical_mass_assigned':None}


def matched_moment_controls():
    rows=[]
    for weights in ((4,1,1),(4,1,1,s.sqrt(3))):
        size=len(weights)+1;L=s.zeros(size)
        for j,w in enumerate(weights):L[j,j+1]=L[j+1,j]=w
        grading=s.diag(*[(-1)**j for j in range(size)])
        if grading*L*grading!=-L:raise ArithmeticError('energy-reflection grading differs')
        moments=[int(s.simplify((L**n)[0,0]))for n in range(9)]
        residue=_kernel_projection(L)[0,0]
        hidden=L[2:,2:];r=s.zeros(size-2,1);r[0]=1
        m0=(r.T*_kernel_projection(hidden)*r)[0]
        assert s.simplify(residue-m0/(16+m0))==0
        rows.append({'dimension':size,'characteristic_coefficients':[int(x)for x in L.charpoly().all_coeffs()],
                     'original_moments_0_through_8':moments,
                     'residual_Gram':F(1),'hidden_second_moment':F((r.T*hidden**2*r)[0]),
                     'hidden_zero_weight':F(m0),'original_zero_weight':F(residue)})
    assert rows[0]['original_moments_0_through_8'][:7]==rows[1]['original_moments_0_through_8'][:7]
    return {'controls':rows,'actual_Hubbard_source_realizations':False,
            'moment_identity_and_chiral_implication_countercontrol':True}


def source_locality_control():
    """Literal whole-CAR source derivative range, before choosing a state."""
    data=composite.source((3,1,1),F(2),F(3,7));a=data['a'];h=data['H'];checks=0
    for x in range(3):
        for y in range(3):
            for spin in range(2):
                v=a[2*y+spin]
                for n in range(3):
                    anti=composite._anti(a[2*x+spin],v.T)
                    if abs(x-y)>n:
                        if np.any(anti):raise ArithmeticError('source moment extends beyond its derivative range')
                        checks+=1
                    v=composite._comm(v,h)
    v=composite._comm(composite._comm(a[4],h),h)
    if not np.any(composite._anti(a[0],v.T)):
        raise ArithmeticError('range-two source propagation was erased')
    gammas=hole_dispersion.gamma_matrices()
    # gamma_matrices uses lexicographically ordered Cartesian parity cells.
    from itertools import product
    points=list(product(range(2),repeat=3))
    for j in range(3):
        gluing=np.diag([(-1)**x[j]for x in points])
        for k,g in enumerate(gammas):
            expected=(-1 if j==k else 1)*g
            if not np.array_equal(gluing@g@gluing,expected):
                raise ArithmeticError('physical-cell Bloch gluing differs')
    return {'source_sites':3,'whole_Fock_dimension':64,
            'exact_outside_range_zero_checks':checks,'range_two_nonzero_checked':True,
            'three_actual_Bloch_boundary_gluings_checked':True,
            'uniform_low_frequency_modulus_inferred':False}


def certificate():
    return {'source_moment_locality_and_chart':source_locality_control(),
            'noncommuting_resolvent_comparison':noncommuting_control(),
            'nodal_kernel_residue_control':kernel_control(),
            'finite_moment_nonselection':matched_moment_controls(),
            'actual_default_memory_floor_at_one_eighth':default_memory_floor(F(1,8)),
            'actual_default_subcontact_band':subcontact_band(),
            'every_momentum_measure_extension':'finite_range_moments_uniform_integrability_and_determinacy_in_local_Bloch_charts',
            'physical_mass_prediction':None,'physical_species':[],
            'default_shrinking_window_weight_proved':False,
            'uniform_many_body_charge_gap_proved':False}
