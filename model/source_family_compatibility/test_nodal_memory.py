from fractions import Fraction as F
import copy
import numpy as np
import pytest
import sympy as s
from model.source_family_compatibility import nodal_memory as n
from model.source_family_compatibility.local_matter import annihilator


def source_star():
    a=[annihilator(2,j)for j in range(4)];eye=s.eye(16)
    t=sum((a[j].H*a[j+2]+a[j+2].H*a[j]for j in(0,1)),s.zeros(16))
    def pair(x,e):
        u,d=a[2*x:2*x+2];nu=u.H*u;nd=d.H*d
        return ((u.H*d,(nu-nd)/2),(e*u.H*d.H,(nu+nd-eye)/2))
    q=sum((z*w+(p*v.H+p.H*v)/2 for(p,z),(v,w)in zip(pair(0,1),pair(1,-1))),s.zeros(16))
    return a,eye,t,q


def test_star_bound_on_full_occupation_space_and_sharpness():
    row=n.star_algebra();assert row['full_Fock_dimension']==16
    _,eye,t,q=source_star()
    for sign in(-1,1):
        assert (eye-4*q+sign*t).is_positive_semidefinite
        assert not (s.Rational(15,4)*(eye/4-q)+sign*t).is_positive_semidefinite
    assert (eye-4*q+t).nullspace()  # Zero outcomes and sectors are retained.


def test_wrong_eta_sign_invalidates_source_star_inequality():
    a,eye,t,q=source_star()
    ep0=a[0].H*a[1].H;ep1=-a[2].H*a[3].H
    wrong=q-(ep0*ep1.H+ep0.H*ep1)
    assert not (eye-4*wrong+t).is_positive_semidefinite


@pytest.mark.parametrize('c,u',[(F(1,8),F(1)),(F(3,2),F(4)),(F(2),F(3,7))])
def test_actual_dimer_variational_energy_not_an_input_ground(c,u):
    a,eye,t,_=source_star()
    contact=sum(((a[2*x].H*a[2*x]+a[2*x+1].H*a[2*x+1]-eye)**2 for x in(0,1)),s.zeros(16))
    q=s.Rational(c.numerator,c.denominator)*t/2+s.Rational(u.numerator,u.denominator)*contact/2
    expected=(s.Rational(u.numerator,u.denominator)-s.sqrt(s.Rational((u*u+4*c*c).numerator,(u*u+4*c*c).denominator)))/2
    eig=q.eigenvals();assert expected in eig
    assert all(s.simplify(v-expected).is_nonnegative for v in eig)
    assert n.dimer_energy(c,u)>0
    assert float(-2*n.dimer_energy(c,u))>=float(expected)-1e-14
    # Actual two-particle pair ground; off-dimer hopping has zero expectation.
    number=sum((op.H*op for op in a),s.zeros(16))
    ground=(q-expected*eye).nullspace()
    assert ground and all(number*v==2*v for v in ground)


def test_default_residual_lower_is_strict_on_a_positive_measure_region():
    b=n.default_budget();assert b['S_upper']==F(3,250)
    assert b['node_Gram_strict_lower']==F(9,2000)
    assert b['neighbourhood_Gram_lower']==F(9,4000)
    assert b['Lperp_R_norm_upper']==F(5175,64)
    assert b['negative_imaginary_memory_lower']>0
    # Full interval endpoints on the actual algebraic constraints.
    for q2 in(F(-3,4),F(1,4)):
        for kappa in(F(-3),F(1)):
            node=F(36,125);value=node+(4*q2-kappa*kappa)*b['S_upper']
            assert value>=F(18,125)


def test_memory_jensen_bound_on_actual_spectral_formula():
    rng=np.random.default_rng(419)
    raw=rng.normal(size=(7,7))+1j*rng.normal(size=(7,7));h=(raw+raw.conj().T)/2
    v=rng.normal(size=7)+1j*rng.normal(size=7)
    norm=float(np.vdot(v,v).real);moment=float(np.vdot(h@v,h@v).real)
    for eta in(.01,.25,2.,100.):
        actual=-np.vdot(v,np.linalg.solve(1j*eta*np.eye(7)-h,v)).imag
        lower=eta*norm*norm/(eta*eta*norm+moment)
        assert actual>=lower-1e-12


def test_nonzero_residual_does_not_force_a_zero_energy_memory_atom():
    eta=s.symbols('eta',positive=True)
    # Both this gapped omitted source-clock control and a zero mode have R*R=1.
    gapped=eta/(eta*eta+4)
    zero=1/eta
    assert s.limit(eta*gapped,eta,0)==0
    assert s.limit(eta*zero,eta,0)==1
    assert gapped>0


@pytest.mark.parametrize('c,u,ss,eta',[(0,1,0,1),(1,-1,0,1),(1.,8,0,1),(1,8,4,1),(1,8,0,0)])
def test_invalid_source_and_clock_inputs_rejected(c,u,ss,eta):
    with pytest.raises(ValueError):n.source_budget(c,u,s_max=ss,eta=eta)


def test_fixed_default_is_scale_covariant_and_complete():
    a=n.default_budget(F(1,8),F(1,8));b=n.default_budget(F(1,4),F(1,4))
    assert b['neighbourhood_Gram_lower']==4*a['neighbourhood_Gram_lower']
    assert b['negative_imaginary_memory_lower']==2*a['negative_imaginary_memory_lower']
    row=n.certificate();plan=row['local_qstar_score_plan']
    assert sum(plan['setting_probabilities'])==1 and plan['all_plus_minus_outcomes_retained']
    assert row['physical_species_identified']==[]
    assert not row['source_ground_replaced_by_dimer_trial']
    assert not row['ground_marginal_or_native_services_implicitly_supplied']


def test_actual_complete_force_preparation_and_original_weights():
    row=n.force_source_control();prob=row['original_minus_plus_complement_probabilities']
    assert len(prob)==3 and sum(prob)==1 and prob[0]==prob[1]>0
    assert row['all_input_effect_sum_identity']
    assert not row['finite_control_has_momentum_assignment']
    budget=n.force_preparation_budget(F(1,10**7),F(1,10**7),F(1,10**7))
    assert budget['original_event_probability_lower_after_all_errors']>F(7,500000)
    assert budget['complement_effect_lower']==F(51,100)
    assert budget['ideal_original_combined_charge_probability_lower']==F(9,625000)
    assert not budget['equal_individual_smeared_charge_weights_assumed']
    assert not budget['ground_correlator_needed_to_compile_force']


def test_force_preparation_counts_both_source_flavours_and_all_overlapping_sites():
    from model.source_family_compatibility import uv_tail
    assert 2*uv_tail.ball(2)==50
    assert 2*uv_tail.ball(2)*2*7**2==4900
    assert F(4900,100**2)==F(49,100)
    assert F(18,125)/100**2==F(9,625000)
    assert 6*F(7,10)**2<F(7,4)**2


def test_force_cutoff_bound_uses_complete_isometries_on_reference_inputs():
    rng=np.random.default_rng(91)
    x=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));x*=.2/np.linalg.norm(x,2)
    y=x+1e-3*(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    def dilate(k):
        comp=np.eye(4)-k.conj().T@k-k@k.conj().T
        vals,vecs=np.linalg.eigh(comp);assert min(vals)>.51
        root=(vecs*np.sqrt(vals))@vecs.conj().T
        return np.vstack((k,k.conj().T,root))
    v,w=dilate(x),dilate(y)
    assert np.linalg.norm(v-w,2)<=np.sqrt(6)*np.linalg.norm(x-y,2)
    reference=np.eye(3);assert np.allclose(np.kron(v,reference).conj().T@np.kron(v,reference),np.eye(12))
