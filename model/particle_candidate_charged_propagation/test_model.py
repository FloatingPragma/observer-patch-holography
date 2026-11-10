"""Exact and independent controls for filtered pair records and their boundaries."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_actual_seed_is_neutral_even_and_keeps_its_compensator():
    control=m.seed_control()
    assert control['global_native_increment']==[0]*5
    assert control['global_CAR_parity']==0
    assert control['owner_support_cells']==2
    rows=m.charged.compensated_control()['rows']
    Q=m.charged.compensated_control()['Qy']
    # A regional charge is nonzero even though the complete physical pair is neutral.
    assert all(m.charged.clean(Q*B-B*Q)==-6*B for B in rows)
    assert not control['global_neutral_pair_is_a_charged_odd_field']
    assert not control['boundary_superselection_or_extended_field_algebra_excluded']


def test_original_ideal_window_precedes_the_clock_allowance():
    member=m.contract.member(0)
    r=m.filtered_original_window(member)
    bare=m.charged.regional.charged.bind_member(member)['result']
    assert r['bare_original_ideal_weight_lower']==F(bare['original_transfer_weight_lower'])
    assert r['bare_original_ideal_weight_lower']>F(37,120)
    assert r['ideal_transfer_interval']==(F(39,4),F(135,4))
    assert r['original_residual_squared_upper']<F(129,400)**2
    assert r['original_dressed_ideal_weight_strict_lower']==F(160801,4000000)
    # The previous recorded event has a different interval and probability.
    prior=m.charged.bind_member(member)['result']['original_compensated_dressed_energy_window']
    assert r['ideal_transfer_interval']!=tuple(map(F,prior['reported_transfer_interval']))
    assert r['original_dressed_ideal_weight_strict_lower']>F(prior['original_dressed_recorded_weight_strict_lower'])


def test_triangle_filter_multiplier_and_rational_band_bound():
    # Independent numerical integration of the normalized triangle on a two-state transition.
    T=1/128;omega=135/4;nodes=np.linspace(-T,T,65537)
    density=(1-np.abs(nodes)/T)/T
    integral=np.trapezoid(density*np.exp(1j*omega*nodes),nodes)
    amplitude=np.sinc(T*omega/(2*np.pi))**2
    assert abs(integral-amplitude)<1e-10
    bound=m.triangle_band_retention(F(1,128),F(135,4))
    assert float(bound)<amplitude**2
    assert F(160801,4000000)*bound>F(39,1000)
    with pytest.raises(ValueError):m.triangle_band_retention(1,3)


def test_complete_multichannel_temporal_analyzer_retains_all_failures():
    # Three outcome columns and three controlled times, on an independent exact two-state source.
    H=s.diag(0,1);B0=s.Matrix([[0,s.Rational(1,2)],[0,0]])
    B1=s.Matrix([[0,0],[s.Rational(1,3),0]])
    failure=s.diag(s.sqrt(8)/3,s.sqrt(3)/2)
    V=(B0,B1,failure);weights=(s.Rational(1,2),s.Rational(1,3),s.Rational(1,6))
    v=s.Matrix([s.sqrt(w) for w in weights]);e=s.eye(3)[:,0]
    W=s.eye(3)-(e-v)*(e-v).H/(1-v[0])
    assert s.simplify(W.H*W)==s.eye(3)
    clocks=(s.eye(2),s.diag(1,s.I),s.diag(1,-1))
    rows=[]
    for B in V:
        nodes=[U.H*B*U for U in clocks]
        for j in range(3):
            rows.append(s.simplify(sum((W[j,k]*v[k]*nodes[k] for k in range(3)),s.zeros(2))))
        assert rows[-3]==s.simplify(sum((weights[k]*nodes[k] for k in range(3)),s.zeros(2)))
    assert s.simplify(sum((K.H*K for K in rows),s.zeros(2)))==s.eye(2)
    rho=s.diag(s.Rational(2,3),s.Rational(1,3))
    probs=[s.simplify(s.trace(rho*K.H*K)) for K in rows]
    assert sum(probs)==1 and 0<probs[0]+probs[3]<1
    purification=s.Matrix([s.sqrt(s.Rational(2,3)),0,0,s.sqrt(s.Rational(1,3))])
    joint=s.zeros(4)
    for K in rows:
        vector=s.kronecker_product(K,s.eye(2))*purification
        joint+=vector*vector.H
    reference=s.Matrix(2,2,lambda i,j:s.simplify(sum(joint[2*k+i,2*k+j] for k in range(2))))
    assert reference==rho
    # The distinguished failure does not complete just the distinguished success rows.
    assert s.simplify(sum((rows[k].H*rows[k] for k in (0,3,6)),s.zeros(2)))!=s.eye(2)


def test_source_filter_compiler_controls_column_errors_without_record_factor():
    member=m.contract.member(0);r=m.filter_compiler(member)
    assert r['inner_time_cells']%2==0
    assert r['continuous_success_column_quadrature_error_upper']<=F(1,65536)
    assert r['one_source_evolution_operator_error_upper']<=F(1,65536)
    assert r['distinguished_original_event_probability_error_upper']<=F(6,65536)
    assert r['fine_outcomes']==97*r['inner_time_cells']
    assert not r['energy_zero_changed'] and not r['source_clock_changed']
    with pytest.raises(ValueError):m.filter_compiler(member,allowance=0)


def test_cayley_error_bound_independently_on_two_sided_spectrum():
    # Diagonal exact unitary evolution versus the literal repeated Cayley service.
    eigenvalues=np.array([-11.,-1.,3.,17.]);T=.07;r=3;K=17
    exact=np.exp(-1j*T*eigenvalues)
    cayley=((1-.5j*T*eigenvalues/r)/(1+.5j*T*eigenvalues/r))**r
    assert np.max(np.abs(cayley-exact))<=K**3*T**3/(12*r*r)


def test_recorded_interval_subtracts_each_clock_and_compiler_error_once():
    r=m.filtered_original_window(m.contract.member(0))
    plan=r['complete_source_clock_plan']
    loss=2*F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['two_complete_clock_error_upper'])
    assert loss==r['complete_two_clock_event_loss_upper']<=F(1,8192)
    expected=r['original_time_filtered_ideal_weight_strict_lower']-loss-r['filter_compiler']['distinguished_original_event_probability_error_upper']
    assert expected==r['original_filtered_recorded_weight_strict_lower']>F(3,80)
    assert r['reported_transfer_interval']==[F(311,32),F(1081,32)]
    assert r['full_clock_probe_records']==97*r['filter_compiler']['inner_time_cells']*int(plan['clock_plan']['clock_size'])**2


def test_member_zero_local_positive_event_is_not_inactive_aperture_weight():
    member=m.contract.member(0);plan=m.joint_plan(member)
    assert not plan['active_source_aperture']
    assert plan['success_amplitude_comparison_upper'] is None
    assert plan['original_probability_comparison_upper'] is None
    assert m.filtered_original_window(member)['original_filtered_recorded_weight_strict_lower']>F(3,80)
    plan=m.joint_plan(m.contract.member(8))
    assert plan['active_source_aperture']
    assert plan['seed_support_owners']==2
    assert not plan['limiting_joint_spectral_weight_lower_proved']
    assert not plan['state_subsequence_convergence_modulus_supplied']
    assert not plan['finite_Fourier_bin_is_a_conserved_momentum_eigenvalue']
    with pytest.raises(ValueError):m.filtered_original_window(m.contract.member(1))


def test_fixed_filter_can_erase_finite_capacity_positive_high_energy_weight():
    # Unit original weight on E=1/a is insufficient for a surviving fixed-width response.
    T=F(1,128);lower=F(39,4)
    values=[m.ultraviolet_filter_upper(F(1,2**n),T,lower) for n in (8,10,12)]
    assert values[0]>values[1]>values[2]>0
    assert values[0]/values[1]==256
    for n,bound in zip((8,10,12),values):
        omega=float(lower*2**n)
        actual=np.sinc(float(T)*omega/(2*np.pi))**4
        assert actual<=float(bound)
    with pytest.raises(ValueError):m.ultraviolet_filter_upper(0)


def test_same_member_identity_and_resealed_physical_promotion_rejected():
    member=m.contract.member(0);value=m.bind_member(member)
    assert all(value['result'][key]==expected for key,expected in member['ids'].items())
    assert value['result']['physical_particle_map']==[]
    assert value['result']['physical_mass_table_GeV']==[]
    changed=deepcopy(value);changed['result']['physical_mass_table_GeV']=[{'mass':1}]
    changed['result_sha256']=m.contract.digest(changed['result'])
    with pytest.raises(ValueError):m.verify_bound(member,changed)
    with pytest.raises(ValueError):m.verify_bound(m.contract.member(1),value)


def test_joint_plan_uses_forward_bond_collar_and_complete_cartesian_bins():
    member=m.contract.member(8);r=m.joint_plan(member,spatial_index=(-2,1,0),time_index=-1)
    p=r['parameters'];a=F(member['source']['spacing'])
    assert r['same_source_bulk_bound']['support_radius']==p['spatial_width']+p['spatial_aperture']+a
    assert r['fine_outcomes_per_origin']==2*r['inner_analyzer_dimension']*r['outer_analyzer_dimension']
    assert r['distinguished_Fourier_bin']==[510,1,0,255]
    assert r['one_common_origin_retained_across_history'] and r['quantum_reference_retained']
    with pytest.raises(ValueError):m.joint_plan(member,spatial_index=(0,0,False))
    with pytest.raises(ValueError):m.joint_plan(member,spatial_width='1/513')


def test_resolved_apparatus_error_is_not_a_state_or_weight_estimate():
    plan=m.joint_plan(m.contract.member(10))
    assert plan['active_source_aperture']
    assert plan['original_probability_comparison_upper']<F(1,50)
    assert not plan['state_subsequence_convergence_modulus_supplied']
    assert not plan['limiting_joint_spectral_weight_lower_proved']
