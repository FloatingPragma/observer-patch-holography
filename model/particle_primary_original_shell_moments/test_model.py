"""Source-bound Ward shells, exact spectral controls and complete clock errors."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import pytest
import sympy as s

from . import model as m


def test_exact_source_shell_geometry_and_all_native_blocks():
    x=m.exact_geometry()
    assert x['source_cube_pair_orbits']['ordered_pair_orbit_sizes']==[8,24,24,8]
    assert x['complete_native_dimension_sum']==15
    assert x['density_first_moment_shell_coefficients_of_negative_bond_expectation']==['3/2','-1/2','0','0']
    assert x['Walsh_eigenmoment_coefficients']==[0,1,2,3]


def test_complete_sixteen_record_menu_for_symbolic_full_occupations():
    cells=m.kernels.CELLS;n=s.symbols('n0:8',real=True)
    successes=[sum((-1)**sum(a*b for a,b in zip(q,x))*n[j] for j,x in enumerate(cells))/240 for q in cells]
    failure_effects=[(1-v*v/900)/8 for v in n]
    assert s.expand(sum(k*k for k in successes)+sum(failure_effects))==1
    # A missing address factor overcounts the complete success effect.
    wrong=s.expand(8*sum(k*k for k in successes)+sum(failure_effects))
    assert wrong!=1
    for bits in product((-30,30),repeat=8):
        values=dict(zip(n,bits))
        assert sum(k.subs(values)**2 for k in successes)==1


def test_independent_energy_eigenbasis_original_kernel():
    x=m.exact_spectral_countercheck()
    assert x['negative_bond_expectation']==F(2,19)
    assert x['four_first_moments']==[F(),F(2,19),F(4,19),F(6,19)]
    assert x['all64_original_spectral_kernel_entries_checked']
    assert not x['control_preparation_replaces_primary_polynomial_density']


def test_density_Ward_identity_with_full_CAR_and_nonconstant_link():
    # Exact second quantization of an arbitrary nonreal two-site block.
    create=[m.kernels.selected.odd.creation_word((i,),4) for i in range(4)]
    N1=sum((a*a.H for a in create[:2]),s.zeros(16));N2=sum((a*a.H for a in create[2:]),s.zeros(16))
    block=s.Matrix([[1,s.I],[s.Rational(1,3),-s.I]])
    hop=sum((block[i,j]*create[i]*create[2+j].H for i,j in product(range(2),repeat=2)),s.zeros(16))
    link=s.Matrix([[0,0,0],[1,0,0],[0,1,0]])
    T=s.kronecker_product(link,hop)+s.kronecker_product(link.H,hop.H)
    H=T+s.kronecker_product(s.diag(0,7,19),s.eye(16))+s.kronecker_product(s.eye(3),N1-2*N2)
    X=s.kronecker_product(s.eye(3),N1);Y=s.kronecker_product(s.eye(3),N2)
    comm=lambda A,B:A*B-B*A
    assert comm(X,comm(H,Y))==T
    assert comm(X,comm(H,X))==-T
    assert comm(X+Y,H)==s.zeros(48)


def test_identical_complete_clock_exchange_keeps_constant_row_mean_zero():
    H=s.diag(-1,0,2);I=s.eye(3);step=s.Rational(1,7)
    U=(I-s.I*step*H/2).inv()*(I+s.I*step*H/2)
    clock=m.spectroscopy.energy_instrument(U,4)
    rho=s.diag(s.Rational(2,3),s.Rational(1,6),s.Rational(1,6))
    K=s.diag(s.Rational(1,3),s.Rational(1,2),1);labels=(0,1,-2,-1)
    probabilities={}
    for a,c in product(range(4),repeat=2):
        J=clock[c]*K*clock[a]
        probabilities[a,c]=s.simplify(s.trace(J*rho*J.H))
    assert all(probabilities[a,c]==probabilities[c,a] for a,c in probabilities)
    assert s.simplify(sum((labels[c]-labels[a])*p for (a,c),p in probabilities.items()))==0
    assert s.simplify(sum(probabilities.values())-s.trace(rho*K*K))==0


def test_signed_mean_error_uses_range_and_keeps_original_weights():
    assert m.finite_record_error(3,F(1,10),F(1,100),F(1,200))==F(22,5)
    # Moving TV mass delta between opposite endpoints changes its signed
    # expectation by twice the absolute label bound, not by that bound.
    R=F(3);delta=F(1,200)
    assert delta*(28*R-(-28*R))==56*R*delta
    with pytest.raises(ValueError):m.finite_record_error(0,1,0,0)


def test_actual_original_member_and_complete_clock_bindings():
    member=m.contract.member(0);receipt=m.bind_member(member);m.verify_bound(member,receipt)
    x=receipt['result'];rows=x['spatial_shell_records'];clock=x['complete_two_energy_clock_plan']
    gamma=F(x['negative_bond_expectation_interval'][0]);unit=gamma/7200
    assert gamma>F(2373,10000)
    assert F(x['finite_signed_mean_error_upper'])<unit/2
    assert F(x['neighboring_nonzero_shell_reported_mean_difference_lower'])>0
    assert x['records_before_source_clocks']==16 and x['probe_CAR_parity']==0
    assert rows[0]['reported_original_signed_transfer_first_moment_interval']==['0','0']
    assert all(F(r['reported_original_signed_transfer_first_moment_interval'][0])>0 for r in rows[1:])
    assert F(clock['energy_resolution'])==unit/16
    assert x['source_input_is_original_selected_density'] and x['all_failures_and_reference_retained']
    assert x['zero_and_nonzero_shell_normalized_ideal_energy_shapes_differ']
    assert not x['nonzero_shell_normalized_energy_shapes_proved_distinct']
    assert not x['odd_G1_H_trion_shell_values_computed'] and x['physical_mass_table_GeV']==[]
    assert x['zero_shell_source_success_weight_lower_factored']==m.regional.common_contract_payload(member,1)['physical_density_floor_factored']
    native=x['native_hopping_sign_and_density_shells']
    assert native['every_native_total_hopping_expectation_strictly_negative']
    assert [r['native_dimension'] for r in native['native_density_channels']]==[3,2,6,3,1]
    for row in native['native_density_channels']:
        d=row['native_dimension'];bound=row['negative_native_bond_expectation_lower_factored']
        assert F(bound['rational_multiplier'])==2*d*F(native['original_selected_polynomial_slope']['modular_slope_lower'])
        assert bound['negative_exponent']==x['zero_shell_source_success_weight_lower_factored']['negative_exponent']
        assert row['four_original_first_moment_factors_of_tau_r']==[str(F(j,32*d*d)) for j in range(4)]
    event=x['complete_density_success_and_conditional_statistics']
    assert F(event['original_complete_density_success_weight_upper'])<F(21,1000)
    assert F(event['conditional_on_density_success_ideal_signed_mean_lower'])>F(19,1000)
    assert F(event['conditional_on_density_success_reported_signed_mean_lower'])>F(15,1000)
    assert not event['conditional_statistic_is_original_event_weight']


def test_native_physical_anchor_norms_and_actual_bipartite_parity():
    c=m.native_sign_control()
    assert c['physical_filled_beta_anchor']['Gauss_generators_checked']==12
    assert c['physical_filled_beta_anchor']['defect_count']==0
    assert c['actual_bipartite_bonds_checked']==12
    for row in c['all_native_first_action_packets']:
        assert row['two_added_defect_packet_norm_squared']==12*row['dimension']
        assert row['unchanged_defect_packet_norm_squared']==0
        assert all(axis[2]==row['dimension'] and axis[0]==0 for axis in row['affected_axis_norms'])


def test_native_sign_from_polynomial_slope_with_all_fock_sectors():
    create=[m.kernels.selected.odd.creation_word((i,),4) for i in range(4)]
    N=[a*a.H for a in create];I=s.eye(16)
    hops=[create[0]*create[2].H+create[2]*create[0].H,
          2*(create[1]*create[3].H+create[3]*create[1].H)]
    H=hops[0]+hops[1]+N[0]-N[2]
    # ||H|| <= 1+2+1 =4; the source analogue uses its bound and selected
    # Taylor polynomial instead of this separate exact control.
    Q=8*I-H;Z=s.trace(Q**4);rho=Q**4/Z;slope=4*4**3/Z
    for r,T in enumerate(hops):
        U=I-2*N[r]
        assert U*H*U==H-2*T
        comm=H*U-U*H;HS=s.trace(comm.H*comm)
        work=s.trace(rho*(U*H*U-H))
        assert work==-2*s.trace(rho*T)>0
        assert work>=slope*HS/2
        assert -s.trace(rho*T)>=slope*s.trace(T*T)


def test_density_success_full_occupation_inequality():
    # Enumerate every pair of positive-beta particles and negative-beta
    # holes. Degeneracy does not change the two commuting count values.
    for p,h in product(range(31),repeat=2):
        assert (p-h)**2<=30*(p+h)
    assert (30-0)**2==30*(30+0)
    member=m.contract.member(0);x=m.bind_member(member)['result']
    event=x['complete_density_success_and_conditional_statistics']
    cap=F(event['source_defect_count_upper'])/240
    assert F(event['original_complete_density_success_weight_upper'])==cap
    ideal=F(x['negative_bond_expectation_interval'][0])/600
    reported=ideal-7*F(x['finite_signed_mean_error_upper'])
    assert F(event['conditional_on_density_success_ideal_signed_mean_lower'])==ideal/cap
    assert F(event['conditional_on_density_success_reported_signed_mean_lower'])==reported/cap


def test_authenticated_scope_and_source_mutations_rejected():
    member=m.contract.member(0);receipt=m.bind_member(member)
    changed=deepcopy(receipt);changed['result']['odd_G1_H_trion_shell_values_computed']=True
    changed['result_sha256']=m.contract.digest(changed['result'])
    with pytest.raises(ValueError):m.verify_bound(member,changed)
    changed=deepcopy(member);changed['source']['electric_coefficient']='106'
    with pytest.raises(ValueError):m.bind_member(changed)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
