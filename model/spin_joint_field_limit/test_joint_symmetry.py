from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.spin_joint_field_limit import joint_symmetry as j


def test_exact_source_algebra_and_fourth_derivative_certificate():
    r=j.algebra_checks()
    assert r['actual_native_two_site_source_attached']
    assert r['all_three_source_double_commutators_checked']
    assert r['cubic_group_commutator_checked']
    assert F(r['fourth_derivative_coefficient_sum'])<4


@pytest.mark.parametrize('a,b',[(.2,.3),(-.4,.1),(.05,-.02),(-.3,-.4)])
def test_cubic_group_commutator_bound_against_independent_exponentials(a,b):
    x=np.array([[0,1],[1,0]],complex)/2
    y=np.array([[0,-1j],[1j,0]])/2
    z=np.diag([1,-1])/2
    actual=expm(1j*a*y)@expm(1j*b*z)@expm(-1j*a*y)@expm(-1j*b*z)
    expansion=expm(-1j*a*b*x)@expm(-.5j*a*a*b*z)@expm(.5j*a*b*b*y)
    assert np.linalg.norm(actual-expansion,2)<=2*(abs(a)+abs(b))**4+1e-14


def test_second_order_C_P_remainder_has_its_source_sign():
    x=np.array([[0,1],[1,0]],complex)/2
    y=np.array([[0,-1j],[1j,0]])/2
    z=np.diag([1,-1])/2
    for n in (2,4,8):
        p=1/(n*n);a=p*p;b=-p
        actual=expm(1j*a*x)@expm(1j*b*y)@expm(-1j*a*x)@expm(-1j*b*y)
        expected=expm(-1j*a*b*z)
        assert np.linalg.norm(actual-expected,2)<=(a*a*abs(b)+abs(a)*b*b)/4+1e-14


def test_joint_family_exact_powers_and_small_kicks():
    for n in (2,5,11,100):
        a=j.auxiliary(n,27,F(2),F(3));c=j.longitudinal(n,27,F(2),F(3),F(5))
        assert a['p']==F(1,n*n) and a['field_over_J']==F(1,n**5)
        assert a['vacuum_norm_squared_upper']==12*27*9/F(n**5)+2*27*4/F(n**6)
        assert a['scaled_kick_kappa_squared']==F(3,2)*27*9/F(n**8)+27*4/F(4*n**9)
        assert c['variance_upper']==12*27*9/F(n)
        assert c['absolute_mean_upper']==3*27*5/F(8*n*n)
        assert not c['unbounded_C_strong_history_convergence_claimed']


def test_all_stage_errors_decay_without_small_C_operator_norm():
    values=[]
    for n in (2,4,8,16,32,64,128):
        a=j.auxiliary(n,27,F(1),F(1));c=j.longitudinal(n,27,F(1),F(1),F(1))
        values.append((a['vacuum_norm_squared_upper'],c['vacuum_norm_squared_upper'],
            j.time_phase_squared(a['scaled_kick_kappa_squared'],1,1,2),
            j.time_phase_squared(c['scaled_kick_kappa_squared'],1,1,2),
            j.C_P_rotation_remainder(n,27,1,1),j.P_P_rotation_remainder(n,27,1,1)))
        assert c['operator_norm_upper']==F(27*n*n,2)
    assert all(all(a>b for a,b in zip(x,y))for x,y in zip(values,values[1:]))


def test_rational_square_root_bound_is_certified():
    for x in (F(),F(1,3),F(17,9),F(10**12),F(1,10**18)):
        upper=j.sqrt_upper(x);assert upper*upper>=x
        if upper:assert (upper-F(1,1<<40))**2<x


@pytest.mark.parametrize('call',[
    lambda:j.auxiliary(0,27,1,1),lambda:j.auxiliary(True,27,1,1),
    lambda:j.longitudinal(2,0,1,1,1),lambda:j.longitudinal(2,27,-1,1,1),
    lambda:j.time_phase_squared(1,-1,1,1),lambda:j.sqrt_upper(1,-1),
    lambda:j.P_P_rotation_remainder(1,27,1,1),
])
def test_unproved_domain_rejected(call):
    with pytest.raises(ValueError):call()


def test_internal_U1_preserves_two_pair_symplectic_form():
    for theta in (.1,.7,1.9):
        rotation=np.array([[math.cos(theta),math.sin(theta)],[-math.sin(theta),math.cos(theta)]])
        block=np.block([[rotation,np.zeros((2,2))],[np.zeros((2,2)),rotation]])
        symplectic=np.block([[np.zeros((2,2)),np.eye(2)],[-np.eye(2),np.zeros((2,2))]])
        assert np.linalg.norm(block.T@symplectic@block-symplectic)<1e-14


def test_certificate_retains_joint_family_and_proof_order():
    r=j.certificate()
    assert r['proof_order']==['auxiliary_transverse_phases','longitudinal_phases','mixed_momentum_Weyl_relation','regularity_and_internal_U1']
    assert r['parent_vacuum_moment_bounds_used_only_for_regularity']
    assert not r['P3_regularity_assumed_before_CCR']
    assert r['common_spatial_action_uses_even_source_translations']
    assert not r['joint_history_state_identified_with_original_ordered_limit']
    assert not r['physical_Lorentz_spin_or_mass_assigned']


def test_rotated_source_momentum_splitting_against_independent_exponentials():
    sy=np.array([[0,-1j],[1j,0]])/2
    sz=np.diag([1,-1])/2
    for a,b in ((.2,.3),(-.4,.1),(.05,-.02)):
        x=1j*a*sy;y=1j*b*sz;c=x@y-y@x
        cubic=(x@c-c@x)/6+(y@c-c@y)/3
        actual=expm(x+y)
        split=expm(x)@expm(y)@expm(-c/2)@expm(cubic)
        assert np.linalg.norm(actual-split,2)<=100*(abs(a)+abs(b))**4+1e-14
