"""Independent finite profiles, selected expectations and scope controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest

from . import model as m


def _grid(level,widths):
    R=2**level;L=2*R*R
    coordinates=np.indices((L,L,L))
    f=np.ones((L,L,L))
    for i,b in enumerate(widths):
        distance=np.minimum(coordinates[i],L-coordinates[i])
        f*=np.maximum(1-distance/(b*R),0)
    return f


def _exact_work_coefficient(level,widths=(1,1,1),amplitude=1):
    a=2.**(-level);f=_grid(level,widths)
    # Evaluate 1-cos(x) as 2 sin²(x/2), avoiding cancellation at small a.
    return sum(np.sum(2*np.sin(np.pi*float(amplitude)*a*(np.roll(f,-1,axis=i)-f)/2)**2)
        for i in range(3))/(3*a)


def test_exact_tent_moments_include_half_period_endpoint():
    for radius in range(1,13):
        for period in (2*radius,2*radius+5):
            t=[max(F(),1-F(min(k,period-k),radius)) for k in range(period)]
            delta=[t[(k+1)%period]-t[k] for k in range(period)]
            actual=m.tent_moments(radius)
            assert actual=={'square_sum':sum(x*x for x in t),
                'fourth_sum':sum(x**4 for x in t),
                'gradient_square_sum':sum(x*x for x in delta),
                'gradient_fourth_sum':sum(x**4 for x in delta)}


def test_profile_geometry_matches_independent_grid_and_explicit_error():
    for n,b,s in ((0,(1,1,1),F(1)),(0,(1,1,1),F(2)),
                  (1,(2,1,1),F(1,2)),(1,(2,2,2),F(1)),
                  (2,(1,2,3),F(3,2)),(3,(2,1,1),F(1))):
        row=m.profile_geometry(n,b,s);f=_grid(n,b)
        for axis in range(3):
            delta=np.roll(f,-1,axis=axis)-f
            assert np.isclose(np.sum(delta**2),float(row['axis_gradient_square_sums'][axis]),rtol=2e-14)
            assert np.isclose(np.sum(delta**4),float(row['axis_gradient_fourth_sums'][axis]),rtol=2e-14)
        exact=_exact_work_coefficient(n,b,s)
        quadratic=np.pi**2*float(row['quadratic_work_over_selected_density_pi2_coefficient'])
        target=np.pi**2*float(row['limiting_work_over_selected_density_pi2_coefficient'])
        assert -1e-12<=quadratic-exact<=float(row['cosine_nonlinearity_error_upper'])+1e-12
        assert abs(exact-target)<=float(row['absolute_error_from_pi2_limit_upper'])+1e-12
        assert float(row['work_over_selected_density_lower'])-1e-12<=exact<=float(row['work_over_selected_density_upper'])+1e-12
    # A genuine coarse-member alias must not be labeled a positive response.
    assert _exact_work_coefficient(0,amplitude=2)<1e-28
    assert not m.profile_geometry(0,amplitude=2)['small_gradient_lower_bound_applies']


def test_exact_spatial_error_orders_and_work_ratios():
    last=None
    for n in (2,3,4,8):
        row=m.profile_geometry(n,(2,1,1));a=F(1,2**n)
        assert row['discretization_pi2_coefficient']/a**2<F(2)
        assert row['cosine_nonlinearity_error_upper']/a**4<F(9,4)
        ratio=m.profile_ratio(n,(2,1,1))
        assert ratio['limiting_work_ratio']==F(3,2)
        if n<=3:
            exact=_exact_work_coefficient(n,(2,1,1))/_exact_work_coefficient(n)
            assert float(ratio['actual_original_work_ratio_lower'])<=exact<=float(ratio['actual_original_work_ratio_upper'])
        width=ratio['actual_original_work_ratio_upper']-ratio['actual_original_work_ratio_lower']
        if last is not None:assert width<last
        last=width
    assert m.profile_ratio(2,(2,2,2))['limiting_work_ratio']==2
    assert m.profile_ratio(2,amplitude=2)['limiting_work_ratio']==4
    assert m.profile_geometry(2,(1,2,3))['three_axis_gradient_square_sum']==m.profile_geometry(2,(3,1,2))['three_axis_gradient_square_sum']


def test_selected_complete_work_identity_and_retained_reference():
    vertices=list(product(range(2),repeat=3));H=np.zeros((8,8))
    for x,u in enumerate(vertices):
        for y,v in enumerate(vertices):
            if sum(a!=b for a,b in zip(u,v))==1:H[x,y]=-1
    Q=np.linalg.matrix_power(4*np.eye(8)-H,3)
    rho=Q@Q;rho/=np.trace(rho)
    rho_root=Q/np.sqrt(np.trace(Q@Q))
    assert np.max(abs(H@rho-rho@H))<1e-14
    h=-np.trace(rho@H)/8
    for s in (F(1,2),F(1),F(2)):
        work=0.;reference=np.zeros((8,8),complex);complete=np.zeros((8,8),complex)
        for X in range(8):
            for sign in (-1,1):
                phases=np.ones(8,dtype=complex);phases[X]=np.exp(1j*sign*np.pi*float(s))
                U=np.diag(phases);K=U/4
                work+=np.trace(rho@(U.conj().T@H@U-H)).real/16
                branch=K@rho_root;reference+=branch.T@branch.conj()
                complete+=K.conj().T@K
        assert np.max(abs(complete-np.eye(8)))<1e-14
        assert np.max(abs(reference-rho.T))<1e-14
        assert abs(work-h*_exact_work_coefficient(0,amplitude=s))<2e-14
        # A different stationary density has its own coefficient; geometry is shared.
        rho2=np.eye(8)/8
        assert abs(np.trace(rho2@H))<1e-15


def test_frequency_sum_rule_ratios_and_alias_are_controlled():
    for n in range(7):
        a=2.**(-n)
        for r in (F(1,2),F(2),F(3),F(7,3)):
            row=m.spatial_frequency_ratio(n,r)
            exact=(np.sin(np.pi*a*float(r)/2)/np.sin(np.pi*a/2))**2
            assert float(row['ratio_lower'])-1e-12<=exact<=float(row['ratio_upper'])+1e-12
    assert (np.sin(np.pi)/np.sin(np.pi/2))**2<1e-28


def test_fixed_profile_work_does_not_force_inelastic_probability():
    for b,s in (((1,1,1),1),((2,1,1),1),((1,1,1),2)):
        target=np.pi**2*float(m.profile_geometry(2,b,s)['limiting_work_over_selected_density_pi2_coefficient'])
        previous=1.
        for n in (8,12,16):
            a=2.**(-n);probability=np.sin(np.sqrt(target*a))**2
            assert probability<previous;previous=probability
            assert abs(probability/a-target)<=target*target*a/3
            # The populated energy is 1/a, outside every fixed finite band eventually.
            assert 1/a>100


def test_fixed_anchor_relative_remainder_cannot_hide_creation():
    energy=np.diag([0.,20.]);T=np.array([[0.,np.sqrt(15)],[np.sqrt(15),0.]])
    for alpha in (.1,.5,.9):
        relative=np.diag([0.,alpha*20]);bounded=T-relative
        assert np.linalg.norm(bounded,2)>=np.sqrt(15)
        assert np.linalg.norm(bounded[:,0])==pytest.approx(np.sqrt(15))
    # With a nonzero cross-block the pure form inequality fails at small t.
    v=np.array([1.,1e-6]);assert abs(v@T@v)>.9*(v@energy@v)


def test_actual_member_bounds_are_source_attached_and_unpromoted():
    for n in (0,1,2,8):
        member=m.contract.member(n);bound=m.bind_member(member);m.verify_bound(member,bound)
        r=bound['result']
        assert F(r['uniform_selected_density_lower'])>F(424,10000)
        assert F(r['selected_density_lower'])==F(r['uniform_selected_density_lower'])+F(r['exact_finite_size_density_floor_correction'])
        assert F(r['exact_finite_size_density_floor_correction'])>=F(r['finite_size_density_floor_correction_lower'])
        assert r['physical_mass_table_GeV']==[]
        assert not r['actual_uniform_dimensionless_gap_proved']
        assert not r['inelastic_fixed_energy_probability_proved']
        assert not r['work_form_transferred_to_limiting_probability_measure']
        assert not r['initial_energy_record_preserves_joint_purification']
        assert r['pure_relative_bond_remainder_squared_norm_lower']==15
        a=F(member['source']['spacing']);hstar=F(r['uniform_selected_density_lower'])
        assert F(r['complete_energy_clock_signed_mean_error_upper'])<=F(11,32)*a*a*F(16,9)*hstar/128
        assert r['complete_fine_record_count']==2*member['source']['side']**3*r['complete_two_energy_clock']['clock_plan']['clock_size']**2
        assert json.loads(json.dumps(bound,sort_keys=True,allow_nan=False))==bound
    bad=deepcopy(bound);bad['result']['selected_density_lower']='1'
    with pytest.raises(ValueError):m.verify_bound(member,bad)
    bad_member=deepcopy(member);bad_member['source']['spacing']='1'
    with pytest.raises(ValueError):m.bind_member(bad_member)


def test_profile_domain_rejects_invalid_or_oversized_settings():
    for args in ((0,(2,1,1),1),(2,(1,0,1),1),(2,(1,1),1),(2,(1,1,1),0)):
        with pytest.raises(ValueError):m.profile_geometry(*args)
