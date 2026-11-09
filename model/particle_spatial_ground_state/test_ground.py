"""Independent checks of source bounds, cooling and exact native stationarity."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
import scipy.linalg as la
import sympy as s

from model.particle_field_state_response import state_response as source
from model.particle_spatial_ground_state import ground as g


def test_spatial_dimensions_and_actual_trace_metric():
    result=g.spatial_source_bounds(side=2,link_degree=3,electric='105')
    assert result['CAR_modes']==480
    assert result['dimension_bits_upper']==480+9*3*24
    assert result['dimension_factored_upper']=={
        'Fock_base':2,'Fock_power':480,'link_base':451,'link_power':72}
    assert result['electric_norm_upper']==F(29,10)*105*24*9
    assert result['kinetic_norm_upper']==360*8
    assert result['spectral_width_upper']==result['spectral_upper']-result['spectral_lower']
    zero=g.spatial_source_bounds(link_degree=0)
    assert zero['electric_norm_upper']==0
    assert zero['dimension_bits_upper']==zero['CAR_modes']


def test_smooth_bound_has_own_coefficients_and_no_large_allocation():
    p=g.spatial_plan(side=4,spacing='1/2',link_degree=2,regulator='smooth',
                     smoothing='1/8',target_excess='1/64')
    assert p['source']['kinetic_norm_upper']==60*64*F(17,8)*2
    assert not p['source']['large_matrix_allocated']
    assert p['cooling']['mean_ground_energy_excess_upper']<=F(1,64)
    assert p['cooling']['state_exactly_stationary']
    assert not p['cooling']['spectral_gap_assumed']


def test_gap_free_plan_succeeds_for_growing_volume():
    for side in (2,4,8):
        source_bounds=g.spatial_source_bounds(side=side,electric='105')
        p=g.cooling_plan(source_bounds['dimension_bits_upper'],
                         source_bounds['spectral_lower'],source_bounds['spectral_upper'],F(1,side))
        assert p['mean_ground_energy_excess_upper']<=F(1,side)
        assert p['dyadic_scale']&(p['dyadic_scale']-1)==0
        assert p['scaled_argument_upper']<=F(1,2)
        assert p['polynomial_degree']==p['Taylor_order']*p['dyadic_scale']
        assert not p['efficient_in_volume_claimed']


def test_scaled_relative_error_matches_independent_two_level_diagonalization():
    h=np.array([[0.,1.],[1.,1.]])
    energies,basis=la.eigh(h)
    b=g.polynomial_bound('8',-1,2,8)
    x=8*(2-energies)/(2*b['dyadic_scale'])
    factors=sum(x**j/math.factorial(j) for j in range(9))
    values=factors**b['dyadic_scale']
    exact=np.exp(8*(2-energies)/2)
    relative=np.max(abs(values/exact-1))
    assert relative<=float(b['relative_amplitude_error_upper'])
    rho=(basis*(values**2/sum(values**2)))@basis.T
    gibbs=(basis*(exact**2/sum(exact**2)))@basis.T
    assert la.norm(rho-gibbs,ord='nuc')/2<=float(b['canonical_state_half_trace_error_upper'])


def test_gap_free_entropy_bound_does_not_require_uniform_gap():
    for gap in (1.,1e-5,1e-12,0.):
        energies=np.array([0.,gap,2.,3.]); beta=7.
        weights=np.exp(-beta*energies);weights/=sum(weights)
        assert weights@energies<=math.log(4)/beta


def test_native_exact_blocks_retain_every_state_and_commute_with_same_H():
    p=g.prepared_control(); phases,h=source.real_computation_chart()
    assert sorted(i for indices,_ in p['amplitude_blocks'] for i in indices)==list(range(208))
    assert {((i//13).bit_count()) for indices,_ in p['amplitude_blocks'] for i in indices}==set(range(5))
    assert p['computational_basis_phases']==phases
    for indices,q in p['amplitude_blocks']:
        hb=h.extract(indices,indices)
        assert q==q.T
        assert q*hb==hb*q
    assert p['mean_energy']<s.Rational(-2)
    assert not p['joint_purification_stationary']


def test_native_exact_energy_and_gibbs_comparison():
    p=g.prepared_control()
    h=np.asarray(source.exact_control()['H'],dtype=complex)
    eigenvalues=la.eigvalsh(h)
    weights=np.exp(-4*(eigenvalues-eigenvalues[0]));weights/=sum(weights)
    exact_gibbs_mean=weights@eigenvalues
    assert abs(float(p['mean_energy'])-exact_gibbs_mean)<=26*float(p['preparation_trace_error_upper'])
    summary=g.native_control_summary()
    lo,hi=summary['mean_energy_interval']
    assert lo<=F(p['mean_energy'])<=hi
    assert 0<summary['mean_ground_energy_excess_upper']<F(157,1000)
    assert summary['ground_energy_interval'][0]<=eigenvalues[0]<=summary['ground_energy_interval'][1]


def test_stationary_ground_inequality_for_noncommuting_native_probe():
    p=g.prepared_control(); _,h=source.real_computation_chart()
    indices,q=next((i,q) for i,q in p['amplitude_blocks'] if len(i)==12)
    hb=h.extract(indices,indices); r=q*q/p['norm_squared']
    # A=|target><source| in one actual reducing block, ||A||=1.
    a=s.zeros(12);a[1,0]=1
    local=s.trace(r*a.T*(hb*a-a*hb))
    eps=F(p['mean_energy'])-g.native_control_summary()['ground_energy_interval'][0]
    assert local.is_real and local>=-_as_sym(eps)
    assert hb*a!=a*hb


def _as_sym(value):return s.Rational(value.numerator,value.denominator)


@pytest.mark.parametrize('kwargs',[{'beta':0},{'order':-1},{'scale':3},{'scale':1},
                                   {'beta':0.5},{'spectral_upper':-5}])
def test_invalid_cooling_data_rejected(kwargs):
    args={'beta':'4','spectral_lower':-4,'spectral_upper':22};args.update(kwargs)
    with pytest.raises(ValueError):g.polynomial_bound(**args)


def test_dyadic_rounding_is_outward_for_both_signs():
    for q in (F(-9,7),F(2,3),F(0),F(2)):
        lo,hi=g.outward_dyadic(q,10)
        assert lo<=q<=hi and hi-lo==F(1,1024)
