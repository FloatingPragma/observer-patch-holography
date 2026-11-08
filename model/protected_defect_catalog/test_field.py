"""Operational odd preparation and quantitative domain controls."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import charge, field


def test_complete_odd_even_source_instrument():
    op = charge.operators(); indices=np.arange(256)
    parity=np.eye(256)[indices ^ 255]
    # Two different amplitudes and sites, not a scalar function of parity.
    angles=np.array([(1-2*(w & 1))/3+(1-2*((w >> 1) & 1))/7 for w in range(256)])
    odd=np.diag(np.sin(angles)); even=np.diag(np.cos(angles))
    assert np.max(np.abs(odd.T@odd+even.T@even-np.eye(256))) < 4e-16
    assert np.array_equal(parity@odd,-odd@parity)
    assert np.array_equal(parity@even,even@parity)
    assert np.array_equal(parity@op['H'],op['H']@parity)


def test_unitary_variance_is_bounded_by_odd_preparation_weight():
    # Independent arbitrary spectral probability laws check the inequality
    # 1-|<exp(2iA)>|^2 <= 4<sin(A)^2> used in the uniform proof.
    rng=np.random.default_rng(7301)
    for _ in range(20):
        angles=rng.uniform(-5,5,9); probabilities=rng.random(9)
        probabilities/=probabilities.sum()
        characteristic=np.dot(probabilities,np.exp(2j*angles))
        weight=np.dot(probabilities,np.sin(angles)**2)
        assert 1-abs(characteristic)**2 <= 4*weight+1e-14


def test_uniform_probability_retains_normalization_cost():
    result=field.preparation_bounds()
    assert 0 < F(result['liminf_original_odd_preparation_probability']) < 1
    assert 0 < F(result['sufficiently_small_scale_filtered_probability_lower']) < F(result['liminf_odd_filtered_original_probability'])
    assert not result['effective_scale_threshold_supplied']
    assert not result['physical_mass_shell_or_species_identified']


@pytest.mark.parametrize('args',[(0,1,1),(1,1,1),(F(1,4),0,1),(F(1,4),1,0)])
def test_invalid_field_filter_domain(args):
    with pytest.raises(ValueError):field.preparation_bounds(*args)


def test_configuration_rotation_and_charge_conjugation_on_source():
    result=field.configuration_charge()
    assert result['source_transverse_charge_orientation_checked']
    assert result['source_combined_translation_rotation_checked']
    assert not result['second_canonical_momentum_pair_attached']
    assert not result['physical_antiparticle_or_Lorentz_spin_identified']


def test_configuration_commutator_bound_on_independent_tensor_spins():
    y=np.array([[0,-1j],[1j,0]])/2; z=np.diag([.5,-.5])
    one=np.eye(2)
    def site(a,x):
        out=np.ones((1,1))
        for j in range(4):out=np.kron(out,a if j==x else one)
        return out
    f=(F(1,2),F(-2,3),F(1,5),F(4,7))
    g=(F(1,3),F(1,2),F(-3,2),F(1,7))
    for p in (F(1,2),F(1,4),F(1,8)):
        a=sum((float(p*p*f[x])*(-1)**x*site(y,x) for x in range(4)))
        b=sum((float(p*p*g[x])*(-1)**x*site(z,x) for x in range(4)))
        bound=field.configuration_commutator_bound(p,sum(abs(u*v) for u,v in zip(f,g)))
        assert np.isclose(np.linalg.norm(a@b-b@a,2),float(bound),rtol=1e-13,atol=1e-17)


@pytest.mark.parametrize('args',[(0,1),(-1,1),(1,-1)])
def test_invalid_configuration_domain(args):
    with pytest.raises(ValueError):field.configuration_commutator_bound(*args)
