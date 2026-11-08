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
