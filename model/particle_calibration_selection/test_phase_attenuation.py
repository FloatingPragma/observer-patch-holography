"""Source comparisons and exact interval controls for phase-selection tests."""
from fractions import Fraction as F

import pytest
import sympy as s

from model.particle_calibration_selection import phase_attenuation as pa
from model.particle_candidate_atlas import attenuation


def test_exact_moments_on_literal_source():
    assert pa.moment_certificate()['literal_three_cycle_replayed']
    r, x, y = '3/5', s.Rational(99, 101), s.Rational(20, 101)
    _, masses = pa.roots_and_masses(x, y, r)
    actual = attenuation.family(r, x+s.I*y)['Y']
    for power in (1, 2, 3):
        assert s.simplify(sum(m**power for m in masses)-s.trace(actual**power)) == 0


@pytest.mark.parametrize('r', ['0', '3/5', '1'])
def test_zero_momentum_full_car_energy_is_phase_flat(r):
    expected = -F(5, 18)*(1+F(r)**2)
    for x, y in ((1, 0), (s.Rational(63, 65), s.Rational(16, 65))):
        energy = pa.free_energy_interval(0, x, y, r)
        assert energy.lo <= expected <= energy.hi
        assert energy.hi-energy.lo < F(1, 10**27)


def test_nonzero_momentum_reverses_the_two_phase_ordering():
    result = pa.finite_momentum_comparisons()
    small, large = result['rows']
    assert F(small['control_minus_zero'][1]) < -F(2, 1000)
    assert F(small['control_minus_endpoint'][1]) < -F(3, 100000)
    assert F(large['control_minus_zero'][0]) > F(5, 1000)
    assert small['interior_minimizer_exists_in_positive_chamber']
    assert result['measured_inputs_used'] == []


def test_phase_conjugation_preserves_the_complete_mass_multiset():
    x, y = s.Rational(63, 65), s.Rational(16, 65)
    _, plus = pa.roots_and_masses(x, y)
    _, minus = pa.roots_and_masses(x, -y)
    assert plus[0] == minus[0]
    assert plus[1] == minus[2]
    assert plus[2] == minus[1]


def test_attenuation_convexity_identity_and_endpoint():
    assert pa.convexity_certificate()['root_slopes_sum'] == '0'
    for p in (F(0), F(1, 256), F(1, 64), F(1)):
        lower = pa.free_energy_interval(p, attenuation='1')
        higher = pa.free_energy_interval(p, attenuation='3/5')
        assert lower.hi < higher.lo


def test_complete_quantum_higgs_sector_and_source_vertex():
    result = pa.quantum_block_certificate()
    assert result['literal_source_full_CAR_states_compared'] == 192
    assert result['complete_zero_charge_block_dimension'] == 307
    assert result['retained_dark_dimension'] == 305
    assert result['full_local_dimension'] == 786432
    assert not result['full_local_ground_identified']
    assert not result['native_controller_or_ground_preparation_selected']


@pytest.mark.parametrize('scalar,mask', [(0, 3), (0, 511), (1, 255), (2, 262143)])
def test_full_vertex_preserves_charge_outside_the_selected_sector(scalar, mask):
    y = s.ImmutableMatrix([[2, s.I, 1], [-s.I, 3, 2*s.I], [1, -2*s.I, 4]])
    charge = mask.bit_count()-2*int(scalar != 0)
    for out_scalar, out_mask in pa.vertex_action(y, scalar, mask):
        assert out_mask.bit_count()-2*int(out_scalar != 0) == charge


def test_complete_sector_contains_unreached_dark_states():
    y = s.diag(1, 2, 3)
    block, basis = pa.zero_charge_block(y, pair_cost=5)
    assert len(basis) == 307
    assert block[0, 0] == 0
    assert block[1:, 1:] == 5*s.eye(306)
    bright = block[:, 0]
    assert (bright.H*bright)[0] == 4*(1+4+9)
    assert sum(value != 0 for value in bright) == 12


@pytest.mark.parametrize('kwargs', [dict(real=1, imag=1), dict(attenuation='-1/10'),
                                    dict(attenuation='11/10'), dict(attenuation=.5)])
def test_invalid_coherent_controls_fail_closed(kwargs):
    with pytest.raises(ValueError):
        pa.roots_and_masses(**kwargs)
