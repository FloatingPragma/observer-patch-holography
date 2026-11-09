"""Native bond algebra, refinement envelopes and complete-record error tests."""
from copy import deepcopy
from fractions import Fraction as F
from math import log

import numpy as np
import pytest

from . import cone


def test_actual_native_bond_and_full_fock_norm():
    rows = cone.native_wilson_bonds()
    assert [r['axis'] for r in rows] == [0, 1, 2]
    assert all(r['two_site_one_body_dimension'] == 120 for r in rows)
    assert all(r['full_Fock_bond_norm_times_spacing'] == 30 for r in rows)
    # Independent full exterior-space matrix for one actual Dirac/internal
    # copy. All256 occupations and the fermionic reordering signs are kept.
    gamma = cone.bounds.field_route.source_data.load()['gamma'][0]
    beta = np.asarray(gamma[3], dtype=complex)
    b = -beta @ (np.eye(4)-np.asarray(gamma[0], dtype=complex))/2
    h = np.block([[np.zeros((4, 4)), b], [b.conj().T, np.zeros((4, 4))]])
    H = np.zeros((256, 256), dtype=complex)
    for ket in range(256):
        for q in range(8):
            if not (ket >> q)&1:
                continue
            after = ket ^ (1 << q)
            signq = (-1)**(ket & ((1 << q)-1)).bit_count()
            for p in range(8):
                if (after >> p)&1:
                    continue
                bra = after | (1 << p)
                signp = (-1)**(after & ((1 << p)-1)).bit_count()
                H[bra, ket] += signq*signp*h[p, q]
    assert np.array_equal(H, H.conj().T)
    actual = np.linalg.eigvalsh(H)
    expected = np.repeat([-2., -1., 0., 1., 2.], [16, 64, 96, 64, 16])
    np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=0)
    assert abs(actual[-1]*15-30) < 2e-13


@pytest.mark.parametrize('kind', ['smooth', 'wilson'])
def test_actual_scaled_members_and_unbounded_electric_coefficients(kind):
    config = cone.family(kind, delta='3/2', magnetic_scale='7/11')
    for a in (F(1), F(1, 3), F(1, 1024)):
        slow = cone.spacing_member(config, a, electric_coefficient='0')
        fast = cone.spacing_member(config, a, electric_coefficient=1/a**9)
        assert slow['rate_times_spacing_upper'] == config['weighted_rate_times_spacing_upper']
        assert fast['weighted_interaction_rate_upper'] == slow['weighted_interaction_rate_upper']
        smaller = cone.spacing_member(config, a, magnetic=F(7, 22)/a)
        assert smaller['rate_times_spacing_upper'] < fast['rate_times_spacing_upper']


def test_superlinear_magnetic_scaling_is_not_silently_admitted():
    config = cone.family('wilson', magnetic_scale='1')
    for a in (F(1, 2), F(1, 7)):
        with pytest.raises(ValueError, match='refinement envelope'):
            cone.spacing_member(config, a, magnetic=1/a**2)


@pytest.mark.parametrize('key,value', [
    ('weighted_rate_times_spacing_upper', F(1)),
    ('log_weight_rational_lower', F(100)),
    ('physical_speed_rational_upper', F(0)),
])
def test_false_rate_or_log_certificate_is_rejected(key, value):
    config = deepcopy(cone.family('wilson'))
    config[key] = value
    with pytest.raises(ValueError, match='source derivation'):
        cone.history_bound(config, '1/128', '1/1000', '4')


def test_alignment_safe_cell_growth_for_offset_cubes():
    for a, radius in ((F(1, 7), F(2, 3)), (F(1, 13), F(1, 10)), (F(1, 8), F(0))):
        for center in (F(0), F(1, 101), F(-3, 17)):
            actual = sum(abs(a*n-center) <= radius for n in range(-100, 101))**3
            assert actual <= cone.support_count(radius, a)
        assert cone.support_count(radius, a) <= (2*radius+3)**3/a**3


def test_wilson_spacing_plan_covers_a_continuum_of_smaller_spacings():
    config = cone.family('wilson')
    plan = cone.spacing_plan(config, '1/1000', '4', '1/10', 20, 2)
    assert config['weighted_rate_times_spacing_upper'] == F(18048, 25)
    assert plan['physical_exponential_margin_lower'] == F(1738, 3125)
    assert plan['all_smaller_spacing_upper'] == F(1, 128)
    assert plan['certified_error_bits'] == 41
    g, q, C = (plan[k] for k in ('physical_exponential_margin_lower',
                                'inverse_spacing_prefactor_power',
                                'support_and_apparatus_prefactor'))
    # Non-dyadic samples check the independent real exponential envelope;
    # monotonicity for every smaller real spacing is proved in PROOF.md.
    for a in (F(1, 128), F(1, 129), F(3, 400), F(1, 200), F(1, 1024)):
        assert a <= g/q
        log_error = log(float(C))-q*log(float(a))-float(g/a)
        assert log_error <= -41*log(2)
    assert cone.history_bound(config, F(1, 64), '1/1000', '4', '1/10', 2)['certified_error_bits'] < 20


def test_growing_apparatus_count_and_fixed_duration_are_included():
    config = cone.family('wilson')
    simple = cone.spacing_plan(config, '1/1000', '4', '1/10', 20, 2, 0)
    growing = cone.spacing_plan(config, '1/1000', '4', '1/10', 20, 2, 5)
    assert growing['inverse_spacing_prefactor_power'] == 8
    assert growing['duration'] == simple['duration'] == F(1, 1000)
    assert growing['all_smaller_spacing_upper'] < simple['all_smaller_spacing_upper']
    assert growing['intermediate_local_dilation_count_upper'] == 2/growing['spacing']**5


def test_inside_cone_is_an_uncertified_margin_not_a_false_zero():
    config = cone.family('wilson')
    result = cone.history_bound(config, '1/128', '1', '1')
    assert result['error_binary_exponent_upper'] is None
    assert result['complete_local_history_half_trace_upper'] == {'base': 2, 'exponent': 0}
    with pytest.raises(ValueError, match='positive physical cone margin'):
        cone.spacing_plan(config, '1', '1')


def test_two_implemented_experiment_errors_both_add():
    assert cone.error_composition(F(1, 16), (F(1, 32),), (F(1, 64), F(1, 128))) == F(15, 128)
    assert cone.error_composition(F(1, 2), (F(1, 2),), (F(1, 2),)) == 1
    with pytest.raises(ValueError):
        cone.error_composition(F(2))


def test_domains_and_complete_scientific_scope_survive_export():
    report = cone.certificate()
    assert report['smooth_fixed_time_example']['duration'] == '1'
    assert report['smooth_fixed_time_example']['certified_error_bits'] == 35
    assert report['physical_support_growth_included']
    assert not report['strict_zero_tail_at_finite_spacing']
    assert not report['Lorentz_covariance_continuum_vacuum_or_chirality_proved']
    assert not report['measured_light_speed_or_particle_mass_inferred']
    for args in [('smooth', '0'), ('wilson', '1', '1/100', '1')]:
        with pytest.raises(ValueError):
            cone.family(*args)
    for bad in (0, -1, 0.1, True):
        with pytest.raises(ValueError):
            cone.history_bound(cone.family('wilson'), bad, '1/1000', '4')
