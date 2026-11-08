"""Independent full-source spectral controls for uniform original energy tails."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import composite, uv_tail as u


@pytest.mark.parametrize('sites', [2, 3])
def test_full_source_spectral_measures_and_all_moment_bounds(sites):
    row = u.finite_control(sites)
    data = composite.source((sites, 1, 1), row['c'], row['U'])
    k = data['H']/data['scale']
    values, vectors = np.linalg.eigh(k)
    ground = vectors[:, abs(values-values[0]) < 1e-10]
    assert ground.shape[1] == row['normalized_full_ground_rank']
    rho = ground@ground.T/ground.shape[1]
    energy = np.maximum(0, values-values[0])
    # A complex nonlocal packet exercises both flavour and intersite Gram terms.
    f = np.arange(1, 2*sites+1)+1j*np.arange(2*sites, 0, -1)
    f /= np.linalg.norm(f)
    operator = sum(z*a for z, a in zip(f, data['a'], strict=True))
    probabilities = []
    for a in (operator, operator.conj().T):
        amplitudes = vectors.conj().T@a@ground
        p = np.sum(abs(amplitudes)**2, axis=1)/ground.shape[1]
        probabilities.append(p)
    p = sum(probabilities)
    assert abs(p.sum()-1) < 1e-13
    degree = 1 if sites == 2 else 2
    for n in range(1, 7):
        actual = float(p@energy**(2*n))
        bound = float(u.even_moment_bound(row['c'], row['U'], degree, n))
        assert actual <= bound*(1+1e-12)
    for cutoff in (F(1, 2), F(2), F(4), F(16), F(32)):
        tail = p[energy > float(cutoff)].sum()
        assert tail <= float(u.tail_budget(row['c'], row['U'], degree, cutoff)
                             ['original_total_charge_tail_upper'])+1e-13
        assert tail <= float(u.exponential_tail(row['c'], row['U'], degree, cutoff))+1e-13
    _, b, _ = u.parameters(row['c'], row['U'], degree)
    slope = 1/(2*b)
    assert p@np.exp(float(slope)*energy) < float(u.exponential_moment(
        row['c'], row['U'], degree, slope))
    assert abs(np.trace(rho)-1) < 1e-13


def test_cubic_ball_cardinality_matches_actual_lattice_points():
    for r in range(7):
        count = sum(abs(x)+abs(y)+abs(z) <= r
                    for x in range(-r, r+1) for y in range(-r, r+1)
                    for z in range(-r, r+1))
        assert u.ball(r) == count


def test_default_alias_budget_and_ballistic_clock_scaling():
    p = F(25, 100352)
    row = u.tail_budget(F(1, 8), 1, 6, 20)
    assert row['moment_order'] == 18
    assert 0 < row['original_total_charge_tail_upper'] < p/16
    assert u.even_moment_bound(F(1, 8), 1, 6, 1) == F(25, 64)
    # Scaling action and energy together preserves the original tail budget.
    for scale in (F(1, 5), F(7), F(10000)):
        other = u.tail_budget(scale/8, scale, 6, 20*scale)
        assert other['original_total_charge_tail_upper'] == row['original_total_charge_tail_upper']
    assert u.tail_budget(1250, 10000, 6, 20)['original_total_charge_tail_upper'] == 1


def test_zero_action_derivatives_and_atomic_tail():
    assert u.local_derivative_bound(1, 0, 0, 6) == 0
    assert u.exponential_tail(1, 0, 0, 1) == 0
    assert u.exponential_moment(1, 0, 0, 100) == 1
    assert u.tail_budget(1, 0, 0, 1)['original_total_charge_tail_upper'] == 0
    # A disconnected centered atom has all original weight at U/2.
    assert u.tail_budget(1, 2, 0, F(1, 2))['original_total_charge_tail_upper'] == 1
    assert 0 < u.tail_budget(1, 2, 0, 4)['original_total_charge_tail_upper'] < 1


@pytest.mark.parametrize('args', [(True, 1, 6), (0, 1, 6), (1, -1, 6),
                                  (1., 1, 6), (1, 1, True), (1, 1, 7)])
def test_reject_unsupported_source_domains(args):
    with pytest.raises(ValueError):
        u.parameters(*args)


def test_strict_energy_and_analytic_strip_boundaries():
    for cutoff in (0, -1, 1., True):
        with pytest.raises(ValueError):
            u.tail_budget(F(1, 8), 1, 6, cutoff)
    for order in (0, 257, True):
        with pytest.raises(ValueError):
            u.tail_budget(F(1, 8), 1, 6, 20, order)
    for slope in (F(1, 2), F(1), -1, True):
        with pytest.raises(ValueError):
            u.exponential_moment(F(1, 8), 1, 6, slope)


def test_exact_psd_control_rejects_negative_and_zero_row_counterexamples():
    for matrix in ([[1, 2], [2, 1]], [[0, 1], [1, 0]], [[-1]]):
        with pytest.raises(ArithmeticError):
            u._positive(matrix)
    u._positive([[F(1, 3), F(2, 3)], [F(2, 3), F(4, 3)]])


def test_cached_controls_detach_all_rational_and_container_slots():
    row = u.finite_control(2)
    expected = u.finite_control(2)
    row['c']._numerator = 100
    row['moments'][0]['uniform_operator_bound']._denominator = 101
    row['moments'].clear()
    assert u.finite_control(2) == expected


def test_scope_never_promotes_moments_to_particles():
    row = u.certificate()
    assert row['separate_charge_measures_determined_by_all_their_moments']
    assert not row['finite_moment_list_identifies_a_pole']
    assert not row['positive_shrinking_window_weight_proved']
    assert row['physical_species'] == row['physical_masses'] == []
