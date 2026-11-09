"""Source-derived finite-horizon and all-horizon capacity-response checks."""
import pytest
import sympy as s

from model.common_source import native_family_record_memory as memory
from model.particle_calibration_selection import capacity_response as c


def test_complete_source_and_reachable_quotient():
    row = c.source()
    owner = memory.architecture()
    assert row['T12'] == owner['T']
    assert row['allocated_raw_labels'] == 128
    assert row['positive_raw_labels'] == 24
    assert row['old_memory_inputs'] == 12
    assert row['T'].shape == (3, 3)
    assert all(sum(row['E12'][:, j]) == s.Rational(1, 12) for j in range(12))
    assert all(sum(row['T12'][:, j]) == 1 for j in range(12))
    assert all(value >= 0 for value in row['T12']-row['E12'])


def test_exact_short_horizons_and_capacity_count():
    assert c.weights(1) == s.Matrix([s.Rational(1, 72), s.Rational(5, 72), 0])
    assert c.weights(2) == s.Matrix([s.Rational(41, 2592), s.Rational(695, 10368), s.Rational(5, 10368)])
    assert c.coefficient(1).rank() == 2
    assert c.coefficient(2).rank() == 3
    for n in (1, 2, 4, 8):
        budget = c.archive_budget(n)
        assert budget['positive_raw_histories'] == 24**n
        assert budget['zero_raw_histories'] + budget['positive_raw_histories'] == 128**n
        assert budget['controller_address_dimension'] == 12
        assert s.trace(c.coefficient(n)) == s.Rational(1, 12)


def test_first_derivative_from_full_twelve_address_channel():
    # Independent response differentiation in all old-memory blocks. This
    # does not use the three-ray running-count recurrence to propagate data.
    row, owner = c.source(), memory.architecture()
    rho = owner['projectors'][2]
    blocks, derivative = [s.zeros(3) for _ in range(12)], [s.zeros(3) for _ in range(12)]
    blocks[2] = rho
    for n in (1, 2, 3):
        next_blocks, next_derivative = [], []
        for new in range(12):
            z, dz = s.zeros(3), s.zeros(3)
            pi = owner['projectors'][new]
            for old in range(12):
                z += row['T12'][new, old]*blocks[old]
                dz += row['T12'][new, old]*derivative[old]
                dz -= s.I*row['E12'][new, old]*(pi*blocks[old]-blocks[old]*pi)
            next_blocks.append(c._clean(z))
            next_derivative.append(c._clean(dz))
        blocks, derivative = next_blocks, next_derivative
        actual = c._clean(sum(derivative, s.zeros(3))/n)
        y = c.coefficient(n)
        assert actual == c._clean(-s.I*(y*rho-rho*y))


def test_all_horizon_invariant_bounds_and_strict_spectral_constant():
    row, bounds = c.source(), c.uniform_spectral_bounds()
    p, a, b = row['p'], row['a'], row['b']
    assert a >= b > 0
    assert min(row['T'][0, :]) >= 1-p
    assert b*a+a*(a*b/(1-a)) == a*b/(1-a)
    assert bounds['smallest_eigenvalue_strict_lower'] == s.Rational(11, 103392)
    assert bounds['middle_to_smallest_strict_upper'] == s.Rational(1171058, 6633)
    assert bounds['largest_to_smallest_strict_upper'] == s.Rational(8616, 11)
    lower = sum((bounds['lower_ray_weights'][i]*row['projectors'][i] for i in range(3)), s.zeros(3))
    epsilon = bounds['smallest_eigenvalue_strict_lower']
    difference = c._clean(lower-epsilon*s.eye(3))
    assert all(s.simplify(difference[:i, :i].det()).is_positive for i in range(1, 4))
    for n in (2, 3, 4, 8, 19, 64):
        w = c.weights(n)
        assert all(w[i] >= bounds['lower_ray_weights'][i] for i in range(3))
        assert w[0]+w[2] <= bounds['middle_eigenvalue_upper']


def test_transient_convergence_bound_keeps_stationary_shape():
    row = c.source()
    stationary = row['E']*row['stationary']
    for n in (1, 2, 5, 12, 64):
        assert sum(abs(v) for v in c.weights(n)-stationary) <= c.finite_coefficient_error_bound(n)
    assert c.finite_coefficient_error_bound(64) < c.finite_coefficient_error_bound(12)
    assert s.trace(c.stationary_coefficient()) == row['p']


def test_complete_finite_channel_is_mixed_without_postselection():
    for n in (1, 2):
        rho = c.complete_unread_density(n, s.Rational(1, 10))
        loss = s.simplify(1-s.trace(rho*rho))
        assert s.trace(rho) == 1
        assert c._clean(rho-rho.H) == s.zeros(3)
        assert loss >= c.purity_loss_floor(n, s.Rational(1, 10)) > 0


def test_comparison_cannot_change_source_and_rejects_invalid_inputs():
    frozen = c.uniform_spectral_bounds()
    result = c.comparison('207', '3500')
    assert result['middle_ratio_above_all_horizon_bound']
    assert result['largest_ratio_above_all_horizon_bound']
    assert c.uniform_spectral_bounds() == frozen
    assert c.comparison('206.7682827', '3477.23')['middle_ratio_above_all_horizon_bound']
    for bad in (0, -1, 1.0, True):
        with pytest.raises(ValueError):
            c.weights(bad)
    for bad in (0, 0.1, s.sqrt(2)):
        with pytest.raises(ValueError):
            c.purity_loss_floor(2, bad)
    with pytest.raises(ValueError):
        c.comparison('10', '2')
