from fractions import Fraction as F
import pytest
from .source import (distances, finite_certificate, minus_log_bounds,
                     reflection_forms)


def test_bond_positivity_distinguishes_valid_probability_from_positive_transfer():
    odd = [[F(1), F(-1)]]
    assert reflection_forms([{1:F(1)}, {0:F(1)}], odd) == (F(1), F(-1))
    lazy = [{0:F(3,4),1:F(1,4)}, {0:F(1,4),1:F(3,4)}]
    assert reflection_forms(lazy, odd) == (F(1), F(1,2))
    row = finite_certificate([(1,), (0,)], lazy, 4)
    assert row['laplacian_gap_lower_exact'] == '2'
    assert row['laplacian_gap_rayleigh_upper_exact'] == '2'
    lo, hi = map(F, (row['event_hamiltonian_gap_lower_exact'],
                    row['event_hamiltonian_gap_upper_exact']))
    assert F(693147180559,10**12) < lo < hi < F(693147180560,10**12)


def test_nonlazy_and_disconnected_sources_rejected():
    with pytest.raises(ValueError): distances([(), ()])
    with pytest.raises(ValueError):
        finite_certificate([(1,), (0,)], [{1:F(1)}, {0:F(1)}], 1)


def test_rational_log_tail_nested_and_identity_exact():
    assert minus_log_bounds(1) == (0, 0)
    a, b = minus_log_bounds(F(2,3), 16)
    c, d = minus_log_bounds(F(2,3), 32)
    assert a < c < d < b
    for value in (0, -1, 2):
        with pytest.raises(ValueError): minus_log_bounds(value)


def test_altered_rate_cannot_reuse_same_transfer_certificate():
    rows = [{0:F(3,4),1:F(1,4)}, {0:F(1,4),1:F(3,4)}]
    with pytest.raises(AssertionError):
        finite_certificate([(1,), (0,)], rows, 5)
