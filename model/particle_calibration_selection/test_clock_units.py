from fractions import Fraction as Q
import pytest
import sympy as s
from . import clock_units as c


def test_actual_complete_source_and_ratios():
    value = c.certificate()
    assert value['complete_input_columns'] == 66
    assert value['raw_seam_coin_outcomes_per_input'] == 60
    for key, interval in value['certified_same_clock_ratios'].items():
        lo, hi = Q(interval['lower']), Q(interval['upper'])
        assert hi-lo < Q(1, 10**22)
        if key == 'middle_over_slow':
            assert 2 < lo < hi < 3
        else:
            assert 2 < lo < hi < 3
    assert not value['physical_particle_identity_proved']
    assert not value['physical_GeV_conversion_derived']


def test_dropping_no_movement_records_breaks_source_dictionary():
    p, b, laplacian = c.source_matrices()
    jumps = s.Matrix(p)
    for i in range(66):
        jumps[i, i] = 0
        total = sum(jumps[i, j] for j in range(66))
        for j in range(66):
            jumps[i, j] /= total
    with pytest.raises(ValueError):
        c.verify_source(jumps, b, laplacian)


def test_log_bounds_and_ratio_error_without_float_premises():
    lo, hi = c.neg_log_bounds(Q(1, 2), Q(1, 2), 40)
    assert Q(693147180559945309417, 10**21) < lo < hi < Q(693147180559945309418, 10**21)
    bound = c.rate_ratio_error(Q(3), Q(2), Q(1, 100), Q(1, 50))
    for ea in (-Q(1, 100), Q(1, 100)):
        for eb in (-Q(1, 50), Q(1, 50)):
            assert abs((3+ea)/(2+eb)-Q(3, 2)) <= bound
    for args in ((Q(0), Q(1)), (0.5, Q(1))):
        with pytest.raises(ValueError):
            c.neg_log_bounds(*args)


def test_independent_clock_rescaling_cannot_cancel():
    # One common factor cancels; two independent factors change the ratio.
    a, b = Q(7, 13), Q(5, 19)
    assert (Q(3)*a)/(Q(3)*b) == a/b
    assert (Q(3)*a)/(Q(2)*b) != a/b
