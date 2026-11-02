from fractions import Fraction as F
import pytest
from . import coefficients as c, envelope


def test_contractivity_is_retained_even_if_midpoint_norm_exceeds_one():
    eps = F(1, 1000)
    # The exact admitted vector is (1,i)/2, whereas these midpoints overshoot.
    rectangles = [((F(1, 2), F(1, 2)+eps/8), (F(), F())),
                  ((F(), F()), (F(1, 2), F(1, 2)+eps/8))]
    result = c.contractive_coefficients(rectangles, eps)
    assert result['l1_norm_upper'] <= 1
    assert result['l1_error_upper'] < eps
    exact_error = abs(result['coefficients'][0][0]-F(1, 2))+abs(result['coefficients'][1][1]-F(1, 2))
    assert exact_error <= result['l1_error_upper']


def test_excess_uncertainty_requires_refinement():
    with pytest.raises(ValueError, match='refine'):
        c.contractive_coefficients([((F(), F(1)), (F(), F()))], F(1, 100))


def test_exact_envelope_spatial_interference_has_the_correct_phase():
    plus = envelope.coefficient_rectangle((1, 0, 0), F(1, 2), 80)
    minus = envelope.coefficient_rectangle((-1, 0, 0), F(1, 2), 80)
    assert plus[0] == minus[0]
    assert plus[1] == (-minus[1][1], -minus[1][0])
    assert plus[0][0] > 0 and plus[1][0] > 0
    result = c.contractive_coefficients([plus, minus], F(1, 10**20), 80)
    assert result['l1_error_upper'] <= F(1, 10**20)


def test_rational_square_root_interval_is_certified_by_squaring():
    for value in (F(), F(4), F(2), F(7, 29), F(10**30)):
        lo, hi = c.sqrt_interval(value, 20)
        assert 0 <= lo <= hi and lo*lo <= value <= hi*hi
        assert hi-lo <= F(1, 1 << 20)
