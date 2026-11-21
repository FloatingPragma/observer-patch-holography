"""Independent synthetic tests for exact radar enclosures; no contact source execution."""

from fractions import Fraction
import json

import mpmath as mp
import pytest

from . import certified_ranges as certified


def _inside(interval, value, scale):
    assert mp.mpf(interval[0]) / scale <= value <= mp.mpf(interval[1]) / scale


def test_directed_real_operations_enclose_exact_rational_samples():
    scale = 16
    boxes = [(-19, -5), (-7, 13), (0, 0), (3, 21)]
    for left in boxes:
        for right in boxes:
            product = certified._mul(left, right, scale)
            difference = certified._sub(left, right)
            for a in (Fraction(left[0], scale), Fraction(sum(left), 2 * scale), Fraction(left[1], scale)):
                for b in (Fraction(right[0], scale), Fraction(sum(right), 2 * scale), Fraction(right[1], scale)):
                    assert Fraction(product[0], scale) <= a * b <= Fraction(product[1], scale)
                    assert Fraction(difference[0], scale) <= a - b <= Fraction(difference[1], scale)
        squared = certified._square(left, scale)
        for point in range(left[0], left[1] + 1):
            assert Fraction(squared[0], scale) <= Fraction(point, scale) ** 2 <= Fraction(squared[1], scale)
    ratio = certified._ratio((3, 17), (7, 23), scale)
    assert Fraction(ratio[0], scale) <= Fraction(3, 23)
    assert Fraction(ratio[1], scale) >= Fraction(17, 7)
    with pytest.raises(ValueError):
        certified._ratio((1, 2), (0, 3), scale)


def test_complex_products_enclose_independent_fraction_calculation():
    scale = 32
    left = ((-17, 21), (4, 13))
    right = ((3, 9), (-11, -5))
    actual = certified._complex_mul(left, right, scale)
    for ar in left[0]:
        for ai in left[1]:
            for br in right[0]:
                for bi in right[1]:
                    real, imaginary = Fraction(ar * br - ai * bi, scale ** 2), Fraction(ar * bi + ai * br, scale ** 2)
                    assert Fraction(actual[0][0], scale) <= real <= Fraction(actual[0][1], scale)
                    assert Fraction(actual[1][0], scale) <= imaginary <= Fraction(actual[1][1], scale)


def test_integer_moments_have_no_similarity_roundoff():
    chart = [[0, 2, 0], [1, 0, 1], [0, 2, 0]]
    moments = certified.integer_moments(chart, [0, 2], 4)
    assert moments[0] == [[1, 0], [0, 1]]
    assert moments[1] == [[0, 0], [0, 0]]
    assert moments[2] == [[2, 2], [2, 2]]
    assert moments[4] == [[8, 8], [8, 8]]
    with pytest.raises(ValueError):
        certified.enclose_propagators(chart, [1, 1, 1], [0, 2], norm_bound=2)
    with pytest.raises(ValueError):
        certified.enclose_propagators(chart, [1, 2, 1], [0, 1], norm_bound=2)
    with pytest.raises(ValueError):
        certified.enclose_propagators(chart, [1, 2, 1], [0, 2], norm_bound=1)


def test_nonnormal_chart_matches_analytic_three_mode_propagator():
    # H has adjacent couplings 1/sqrt(2), with eigenvalues -1, 0 and 1.
    chart = [[0, 2, 0], [1, 0, 1], [0, 2, 0]]
    kernel = certified.enclose_propagators(chart, [1, 2, 1], [0, 2], degree=120, norm_bound=2)
    scale = kernel['dyadic_scale']
    with mp.workdps(100):
        for n, table in enumerate(kernel['amplitudes']):
            diagonal = (1 + mp.cos(mp.mpf(n) / 2)) / 2
            off_diagonal = (mp.cos(mp.mpf(n) / 2) - 1) / 2
            for target in (0, 1):
                for source in (0, 1):
                    _inside(table[target][source][0], diagonal if target == source else off_diagonal, scale)
                    _inside(table[target][source][1], mp.mpf(0), scale)
                    assert table[target][source][0][1] - table[target][source][0][0] < 200


def test_two_mode_first_detection_matches_geometric_hit_law():
    kernel = certified.enclose_propagators([[0, 2], [2, 0]], [1, 1], [0, 1], degree=120, norm_bound=1)
    result = certified.renewal_enclosures(kernel)
    scale = result['dyadic_scale']
    with mp.workdps(100):
        hit = mp.sin(mp.mpf(1) / 2) ** 2
        miss = 1 - hit
        probabilities = [hit * miss ** n for n in range(12)]
        success = 1 - miss ** 12
        mean = mp.fsum((n + 1) * p for n, p in enumerate(probabilities)) / success
        for source, target in ((0, 1), (1, 0)):
            for n, p in enumerate(probabilities):
                _inside(result['first_hit_probability'][source][target][n], p, scale)
            _inside(result['leg_success_probability'][source][target], success, scale)
            _inside(result['leg_timeout_probability'][source][target], miss ** 12, scale)
            _inside(result['leg_conditional_mean_ticks'][source][target], mean, scale)
            _inside(result['distance'][source][target], mean / 2, scale)
            _inside(result['round_trip_success_probability'][source][target], success ** 2, scale)
            assert result['distance_defined'][source][target]
        exported = certified.export_grid(result)
        _inside((exported['distance_lower'][0][1], exported['distance_upper'][0][1]), mean / 2,
                exported['scale'])
        assert exported['distance_upper'][0][1] - exported['distance_lower'][0][1] <= 1


def test_disconnected_pair_retains_exact_timeout_and_undefined_mean():
    kernel = certified.enclose_propagators([[2, 0], [0, 4]], [1, 1], [0, 1], degree=60, norm_bound=2)
    result = certified.renewal_enclosures(kernel)
    scale = result['dyadic_scale']
    assert result['leg_success_probability'][0][1] == (0, 0)
    assert result['leg_timeout_probability'][0][1] == (scale, scale)
    assert not result['leg_positive_success_proved'][0][1]
    assert not result['distance_defined'][0][1]
    assert result['leg_conditional_mean_ticks'][0][1] == (scale, 12 * scale)
    assert result['distance'][0][1] == (scale // 2, 6 * scale)


def test_integral_remainder_covers_a_deliberately_short_series():
    kernel = certified.enclose_propagators([[2]], [1], [0], ticks=3, dt=Fraction(3, 2),
                                          degree=2, norm_bound=1)
    with mp.workdps(100):
        for n, table in enumerate(kernel['amplitudes']):
            value = mp.exp(-1j * mp.mpf(3 * n) / 2)
            _inside(table[0][0][0], value.real, kernel['dyadic_scale'])
            _inside(table[0][0][1], value.imag, kernel['dyadic_scale'])
    assert Fraction(*kernel['maximum_taylor_remainder']) == Fraction(9, 2) ** 3 / 6


def test_all_exported_numeric_endpoints_are_exact_integers():
    kernel = certified.enclose_propagators([[0, 2], [2, 0]], [1, 1], [0, 1], ticks=3,
                                          degree=60, norm_bound=1)
    exported = certified.export_grid(certified.renewal_enclosures(kernel))
    def inspect(value):
        if isinstance(value, dict):
            for child in value.values():
                inspect(child)
        elif isinstance(value, list):
            for child in value:
                inspect(child)
        else:
            assert type(value) in (int, bool)
    inspect(exported)
    assert json.loads(json.dumps(exported)) == exported
    scale = exported['scale']
    assert exported['leg_success_probability_lower'][0][1] + exported['leg_timeout_probability_upper'][0][1] == scale
    with pytest.raises(ValueError):
        certified.enclose_propagators([[0.0]], [1], [0], ticks=1)
    with pytest.raises(ValueError):
        certified.enclose_propagators([[0]], [1], [0], ticks=1, dt=.5)
