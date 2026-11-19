"""Synthetic exact metric controls, without reading the scientific distances."""

from fractions import Fraction
from itertools import permutations

import numpy as np
import pytest

from .metric_certificate import (THRESHOLDS, certify, determinant, five_point_box,
                                 gram_integer, rational)


def _permutation_determinant(a):
    result = 0
    for p in permutations(range(len(a))):
        sign = (-1) ** sum(p[i] > p[j] for i in range(len(a)) for j in range(i + 1, len(a)))
        term = sign
        for i, j in enumerate(p):
            term *= int(a[i][j])
        result += term
    return result


def test_bareiss_matches_independent_permutation_formula():
    rng = np.random.default_rng(80)
    for n in range(1, 5):
        for _ in range(40):
            a = rng.integers(-5, 6, (n, n)).tolist()
            assert determinant(a) == _permutation_determinant(a)
    assert determinant([[0, 2], [3, 4]]) == -6
    assert determinant([[1, 2], [2, 4]]) == 0


def test_regular_simplex_excluded_with_quantitative_margin():
    result = five_point_box([Fraction(9, 10)] * 10, [Fraction(11, 10)] * 10)
    assert result["excludes_3d_fit"]
    assert result["corners_checked"] == 1024
    assert all(int(x) > 0 for x in result["minimum_integral_leading_minors"])
    assert Fraction(result["uniform_gram_eigenvalue_lower_bound"]) > 0
    # The exact regular 4-simplex has anchored eigenvalues 1/2,1/2,1/2,5/2.
    fixed = five_point_box([1] * 10, [1] * 10)
    assert 0 < Fraction(fixed["uniform_gram_eigenvalue_lower_bound"]) <= Fraction(1, 2)


def test_inner_segment_certifies_only_a_subset():
    result = five_point_box([Fraction(1, 4)] * 10, [Fraction(9, 4)] * 10)
    assert result["status"] == "subset_compatible"
    assert result["subset_fit_witness"]["kind"] == "positive_definite_midpoint_to_negative_direction_corner"
    assert not result["excludes_3d_fit"]


def _rational_three_dimensional_distances():
    # The Euler brick has face diagonals 125, 244 and 267. All distances
    # among these six points are integers, and the full set spans R^3.
    points = np.array([[0, 0, 0], [44, 0, 0], [-44, 0, 0],
                       [0, 117, 0], [0, -117, 0], [0, 0, 240]], dtype=np.int64)
    squares = ((points[:, None] - points[None, :]) ** 2).sum(axis=-1)
    distances = np.rint(np.sqrt(squares)).astype(np.int64)
    assert np.array_equal(distances ** 2, squares)
    return distances


def test_nondegenerate_euclidean_data_are_never_excluded():
    result = certify(_rational_three_dimensional_distances(), numerical_guard=0, thresholds=(0, Fraction(1, 20)))
    for row in result["thresholds"]:
        assert not row["global_3d_fit_excluded"]
        assert row["compatible_subsets"] == 6
    assert "need not share" in result["no_global_positive_claim"]


def test_unknown_data_remain_unresolved_and_envelopes_are_usable():
    unknown = np.full((6, 6), np.nan)
    np.fill_diagonal(unknown, 0)
    result = certify(unknown, thresholds=(Fraction(1, 20),))
    assert result["thresholds"][0]["unresolved_subsets"] == 6
    assert not result["thresholds"][0]["global_3d_fit_excluded"]
    lower, upper = .99 * (1 - np.eye(6)), 1.01 * (1 - np.eye(6))
    bands = certify(unknown, lower, upper, thresholds=(0,), confidence_failure_probability=1e-6)
    assert bands["thresholds"][0]["global_3d_fit_excluded"]
    assert bands["confidence_failure_probability"] == 1e-6
    broad = certify(unknown, .1 * (1 - np.eye(6)), 10 * (1 - np.eye(6)), thresholds=(Fraction(1, 20),))
    assert not broad["thresholds"][0]["global_3d_fit_excluded"]


def test_outward_rounding_contains_binary_float_and_inner_rounding_is_required():
    value = rational(.1)
    assert value == Fraction.from_float(.1)
    result = five_point_box([value] * 10, [value] * 10, denominator=10)
    lows = list(map(int, result["outer_squared_lower_numerators"]))
    highs = list(map(int, result["outer_squared_upper_numerators"]))
    assert all(Fraction(a, 10) <= value <= Fraction(b, 10) for a, b in zip(lows, highs))
    # A very coarse outer box can meet the singular boundary although the
    # original point is a strictly positive definite simplex. No fit follows.
    coarse = five_point_box([Fraction(1, 3)] * 10, [Fraction(1, 3)] * 10, denominator=1)
    assert coarse["status"] == "unresolved" and coarse["subset_fit_witness"] is None


def test_fixed_thresholds_and_exact_anchored_matrix():
    assert THRESHOLDS == (Fraction(0), Fraction(3, 100), Fraction(1, 20), Fraction(3, 40), Fraction(1, 10))
    assert gram_integer([1] * 10) == [[2 if i == j else 1 for j in range(4)] for i in range(4)]
    simplex = 1 - np.eye(6, dtype=np.int64)
    result = certify(simplex, numerical_guard=0)
    exclusions = [row["global_3d_fit_excluded"] for row in result["thresholds"]]
    assert exclusions[0]
    assert exclusions == sorted(exclusions, reverse=True)
    with pytest.raises(ValueError, match="symmetric"):
        simplex[0, 1] = 2
        certify(simplex, thresholds=(0,))
