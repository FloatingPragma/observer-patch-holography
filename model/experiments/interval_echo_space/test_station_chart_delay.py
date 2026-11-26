from fractions import Fraction
from itertools import combinations
from math import isqrt

from .station_chart_delay import (
    COMMON_DELAY, MILLI_COORDINATES, archived_degrees, archived_ranges,
    delay_polynomial, ranges_within_tolerance,
)


def test_published_ranges_have_no_exact_common_delay_three_space_fit():
    assert archived_degrees() == {20: 0, 80: 0, 320: 0}


def test_euclidean_distances_with_a_known_delay_share_the_physical_root():
    # Orthogonal axes of an Euler brick give 3D stations with integer pair ranges.
    sites = [(0, 0, 0), (44, 0, 0), (0, 117, 0), (0, 0, 240),
             (-44, 0, 0), (0, -117, 0), (0, 0, -240)]
    ranges = {}
    for i, j in combinations(range(7), 2):
        squared = sum((x-y)**2 for x, y in zip(sites[i], sites[j], strict=True))
        distance = isqrt(squared)
        assert distance**2 == squared
        ranges[i, j] = distance + 3
    polynomial = delay_polynomial(ranges)
    assert polynomial.eval(3) == 0 and polynomial.eval(2) != 0


def test_shared_delay_and_explicit_points_fit_every_archived_population_dyad():
    # Enclose each maximum relative error to 0.001 percentage point, without
    # a square root, optimization or any binary floating point arithmetic.
    bounds = {20: (4974, 4975), 80: (4978, 4979), 320: (4424, 4425)}
    for cells, ranges in archived_ranges():
        assert len(ranges) == 21
        points = MILLI_COORDINATES[cells]
        lo, hi = (Fraction(edge, 100000) for edge in bounds[cells])
        assert ranges_within_tolerance(ranges, points, COMMON_DELAY, hi)
        assert not ranges_within_tolerance(ranges, points, COMMON_DELAY, lo)
