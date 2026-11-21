"""Exact five-point obstructions to a three-dimensional distance fit.

Ten squared distances determine a 4 by 4 anchored Gram matrix affinely.
Positive definiteness at all 1024 box corners proves positive definiteness
throughout the box, by convexity. One such five-point obstruction excludes
a three-dimensional fit of all six stations. Failure of this sufficient
test does not establish a fit.

An inner-box positive definite midpoint and an indefinite corner give a
different certificate: the first loss of positive definiteness along their
segment is positive semidefinite of rank at most three. This proves a fit
for that subset within that box, without combining different subsets. For
confidence envelopes it proves compatibility with some values in the
envelope, not a fit of every possible underlying distance matrix.

All decisions use integer determinants. The caller supplies an absolute
numerical enclosure for computed distances; this code does not certify
the floating propagation that produced them.
"""

from fractions import Fraction
from itertools import combinations
from math import isfinite

import numpy as np

def require(condition, message):
    if not condition:
        raise ValueError(message)


THRESHOLDS = (Fraction(0), Fraction(3, 100), Fraction(1, 20), Fraction(3, 40), Fraction(1, 10))
GRID = 10 ** 12
NUMERICAL_GUARD = Fraction(1, 10 ** 10)
PAIRS = tuple(combinations(range(5), 2))


def rational(value):
    """Convert the represented value exactly; missing data remain missing."""
    if value is None:
        return None
    if isinstance(value, (float, np.floating)):
        return Fraction.from_float(float(value)) if isfinite(float(value)) else None
    require(not isinstance(value, (bool, np.bool_)), "boolean distance or bound")
    return Fraction(int(value)) if isinstance(value, np.integer) else Fraction(value)


def floor(value):
    return value.numerator // value.denominator


def ceil(value):
    return -floor(-value)


def determinant(matrix):
    """Exact fraction-free elimination, with row pivoting and divisibility checks."""
    require(all(x == int(x) for row in matrix for x in row), "determinant entries must be integral")
    a = [[int(x) for x in row] for row in matrix]
    n = len(a)
    require(all(len(row) == n for row in a), "determinant requires a square integer matrix")
    if not n:
        return 1
    sign, previous = 1, 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss division is not exact")
                a[i][j] = numerator // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def gram_integer(squares):
    """Return 2Q times the anchored Gram for squared distances stored over Q."""
    require(len(squares) == 10, "five points require ten squared distances")
    require(all(x == int(x) for x in squares), "squared-distance numerators must be integral")
    d = [[0] * 5 for _ in range(5)]
    for (i, j), value in zip(PAIRS, squares):
        d[i][j] = d[j][i] = int(value)
    return [[d[0][i] + d[0][j] - d[i][j] for j in range(1, 5)] for i in range(1, 5)]


def leading_minors(matrix):
    return [determinant([row[:k] for row in matrix[:k]]) for k in range(1, 5)]


def principal_minors(matrix):
    return [determinant([[matrix[i][j] for j in indices] for i in indices])
            for k in range(1, 5) for indices in combinations(range(4), k)]


def _corner(lower, upper, mask):
    return [upper[k] if mask & (1 << k) else lower[k] for k in range(10)]


def five_point_box(lower_squared, upper_squared, denominator=GRID):
    """Certify one rational squared-distance box, with outward and inward grids.

    Input endpoints are exact Fractions. Multiplication, floor and ceiling
    occur in rational arithmetic, avoiding floating inward rounding.
    """
    lower = [rational(x) for x in lower_squared]
    upper = [rational(x) for x in upper_squared]
    require(type(denominator) is int and denominator > 0 and len(lower) == len(upper) == 10,
            "invalid five-point rational box")
    require(all(a is not None and b is not None and 0 <= a <= b for a, b in zip(lower, upper)),
            "invalid squared-distance bounds")
    outer_lo = [floor(x * denominator) for x in lower]
    outer_hi = [ceil(x * denominator) for x in upper]
    minima = [None] * 4
    failed, first_bad, eigen_bound = 0, None, None
    for mask in range(1 << 10):
        gram = gram_integer(_corner(outer_lo, outer_hi, mask))
        minors = leading_minors(gram)
        minima = [x if old is None else min(old, x) for old, x in zip(minima, minors)]
        if min(minors) <= 0:
            failed += 1
            if first_bad is None:
                first_bad = mask
        else:
            # For positive eigenvalues, product of any three <= (trace/3)^3.
            bound = Fraction(27 * minors[-1], sum(gram[i][i] for i in range(4)) ** 3 * 2 * denominator)
            eigen_bound = bound if eigen_bound is None else min(eigen_bound, bound)
    result = dict(
        status="excluded" if failed == 0 else "unresolved",
        excludes_3d_fit=failed == 0, squared_distance_denominator=denominator,
        pair_order=[list(pair) for pair in PAIRS],
        outer_squared_lower_numerators=[str(x) for x in outer_lo],
        outer_squared_upper_numerators=[str(x) for x in outer_hi],
        corners_checked=1 << 10, nonpositive_definite_corners=failed,
        minimum_integral_leading_minors=[str(x) for x in minima],
        minor_denominators=[str((2 * denominator) ** k) for k in range(1, 5)],
        first_nonpositive_corner=first_bad,
        uniform_gram_eigenvalue_lower_bound=str(eigen_bound) if failed == 0 else None,
        subset_fit_witness=None,
    )
    if failed == 0:
        return result

    inner_lo = [ceil(x * denominator) for x in lower]
    inner_hi = [floor(x * denominator) for x in upper]
    if any(a > b for a, b in zip(inner_lo, inner_hi)):
        result["unresolved_reason"] = "no_inner_grid_box; failure of the outer test does not prove a fit"
        return result
    # Midpoint squared distances have denominator 2Q, so this matrix is 4Q G.
    middle = gram_integer([a + b for a, b in zip(inner_lo, inner_hi)])
    center_minors = principal_minors(middle)
    if min(center_minors) >= 0 and center_minors[-1] == 0:
        witness = dict(kind="positive_semidefinite_inner_midpoint", rank_at_most=3,
                       midpoint_integral_principal_minors=[str(x) for x in center_minors])
    elif min(leading_minors(middle)) > 0:
        witness = None
        for mask in range(1 << 10):
            gram = gram_integer(_corner(inner_lo, inner_hi, mask))
            minors = principal_minors(gram)
            if min(minors) >= 0 and minors[-1] == 0:
                witness = dict(kind="positive_semidefinite_inner_corner", corner=mask, rank_at_most=3)
                break
            if min(minors) < 0:
                witness = dict(kind="positive_definite_midpoint_to_negative_direction_corner",
                               corner=mask, negative_principal_minor=str(min(minors)), rank_at_most=3,
                               parameter="first loss of positive definiteness on the segment, in (0,1)")
                break
    else:
        witness = None
    if witness is not None:
        witness.update(inner_squared_lower_numerators=[str(x) for x in inner_lo],
                       inner_squared_upper_numerators=[str(x) for x in inner_hi],
                       squared_distance_denominator=denominator)
        result.update(status="subset_compatible", subset_fit_witness=witness)
    else:
        result["unresolved_reason"] = "corner test does not exclude and no subset fit is certified"
    return result


def certify(distances, lower=None, upper=None, *, numerical_guard=NUMERICAL_GUARD,
            thresholds=THRESHOLDS, denominator=GRID, confidence_failure_probability=None):
    """Check all six five-station subsets at predeclared relative tolerances.

    Relative error means |fitted-true| <= epsilon*true on every pair. Supplied
    simultaneous confidence intervals are expanded by the absolute numerical
    guard before forming the outer fit box. Statistical coverage is inherited
    from those intervals; no independence among distance estimates is needed.
    """
    data = np.asarray(distances, dtype=object)
    require(data.shape == (6, 6), "six-station distance matrix required")
    guarded = rational(numerical_guard)
    require(guarded is not None and guarded >= 0, "invalid absolute numerical guard")
    require(confidence_failure_probability is None or 0 < confidence_failure_probability < 1,
            "invalid confidence failure probability")
    has_intervals = lower is not None or upper is not None
    require(not has_intervals or (lower is not None and upper is not None), "both distance envelopes required")
    lo = np.asarray(lower if has_intervals else data, dtype=object)
    hi = np.asarray(upper if has_intervals else data, dtype=object)
    require(lo.shape == hi.shape == (6, 6), "distance envelope shape differs")
    bounds = {}
    for i, j in combinations(range(6), 2):
        values = [rational(x) for x in (lo[i, j], lo[j, i], hi[i, j], hi[j, i])]
        if any(x is None for x in values):
            bounds[i, j] = None
            continue
        a, ar, b, br = values
        require(a == ar and b == br and 0 <= a <= b, "distance envelopes must be symmetric and nonnegative")
        bounds[i, j] = (max(Fraction(0), a - guarded), b + guarded)
    reports = []
    for threshold in thresholds:
        error = rational(threshold)
        require(error is not None and 0 <= error < 1, "relative tolerance must lie in [0,1)")
        subsets = []
        for stations in combinations(range(6), 5):
            selected = [bounds[stations[i], stations[j]] for i, j in PAIRS]
            if any(value is None for value in selected):
                result = dict(status="unresolved", excludes_3d_fit=False,
                              unresolved_reason="missing_distance_bound", subset_fit_witness=None)
            else:
                squared_low = [((1 - error) * a) ** 2 for a, _ in selected]
                squared_high = [((1 + error) * b) ** 2 for _, b in selected]
                result = five_point_box(squared_low, squared_high, denominator)
            result["stations"] = list(stations)
            subsets.append(result)
        reports.append(dict(relative_tolerance=str(error), relative_tolerance_float=float(error),
                            global_3d_fit_excluded=any(x["excludes_3d_fit"] for x in subsets),
                            excluded_subsets=sum(x["excludes_3d_fit"] for x in subsets),
                            compatible_subsets=sum(x["status"] == "subset_compatible" for x in subsets),
                            unresolved_subsets=sum(x["status"] == "unresolved" for x in subsets),
                            subsets=subsets))
    return dict(method="exact five-point Gram box corners and integral principal minors",
                distance_numerical_guard=str(guarded), rational_grid_denominator=denominator,
                numerical_scope="conditional on the supplied absolute distance enclosure; all determinant decisions are exact",
                data_domain="simultaneous distance envelope" if has_intervals else "enclosed Born-distance values",
                confidence_failure_probability=confidence_failure_probability if has_intervals else None,
                confidence_statement="exclusions hold whenever all supplied distance intervals cover" if has_intervals else None,
                no_global_positive_claim="compatible five-point subsets need not share a six-point realization",
                thresholds=reports)
