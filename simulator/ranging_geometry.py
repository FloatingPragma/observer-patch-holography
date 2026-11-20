"""Six-station geometry from two local first-click records per trial.

The stations compare the outward and return records after the experiment.
On a confirmed round trip the return time excludes the recorded latch dwell.
The time indices start at one. Index zero on either leg records a timeout.
"""

from fractions import Fraction
from itertools import combinations, product
from math import isqrt

import numpy as np

from .metric_certificate import THRESHOLDS, certify


N_STATIONS = 6
HORIZON = 12
RANGE_TRIALS = 2**34
LATTICE = 2**48
RANGE_TV = Fraction((HORIZON + 1)**2, LATTICE)
PAIRS = tuple(combinations(range(N_STATIONS), 2))


def _sampling_radius(n):
    """Round sqrt(12/n) upwards to a rational with denominator 2^40."""
    grid = 2**40
    return Fraction(isqrt(12 * grid**2 // n) + 1, grid)


def distance_bands(records, *, trials=RANGE_TRIALS, tv=RANGE_TV):
    """Simultaneous bands for the symmetrized conditional mean half-return.

    Given k completions, their durations are iid from the conditional
    acquisition distribution. A second Hoeffding bound lower-bounds the
    acquired success probability q. Total variation d between source and
    acquisition changes a conditional mean on [1,H] by at most
    (H-1)*min(1,2*d/q). When positivity cannot be established the pair
    contributes no source range. Counts on one direction carry no weight
    in the mean of the reverse direction.
    """
    counts = np.asarray(records)
    if counts.shape != (2, 6, 6, 13, 13) or counts.dtype != np.int64:
        raise ValueError("two six-station integer first-click histograms required")
    if type(trials) is not int or not 0 < trials <= 2**40 or not isinstance(tv, Fraction) or not 0 <= tv < 1:
        raise ValueError("invalid trial count or acquisition error")
    if np.any(counts < 0) or np.any(counts > trials):
        raise ValueError("invalid counts")
    totals = counts.sum(axis=(-2, -1))
    if not np.all(totals == trials * (1 - np.eye(6, dtype=np.int64))):
        raise ValueError("missing trials or populated self-pair")
    # 60 directed menus: 60 two-sided mean and 60 one-sided success bounds.
    # Their total failure probability is <= 180*exp(-24) < 1e-8.
    menus = []
    for block in counts:
        directed = {}
        for i, j in product(range(6), repeat=2):
            if i == j:
                continue
            joint = block[i, j]
            k = int(joint[1:, 1:].sum())
            band = None
            if k:
                success_low = Fraction(k, trials) - _sampling_radius(trials)
                if tv == 0 or success_low > tv:
                    bias = min(1, 2*tv/success_low) if tv else 0
                    numerator = sum((u + v) * int(joint[u, v]) for u, v in
                                    product(range(1, HORIZON+1), repeat=2))
                    mean = Fraction(numerator, 2*k)
                    radius = (HORIZON-1) * (_sampling_radius(k) + bias)
                    band = (max(1, mean-radius), min(HORIZON, mean+radius))
            directed[i,j] = band
        symmetric = {}
        for i,j in PAIRS:
            forward, reverse = directed[i,j], directed[j,i]
            symmetric[i,j] = (tuple((a+b)/2 for a,b in zip(forward, reverse))
                              if forward is not None and reverse is not None else None)
        menus.append(symmetric)
    return menus


def _squared_bounds(band, tolerance):
    return ((1-tolerance)*band[0])**2, ((1+tolerance)*band[1])**2


def _euclidean_obstruction(bands, tolerance):
    """A positive lower bound on x^T D^2 x for any zero-sum x excludes every Euclidean dimension."""
    for x in product((-1, 0, 1), repeat=6):
        if sum(x) or next((v for v in x if v), 0) != 1:
            continue
        lower = Fraction(0)
        for pair in PAIRS:
            i, j = pair
            c = 2*x[i]*x[j]
            if not c:
                continue
            if bands[pair] is None:
                break
            lo, hi = _squared_bounds(bands[pair], tolerance)
            lower += c * (lo if c > 0 else hi)
        else:
            if lower > 0:
                return {"zero_sum_coefficients": list(x), "quadratic_form_lower": str(lower)}
    return None


def _octahedron(bands, tolerance, opposite_pairs):
    """One shared three-dimensional placement that fits every range in its band."""
    lower, upper = Fraction(0), None
    for pair, band in bands.items():
        if band is None:
            return None
        factor = 4 if pair in opposite_pairs else 2
        left, right = band
        lower = max(lower, ((1-tolerance)*right)**2 / factor)
        cap = ((1+tolerance)*left)**2 / factor
        upper = cap if upper is None else min(upper, cap)
    if lower <= upper and upper > 0:
        return {"radius_squared": str((lower+upper)/2), "opposite_station_pairs": list(opposite_pairs)}
    return None


def classify(bands, opposite_pairs):
    """All three tests consume the same simultaneous distance bounds."""
    if set(bands) != set(PAIRS):
        raise ValueError("all fifteen pair entries required")
    if (len(opposite_pairs) != 3 or any(pair not in PAIRS for pair in opposite_pairs)
            or sorted(i for pair in opposite_pairs for i in pair) != list(range(6))):
        raise ValueError("three disjoint, ordered station pairs required")
    lower = np.full((6,6), None, dtype=object)
    upper = lower.copy()
    for (i,j), band in bands.items():
        if band is None:
            continue
        if len(band) != 2 or not all(isinstance(v, (Fraction,int)) and not isinstance(v,bool) for v in band):
            raise ValueError("rational bounds required")
        lo, hi = map(Fraction, band)
        if not 0 < lo <= hi:
            raise ValueError("invalid range interval")
        lower[i,j] = lower[j,i] = lo
        upper[i,j] = upper[j,i] = hi
    result = certify(lower, lower, upper, numerical_guard=0, confidence_failure_probability=1e-8)
    for tolerance, row in zip(THRESHOLDS, result["thresholds"]):
        any_dim = _euclidean_obstruction(bands, tolerance)
        fit = _octahedron(bands, tolerance, opposite_pairs)
        excluded = row["global_3d_fit_excluded"] or any_dim is not None
        if excluded and fit is not None:
            raise ArithmeticError("contradictory Euclidean certificates")
        row.update(any_dimension_obstruction=any_dim, octahedron=fit,
                   six_station_verdict="excluded" if excluded else "fit" if fit else "unresolved")
    return result
