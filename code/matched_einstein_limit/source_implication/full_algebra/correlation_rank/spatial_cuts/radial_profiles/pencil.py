"""Uniform rank and hidden-direction bounds for two affine charts of a real pencil."""

from math import ceil, floor, isqrt
from fractions import Fraction

import numpy as np

from .. import measure as cuts

UNIT = 1 << 80


def integer_boxes(array):
    """Centre and radius integers divided by UNIT; all rounding is outward."""
    lo, hi = (
        np.vectorize(
            lambda x: f(cuts.survey.dyadic(x._mpi_[k]) * UNIT), otypes=[object]
        )(array)
        for k, f in ((0, floor), (1, ceil))
    )
    center = (lo + hi) // 2
    return center, hi - center


def integer_proposal(array):
    # Exact rounding of the proposal, independently of the numerical optimizer.
    return np.vectorize(
        lambda x: round(Fraction.from_float(float(x)) * UNIT), otypes=[object]
    )(array)


def frobenius_upper(array):
    return 1 + isqrt(sum(x * x for x in array.flat))


def audit(base, linear, area_scale, maximum_depth=25):
    """Enclose all B(t), -1 <= t <= 1, by subdividing closed dyadic intervals."""
    b, br = integer_boxes(base)
    velocity, vr = integer_boxes(linear)
    count = deepest = 0
    min_signal = None
    largest_defect = Fraction(0)
    todo = [(0, 0)]
    while todo:
        # t in [(n-1)/2^d, (n+1)/2^d].
        n, d = todo.pop()
        den = 1 << d
        cm = den * b + n * velocity
        cr = den * br + (abs(n) + 1) * vr
        candidate = np.linalg.pinv(np.asarray(cm, dtype=float) / (UNIT * den))
        right = integer_proposal(candidate)
        size = cm.shape[0]
        scale = den * UNIT**2
        # Use |L R| for the one common parameter displacement, before taking
        # entrywise absolute values. The other uncertainties are independent.
        dom = (
            abs(np.eye(size, dtype=object) * scale - cm @ right)
            + abs(velocity @ right)
            + cr @ abs(right)
        )
        defect = Fraction(frobenius_upper(dom), scale)
        if defect >= Fraction(1, 2):
            if d >= maximum_depth:
                raise ValueError("could not certify the augmented material matrix")
            todo.extend(((2 * n - 1, d + 1), (2 * n + 1, d + 1)))
            continue
        count += 1
        deepest = max(deepest, d)
        largest_defect = max(largest_defect, defect)
        # W R has a genuine inverse. Applying R(WR)^-1 to the last
        # unit vector gives a hidden tangent with scaled dA = +1.
        norm_r = Fraction(frobenius_upper(right), UNIT)
        norm_last = Fraction(frobenius_upper(right[:, -1]), UNIT)
        upper = norm_last + norm_r * defect / (1 - defect)
        signal = Fraction("23.1042361") / (area_scale * upper)
        min_signal = signal if min_signal is None else min(min_signal, signal)
    return dict(
        certified_boxes=count,
        max_depth=deepest,
        maximum_defect_ppm=ceil(largest_defect * 10**6),
        minimum_entropy_slope_piconats=floor(min_signal * 10**12),
    )
