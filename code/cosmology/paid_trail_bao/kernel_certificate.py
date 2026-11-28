"""Rational verification of mixtures and dual work bounds on measured annular maps."""

from bisect import bisect_left, bisect_right
from collections import defaultdict
from fractions import Fraction as Q

from . import well_certificate as kernel


def work_factor():
    def atan(x):
        series = sum((-1) ** k * x ** (2 * k + 1) / (2 * k + 1) for k in range(16))
        return series, series + x**33 / 33

    a, b = atan(Q(1, 5))
    e, f = atan(Q(1, 239))
    pi_lo, pi_hi = 16 * a - 4 * f, 16 * b - 4 * e
    x = Q(1, 50)
    log_lo = 2 * sum(x ** (2 * k + 1) / (2 * k + 1) for k in range(8))
    log_hi = log_lo + 2 * x**17 / (17 * (1 - x * x))
    return 4 * pi_lo * log_lo, 4 * pi_hi * log_hi


class AnnularMap:
    def __init__(self, annuli):
        self.edges, weights = ([kernel.decimal(x) for x in a] for a in annuli)
        self.slopes = [
            w / (b * b - a * a)
            for a, b, w in zip(self.edges[:-1], self.edges[1:], weights, strict=True)
        ]
        self.values = [Q(0)]
        for a, b, s in zip(self.edges[:-1], self.edges[1:], self.slopes, strict=True):
            self.values.append(self.values[-1] + s * (b - a))

    def primitive(self, r):
        i = bisect_right(self.edges, r) - 1
        if i < 0:
            return Q(0)
        if i >= len(self.slopes):
            return self.values[-1]
        return self.values[i] + (r - self.edges[i]) * self.slopes[i]

    def density(self, r, h):
        if h:
            return (self.primitive(r + h) - self.primitive(abs(r - h))) / (r * h)
        i = bisect_left(self.edges, r)
        if i < len(self.edges) and self.edges[i] == r:
            return (
                (self.slopes[i - 1] if i else 0)
                + (self.slopes[i] if i < len(self.slopes) else 0)
            ) / r
        return 2 * self.slopes[i - 1] / r if 0 < i < len(self.edges) else Q(0)


def intervals(rows):
    for r, v, e, b in zip(
        *(rows[k] for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")),
        strict=True,
    ):
        r, v, e, b = map(kernel.decimal, (r, v, e, b))
        yield r, max(Q(0), v - 2 * e) ** 2 - b * b, (v + 2 * e) ** 2 - b * b


def witness(rows, annuli, mass, fit):
    table = AnnularMap(annuli)
    radii = [kernel.decimal(x) for x in fit["offset_kpc"]]
    weights = [kernel.decimal(x) for x in fit["weights"]]
    material = kernel.decimal(mass)
    f_lo, f_hi = work_factor()
    margins = []
    for r, lo, hi in intervals(rows):
        minus, plus = Q(49, 50) * r, Q(51, 50) * r
        value = material * sum(
            w * (table.density(minus, h) - table.density(plus, h))
            for h, w in zip(radii, weights, strict=True)
        )
        limits = sorted((value / f_lo, value / f_hi))
        margins.append(min(limits[0] - lo, hi - limits[1]))
    if min(margins) <= 0:
        raise ArithmeticError("Deposited mixture is outside an original speed interval")
    strength = max(sum(w for w in weights if w > 0), -sum(w for w in weights if w < 0))
    return dict(
        verified_rings=len(margins),
        minimum_margin_speed2=float(min(margins)),
        upper=float(strength),
    )


def dual(rows, annuli, mass, coefficients, signed):
    y = [Q(x) for x in coefficients]
    edges, weights = ([kernel.decimal(x) for x in a] for a in annuli)
    material = kernel.decimal(mass)
    events = defaultdict(Q)
    objective = Q(0)
    for wi, (r, lo, hi) in zip(y, intervals(rows), strict=True):
        if not wi:
            continue
        objective += wi * (lo if wi > 0 else hi)
        for factor, rad in ((wi, Q(49, 50) * r), (-wi, Q(51, 50) * r)):
            for knot, slope in kernel.slope_events(
                rad, edges, weights, [(1, Q(1))]
            ).items():
                events[knot] += factor * material * slope
    if objective <= 0:
        raise ArithmeticError("No positive dual lower bound")
    knots = sorted(events)
    value = slope = peak = Q(0)
    trough = Q(0)
    for a, b in zip(knots, knots[1:]):
        slope += events[a]
        if a == 0:
            peak = max(peak, slope)
            trough = min(trough, slope)
        value += slope * (b - a)
        peak = max(peak, value / b)
        trough = min(trough, value / b)
    assert value == 0
    denom = peak - trough if signed else peak
    assert denom > 0
    # The exact lower bound follows from rational geometry and the lower
    # endpoint of the independently enclosed positive work factor.
    lower = work_factor()[0] * objective / denom
    return dict(
        lower=float(lower),
        knots=len(knots),
        positive_rings=sum(x > 0 for x in y),
        negative_rings=sum(x < 0 for x in y),
    )
