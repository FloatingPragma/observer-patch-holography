"""Planar annular convolution of the finite tide and its continuous-scale bounds."""

from fractions import Fraction as Q
from functools import lru_cache
import math

import numpy as np

ANGLE = 2**22
WEIGHT = 2**24
SUB = 8
LENGTH = 2**15


def shape_data(exterior, signal, donor):
    knots = list(map(int, exterior["shells"]))
    assert knots == list(range(1, max(knots) + 1))

    def profile(name):
        return [
            -Q(
                row["comparisons"][signal]["distant_equal_range"]["contrasts"][name][
                    "projected"
                ]["mean"]
            )
            for row in exterior["shells"].values()
        ]

    values = profile(signal)
    ordinary = [x**3 * v for x, v in zip(knots, profile(donor), strict=True)]
    sign = 1 if signal == "source" else -1
    assert all(sign * v >= 0 for v in values) and values[0] == values[-1] == 0
    return knots, values, {k: v for k, v in zip(knots, ordinary, strict=True) if v > 0}


def integral(knots, values, power, lo, hi):
    """Exact moment of a piecewise-linear compact kernel in native hop units."""
    result = Q(0)
    for x, y, f, g in zip(knots[:-1], knots[1:], values[:-1], values[1:], strict=True):
        a, b = max(Q(x), lo), min(Q(y), hi)
        if a >= b:
            continue
        slope = (g - f) / Q(y - x)
        result += (f - x * slope) * (b ** (power + 1) - a ** (power + 1)) / (power + 1)
        result += slope * (b ** (power + 2) - a ** (power + 2)) / (power + 2)
    return result


def kernel_bins(knots, values):
    mass = []
    bins = []
    for i in range(knots[0] * SUB, knots[-1] * SUB):
        weight = integral(knots, values, 2, Q(i, SUB), Q(i + 1, SUB))
        if weight:
            bins.append(i)
            mass.append(abs(weight))
    moment = sum(mass)
    shares = [m * WEIGHT / moment for m in mass]
    return (
        np.array(bins),
        np.array(list(map(math.floor, shares))),
        np.array(list(map(math.ceil, shares))),
        moment,
    )


def arc_angle(radius, boundary, separation):
    """Angular first moment of a disk at separation s from an off-center observer."""
    s = np.asarray(separation, float)
    if boundary == 0:
        return np.zeros_like(s)
    if radius == boundary:
        return np.sqrt(np.maximum(0, 4 - (s / radius) ** 2))
    low, high = abs(radius - boundary), radius + boundary
    safe = np.where((s > low) & (s < high), s, high)
    return np.sqrt((high - safe) * (high + safe) * (safe - low) * (safe + low)) / (
        radius * safe
    )


@lru_cache(maxsize=4)
def quadrature(order):
    return np.polynomial.legendre.leggauss(order)


def annular_force(radius, hop, knots, values, order=12):
    """Mean inward acceleration per unit mass, before the ordinary conversion.

    Integrate each disk edge separately. A sin² substitution removes the
    square-root endpoint slopes at both tangencies. The returned columns
    are uniform-area annuli. Lengths are in kpc.
    """
    leg, weights = quadrature(order)
    phase = np.pi * (leg + 1) / 2
    node = (1 - np.cos(phase)) / 2
    weights = weights * np.pi * np.sin(phase) / 4
    density = 1 / (np.pi * np.diff(np.r_[0, radius**2]))
    knots = np.asarray(knots, float)
    values = np.asarray(values, float)
    out = np.zeros((len(radius), len(radius)))
    for x0, x1, f, g in zip(
        knots[:-1], knots[1:], values[:-1], values[1:], strict=True
    ):
        if f == g == 0:
            continue
        x = radius[:, None]
        R = radius[None, :]
        left = np.maximum(abs(x - R) / hop, x0)
        right = np.maximum(left, np.minimum((x + R) / hop, x1))
        width = right - left
        for u, w in zip(node, weights, strict=True):
            s = left + width * u
            rho = hop * s
            upper = (x + R) ** 2 - rho**2
            lower = rho**2 - (x - R) ** 2
            angle = np.sqrt(np.maximum(0, upper) * np.maximum(0, lower)) / (x * rho)
            out += w * width * s * s * (f + (g - f) * (s - x0) / (x1 - x0)) * angle
    return np.diff(out, axis=1, prepend=np.zeros((len(radius), 1))) * density


def angle_floor(x, R, s):
    """All lengths integer; bound the square root with an integer square root."""
    if s == 0:
        return 2 * ANGLE if R == x else 0
    if not abs(x - R) < s < x + R:
        return 0
    numer = ((x + R) ** 2 - s * s) * (s * s - (x - R) ** 2) * ANGLE**2
    return math.isqrt(numer // ((x * s) ** 2))


def angular_boxes(x, edges):
    """Exact trajectory-independent extrema enclosures for each arc."""
    end = x + int(edges[-1])
    step = max(1, end // 360)
    grid = np.unique(np.r_[np.arange(0, end + 1, step), abs(x - edges), x + edges, end])
    left, right = grid[:-1], grid[1:]
    low = np.zeros((len(left), len(edges)), np.int64)
    high = np.zeros_like(low)
    for j, R in enumerate(map(int, edges)):
        y = np.array([angle_floor(x, R, int(s)) for s in grid])
        low[:, j] = np.minimum(y[:-1], y[1:])
        high[:, j] = np.minimum(2 * ANGLE, np.maximum(y[:-1], y[1:]) + 1)
        peak = abs(x * x - R * R)
        on = (left**2 <= peak) & (peak <= right**2)
        high[on, j] = (2 * ANGLE * min(x, R) + x - 1) // x
        outside = (right <= abs(x - R)) | (left >= x + R)
        high[outside, j] = 0
    inner_lo = np.pad(low[:, :-1], ((0, 0), (1, 0)))
    inner_hi = np.pad(high[:, :-1], ((0, 0), (1, 0)))
    return grid, low - inner_hi, high - inner_lo


def density_ticks(radius, stocks):
    radii = [Q(str(x)) for x in radius]
    areas = [b * b - a * a for a, b in zip([Q(0)] + radii[:-1], radii, strict=True)]
    density = []
    for s in stocks:
        cum = [Q(str(v)) for v in s]
        density.append(
            [(b - a) / A for a, b, A in zip([Q(0)] + cum[:-1], cum, areas, strict=True)]
        )
    unit = Q(2) ** (math.ceil(math.log2(float(max(max(row) for row in density)))) - 24)
    low = np.array([[math.floor(v / unit) for v in row] for row in density], np.int64)
    high = np.array([[math.ceil(v / unit) for v in row] for row in density], np.int64)
    assert np.all(low >= 0) and high.max() <= 2**24
    return low, high, unit


def length_edges(maxr):
    # In kpc; the first interval starts at zero. Outside the final interval
    # no kernel point can reach the disk. All finite endpoints are dyadic.
    stop = math.ceil(math.log2(2 * maxr))
    return np.array(
        [0]
        + [2 ** (k + 12) * (8 + j) for k in range(-12, stop) for j in range(8)]
        + [2 ** (stop + 15)],
        np.int64,
    )


def max_ranges(v, a, b):
    """Max per row of v on all closed index ranges [a,b]."""
    length = b - a + 1
    spans = np.floor(np.log2(length)).astype(int)
    result = np.empty((len(v), *a.shape), np.int64)
    table = v
    for k in range(int(spans.max()) + 1):
        use = spans == k
        if np.any(use):
            result[:, use] = np.maximum(
                table[:, a[use]], table[:, b[use] - (1 << k) + 1]
            )
        jump = 1 << k
        if table.shape[1] > jump:
            table = np.maximum(table[:, :-jump], table[:, jump:])
    return result


def disk_boxes(x, edges, dlo, dhi, scale, ticks, kernel, *, keep_annuli=False):
    bins, wlo, whi, _ = kernel
    grid, alo, ahi = angular_boxes(x, edges)

    def upper(weights):
        return np.maximum(weights, 0) @ dhi.T + np.minimum(weights, 0) @ dlo.T

    # pi times surface density in dyadic density units. Rounding and
    # sign are outward before the subsequent weighted integration.
    # Lower is floor(-upper(-alo)/ANGLE); upper is ceil(upper(ahi)/ANGLE).
    lower = (-upper(-alo)) // ANGLE
    higher = -((-upper(ahi)) // ANGLE)
    # Append the zero-response region beyond the largest possible separation.
    lower = np.pad(lower.T, ((0, 0), (0, 1)))
    higher = np.pad(higher.T, ((0, 0), (0, 1)))
    den = SUB * LENGTH
    left = (scale[:-1, None] * bins * ticks) // den
    right = -((-(scale[1:, None] * (bins + 1) * ticks)) // den)
    a = np.clip(np.searchsorted(grid, left, side="left") - 1, 0, len(grid) - 1)
    b = np.clip(np.searchsorted(grid, right, side="right") - 1, 0, len(grid) - 1)

    def integrated(weights):
        return np.maximum(weights, 0) @ whi + np.minimum(weights, 0) @ wlo

    upper_grid = integrated(max_ranges(higher, a, b))
    lower_grid = -integrated(max_ranges(-lower, a, b))

    # Each annulus can choose its own distribution of lengths. Integrate
    # all native radii at the same length before taking the maximum.
    # Round back to ANGLE units before multiplying by density to keep
    # the interval arithmetic inside int64. Include the zero-field tail.
    def per_annulus(bounds):
        angles = np.pad(bounds.T, ((0, 0), (0, 1)))
        sums = integrated(max_ranges(angles, a, b))
        per_length = -((-sums) // WEIGHT)
        return np.pad(
            np.minimum(np.maximum(0, bounds.max(0))[:, None], per_length),
            ((0, 0), (0, 1)),
        )

    ann_lo, ann_hi = -per_annulus(-alo), per_annulus(ahi)
    if keep_annuli:
        return ann_lo, ann_hi
    independent_hi = -((-upper(ann_hi.max(1))) // ANGLE) * WEIGHT
    independent_lo = (-upper((-ann_lo).max(1))) // ANGLE * WEIGHT
    return lower_grid, upper_grid, independent_lo, independent_hi
