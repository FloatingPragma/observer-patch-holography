"""Capped two-leg radar on the retained joint strain and both original heads.

The fresh mobile is localized at the sending port. The strain and both heads
remain live during both legs, including the interval following a relay hit.
Axes on the retained source are (strain, left coin, right coin). Every miss
or unfinished return receives the same deadline range as a return at the cap.
"""

from functools import lru_cache
from itertools import combinations, permutations
from operator import index

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .source import measure as source

TICKS = 12
STEP = 0.5
WIDTH = 8


def station_rows(board, stations):
    """Packed (station, strain, left coin, right coin) addresses in the full mobile."""
    return (
        np.asarray(stations)[:, None, None] * 4
        + np.arange(2)[None, :, None] * (4 * board.n)
        + np.arange(4)[None, None, :]
    ).reshape(-1, WIDTH)


def station_kernels(board, route, charge, coupled, stations):
    """Unmonitored station-to-station amplitudes at each half-tick."""
    integer_h = source.system(board, route, charge, coupled=coupled)[0]
    generator = (-1j / 40) * integer_h
    rows = station_rows(board, stations)
    size = rows.size
    launch = np.zeros((generator.shape[0], size), dtype=complex)
    launch[rows.ravel()] = np.eye(size)
    sampled = expm_multiply(
        generator,
        launch,
        start=0,
        stop=TICKS * STEP,
        num=TICKS + 1,
        traceA=generator.diagonal().sum(),
    )
    n = len(stations)
    return (
        sampled[:, rows.ravel(), :]
        .reshape(TICKS + 1, n, WIDTH, n, WIDTH)
        .transpose(0, 1, 3, 2, 4)
    )


def first_arrivals(kernels):
    """First-hit Kraus maps indexed by (sender, target, tick, output, input)."""
    steps, n, _, width, _ = kernels.shape
    hits = np.zeros((n, n, steps, width, width), complex)
    for sender, target in permutations(range(n), 2):
        for t in range(1, steps):
            hits[sender, target, t] = kernels[t, target, sender]
            if t > 1:
                prior = (
                    kernels[t - 1 : 0 : -1, target, target] @ hits[sender, target, 1:t]
                )
                hits[sender, target, t] -= prior.sum(axis=0)
    return hits


def gram(amplitude):
    return amplitude.conj().swapaxes(-1, -2) @ amplitude


def route_events(kernels):
    """Inclusive effects for first relay and first completed return at each tick."""
    hits = first_arrivals(kernels)
    n, _, steps, _, _ = hits.shape
    outward = gram(hits)
    returns = np.zeros_like(outward)
    for p, q in permutations(range(n), 2):
        for t in range(2, steps):
            # The measured relay ticks are distinct branches. The full source
            # matrix from the first leg enters the return Kraus map unchanged.
            histories = hits[q, p, t - 1 : 0 : -1] @ hits[p, q, 1:t]
            returns[p, q, t] = gram(histories).sum(axis=0)
    outward.setflags(write=False)
    returns.setflags(write=False)
    return outward, returns


@lru_cache(maxsize=16)
def _directed_events(level, charge, coupled, stations):
    board = source.formation.federation(level)
    route = source.formation.region(board)
    if len(stations) < 2 or len(set(stations)) != len(stations):
        raise ValueError("enroll at least two distinct ports")
    if min(stations) < 0 or max(stations) >= board.n:
        raise ValueError("a radar station lies outside the finite carrier")
    if charge not in (-1, 1):
        raise ValueError("the fresh mobile charge must be −1 or +1")
    return route_events(station_kernels(board, route, charge, coupled, stations))


def directed_events(level=1, charge=-1, coupled=True, stations=None):
    """Return (first relay, first return) effects for ordered stations; unused diagonal is zero."""
    if stations is None:
        board = source.formation.federation(level)
        stations = source.formation.region(board)
    return _directed_events(level, charge, coupled, tuple(index(p) for p in stations))


def range_effects(level=1, charge=-1, coupled=True, stations=None, deadline=TICKS):
    """Orientation-averaged effects for each lexicographic unordered pair of station indices.

    With k half-tick appointments the maximum range is k/4 in native exchange
    time units. Both timeout outcomes are assigned that value.
    """
    k = index(deadline)
    if not 1 <= k <= TICKS:
        raise ValueError("the deadline must be an integer in 1..12")
    _, returns = directed_events(level, charge, coupled, stations)
    rebate = (k - np.arange(k)) * STEP / 2
    capped = k * STEP * np.eye(WIDTH) / 2 - np.einsum(
        "t,pqtij->pqij", rebate, returns[:, :, :k]
    )
    return np.stack(
        [
            (capped[p, q] + capped[q, p]) / 2
            for p, q in combinations(range(len(capped)), 2)
        ]
    )


@lru_cache(maxsize=1)
def _chebyshev_weights():
    from ...source_implication.isometric_collar import certify

    return certify.weights()


def _fixed_kernels(board, route, charge, coupled, stations):
    """Fixed-point Chebyshev station blocks of exp(-it A/40), with A integer.

    Reuse the published 32t Taylor/Chebyshev weights. Here the contraction is
    A/1280, including the strain and both head contacts. Two limbs implement
    an exact floor division by 1280 or 640 after every sparse product.
    """
    from fractions import Fraction
    from math import isqrt

    from ...source_implication.isometric_collar import certify as fixed

    a = source.system(board, route, charge, coupled=coupled)[0]
    row = int(abs(a).sum(axis=1).max())
    if row > 1280 or (a != a.T).nnz:
        raise ValueError(
            "the native generator exceeds the certified Hermitian contraction"
        )
    degree, unit, base, frac = fixed.DEGREE, fixed.UNIT, fixed.BASE, fixed.FRAC
    indices = station_rows(board, stations).ravel()
    width = len(indices)
    co = _chebyshev_weights()
    high = np.zeros((a.shape[0], width), dtype=np.int64)
    high[indices] = (unit // base) * np.eye(width, dtype=np.int64)
    low = np.zeros_like(high)
    previous_high = previous_low = None
    out = np.zeros((2, TICKS, width, width), dtype=object)
    for k in range(degree + 1):
        if row * (int(abs(high).max()) + 1) >= 2**62 or row * base >= 2**62:
            raise OverflowError(
                "fixed-point native sparse product exceeds its checked limb bound"
            )
        center = high[indices].astype(object) * base + low[indices]
        out[k % 2] += (
            (1 if k % 4 in (0, 3) else -1) * co[:, k, None, None] * center[None, :, :]
        )
        if k == degree:
            break
        divisor = 1280 if k == 0 else 640
        upper = a @ high
        lower = a @ low
        upper += lower >> fixed.LIMB
        # The remainder numerator is nonnegative and <1280*2^50<2^61.
        new_low = ((upper % divisor) * base + (lower & fixed.MASK)) // divisor
        new_high = upper // divisor
        if k:
            new_low -= previous_low
            new_high += (new_low >> fixed.LIMB) - previous_high
            new_low &= fixed.MASK
        previous_high, previous_low, high, low = high, low, new_high, new_low
    # Each floor incurs a column error <=sqrt(dim)/unit. The Chebyshev second
    # kind polynomials propagate a perturbation across j steps with norm <=j+1.
    # Rounded weights have absolute value <=3; each exact column has norm <=1.
    column = Fraction(
        3 * (1 + isqrt(a.shape[0])) * degree * (degree + 1) * (degree + 2), 6 * unit
    )
    column += Fraction(degree + 1, frac) + Fraction(1, 10**59)
    # Eight real and eight imaginary entries of a column are floored again.
    # The operator norm of an 8-by-8 error is at most sqrt(8) times the maximum
    # column error. Integer upper factors here cover both square roots.
    error = WIDTH * (column + Fraction(2 * WIDTH, frac))
    n = len(stations)
    data = (
        (out // unit).reshape(2, TICKS, n, WIDTH, n, WIDTH).transpose(0, 1, 2, 4, 3, 5)
    )
    return data, error


@lru_cache(maxsize=8)
def _fixed_events(level, charge, coupled, stations):
    from fractions import Fraction

    from ...source_implication.isometric_collar import certify as fixed

    board = source.formation.federation(level)
    if (
        charge not in (-1, 1)
        or len(stations) < 2
        or len(set(stations)) != len(stations)
        or min(stations) < 0
        or max(stations) >= board.n
    ):
        raise ValueError("invalid charge or list of enrolled stations")
    raw, error = _fixed_kernels(
        board, source.formation.region(board), charge, coupled, stations
    )
    n = len(stations)
    hits = np.zeros((n, n, TICKS + 1, 2, WIDTH, WIDTH), dtype=object)
    returns = np.zeros_like(hits)
    for p, q in permutations(range(n), 2):
        for t in range(1, TICKS + 1):
            hit = raw[:, t - 1, q, p].copy()
            for j in range(1, t):
                hit -= fixed.multiply(raw[:, t - j - 1, q, q], hits[p, q, j])
            hits[p, q, t] = hit
    for p, q in permutations(range(n), 2):
        for t in range(2, TICKS + 1):
            for j in range(1, t):
                branch = fixed.multiply(hits[q, p, t - j], hits[p, q, j])
                returns[p, q, t] += fixed.multiply(fixed.adjoint(branch), branch)
    # Each exact restricted propagator or first-hit Kraus map is a contraction.
    # A floored product has complex operator error <=sqrt(2)*WIDTH/frac.
    rounding = Fraction(2 * WIDTH, fixed.FRAC)
    first = [Fraction(0)]
    events = [Fraction(0)]
    for t in range(1, TICKS + 1):
        first.append(
            error + sum(error + (1 + error) * first[j] + rounding for j in range(1, t))
        )
        total = Fraction(0)
        for j in range(1, t):
            branch = first[j] + first[t - j] + first[j] * first[t - j] + rounding
            total += 2 * branch + branch**2 + rounding
        events.append(total)
    return returns, events


def certified_ranges(level=1, charge=-1, coupled=True, stations=None, deadline=TICKS):
    """Return (integer centers, denominator, absolute operator-norm radius).

    Centers have shape (unordered edge, real/imaginary, 8, 8) and are exactly
    Hermitian. The radius controls each true expectation on every joint density
    of trace one without first converting any integer to floating point.
    """
    from fractions import Fraction

    from ...source_implication.isometric_collar import certify as fixed

    k = index(deadline)
    if not 1 <= k <= TICKS:
        raise ValueError("the deadline must be an integer in 1..12")
    if stations is None:
        stations = source.formation.region(source.formation.federation(level))
    stations = tuple(index(p) for p in stations)
    returns, bounds = _fixed_events(level, charge, coupled, stations)
    centers = []
    for p, q in combinations(range(len(stations)), 2):
        value = np.zeros((2, WIDTH, WIDTH), dtype=object)
        value[0] = 2 * k * fixed.FRAC * np.eye(WIDTH, dtype=object)
        for t in range(2, k):
            value -= (k - t) * (returns[p, q, t] + returns[q, p, t])
        # The orientation average above has denominator 8F. Orthogonal
        # Hermitian symmetrization doubles the denominator and cannot increase
        # the operator error about the true Hermitian result.
        centers.append(value + fixed.adjoint(value))
    radius = sum(
        (Fraction(k - t, 4) * bounds[t] for t in range(2, k)), start=Fraction(0)
    )
    return np.stack(centers), 16 * fixed.FRAC, radius
