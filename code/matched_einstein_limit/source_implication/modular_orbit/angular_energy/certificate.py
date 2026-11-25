"""Enclose a star marginal and its native-energy representation from its preparation."""

import numpy as np
from mpmath import iv

from model.experiments.interval_echo_angular import modular
from model.experiments.interval_echo_einstein_station import cap

from .stationary import dephasing_bound

DEN = 2**75


def filter_bounds():
    # The eight-witness offer has mode odds (1/sqrt(2))^8 and ^9.
    assert cap.standard_count() == 8
    even = iv.mpf(2) ** (-iv.mpf(8) / 32)
    odd = iv.mpf(2) ** (-iv.mpf(9) / 32)
    bounds = []
    for value in ((even + odd) / 2, (even - odd) / 2):
        # Floating-point guesses are certified by interval containment.
        floor = int(float(value.mid) * DEN)
        margin = 2**24
        low, high = floor - margin, floor + margin
        assert 0 < low and iv.mpf(low) / DEN < value < iv.mpf(high) / DEN
        bounds.append((low, high))
    return bounds


def reference_matrix(upper, coefficients):
    def rounding(a):
        return (a + (DEN - 1 if upper else 0)) // DEN

    diagonal, off = [pair[int(upper)] for pair in coefficients]
    word_count = len(cap.WORDS)
    forward = np.eye(word_count, dtype=object) * DEN
    for edge in cap.cell().seams:
        left, right = cap.seam_indices(edge)
        a, b = forward[left], forward[right]
        forward[left] = rounding(diagonal * a + off * b)
        forward[right] = rounding(off * a + diagonal * b)
    state = rounding(forward.T @ forward)
    for _ in range(3):
        state = rounding(state @ state)
    # The local charge subtracts its vacuum value. Positive normalizers of
    # the full state and the rotation average therefore cancel.
    total = np.zeros_like(state)
    rotations = modular.rotations()
    for g in rotations:
        idx = [cap.WORD_INDEX[tuple(sorted(g[i] for i in word))] for word in cap.WORDS]
        total[np.ix_(idx, idx)] += state
    return total


def star_marginal(total):
    star = tuple(sorted((0, *cap.cell().neighbours[0])))
    region = cap.CapRestriction(cap.WORDS, star)
    # The repository's partial trace preserves the integer representation.
    marginal = region.reduce(total)
    ring = list(star[1:])
    x = ring[0]
    near = next(y for y in ring if y in cap.cell().neighbours[x])
    far = next(y for y in ring if x != y and y not in cap.cell().neighbours[x])
    idx = {q: star.index(q) + 1 for q in star}
    return tuple(
        int(marginal[i, j])
        for i, j in (
            (0, 0),
            (1, 1),
            (1, idx[x]),
            (idx[x], idx[x]),
            (idx[x], idx[near]),
            (idx[x], idx[far]),
        )
    )


def star_coefficients(stationary):
    weights = filter_bounds()
    lower, upper = (reference_matrix(side, weights) for side in (False, True))
    low, high = map(star_marginal, (lower, upper))
    vacuum, u, v, z, w, y = [iv.mpf([a, b]) / DEN for a, b in zip(low, high)]
    five = iv.sqrt(5)
    # One pole and the normalized uniform ring give a two-dimensional block.
    e = z + 2 * w + 2 * y
    off = five * v
    gap = iv.sqrt((u - e) ** 2 + 4 * off**2)
    small, large = ((u + e - gap) / 2, (u + e + gap) / 2)
    assert small > 0
    slope = (iv.ln(vacuum / large) - iv.ln(vacuum / small)) / gap
    scalar = iv.ln(vacuum / small) - slope * small
    L00, L11, L01 = scalar + slope * u, scalar + slope * e, slope * off
    # The projections of the constant and axial triplet on the block are
    # (1,sqrt(5))/sqrt(12) and (1,1)/2. A transverse triplet has half
    # its total norm squared in the first Fourier ring doublet.
    a = (L00 + 2 * five * L01 + 5 * L11) / 12
    b = (L00 + (1 + five) * L01 + five * L11) / (2 * iv.sqrt(12))
    axial = (L00 + 2 * L01 + L11) / 4
    ring1 = z + (five - 1) * w / 2 - (five + 1) * y / 2
    assert ring1 > 0
    c = iv.ln(vacuum / ring1) / 2
    d = axial - c
    error = iv.mpf(0)
    trace_error = iv.mpf(0)
    if stationary:
        trace_error = dephasing_bound(lower, upper, DEN)
        ring2 = z - (five + 1) * w / 2 + (five - 1) * y / 2
        floor = min(small.a, ring1.a, ring2.a)
        assert floor > trace_error and vacuum.a > trace_error
        error = trace_error / (floor - trace_error) + trace_error / (
            vacuum.a - trace_error
        )
    perturb = iv.mpf([-1, 1]) * error
    return (
        a + perturb,
        a - c - d / 3 + 2 * perturb,
        b + perturb,
        d + 2 * perturb,
        error,
        trace_error / vacuum,
    )


def solve2(a, b, c, d, y, z):
    determinant = a * d - b * c
    assert determinant > 0 or determinant < 0
    return (y * d - b * z) / determinant, (a * z - y * c) / determinant


def verify(stationary=False):
    iv.dps = 40
    ka, km, kb, kd, log_error, trace_error_ratio = star_coefficients(stationary)
    q = 1 / iv.sqrt(5)

    def energy(x, y):
        a = iv.mpf(4) / 3
        trace = (17 - q) / 4
        axis = (iv.mpf(17) * (x**2 + y**2) / 2 - x * y) / 4
        return a - trace / 3, 8 * (x + y) / iv.sqrt(48), (3 * axis - trace) / 2

    edge = [energy(x, y) for x, y in ((1, q), (q, q), (q, -q), (-q, -q), (-1, -q))]
    (ms, bs, ds), (mr, br, dr) = edge[:2]
    spoke, rim = solve2(5 * bs, 5 * br, 5 * ds, 5 * dr, kb, kd)
    # A linear functional vanishing on each internal spoke and rim.
    dipole_weight, quadrupole_weight = solve2(bs, ds, br, dr, ms, mr)
    signal = km - dipole_weight * kb - quadrupole_weight * kd
    other = [m - dipole_weight * b - quadrupole_weight * d for m, b, d in edge[2:]]
    assert spoke > 0 and rim > 0 and signal > 0 and all(t < 0 for t in other)
    common = signal / (30 * ms)
    scalar = ka - (30 * common + 5 * (spoke + rim)) * iv.mpf(4) / 3
    assert common < 0
    return dict(
        log_change_norm_bound=str(log_error),
        trace_change_over_unnormalized_local_vacuum=str(trace_error_ratio),
        monopole=str(km),
        dipole=str(kb),
        quadrupole=str(kd),
        spoke_weight=str(spoke),
        rim_weight=str(rim),
        rim_excess_over_constant_interior=str(rim - spoke),
        total_energy_weight=str(common),
        identity=str(scalar),
        separating_dipole_weight=str(dipole_weight),
        separating_quadrupole_weight=str(quadrupole_weight),
        modular_witness=str(signal),
        other_seam_witness=[str(x) for x in other],
    )
