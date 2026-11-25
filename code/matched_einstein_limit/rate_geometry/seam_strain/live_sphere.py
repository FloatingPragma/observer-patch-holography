"""Fit a round three-sphere to five clock stations and predict held ranges.

Rows are ordered (0, 1, 38, 39, 5, 7, 40, 42). This is a comparison
geometry with minimizing geodesic distances, evaluated at population means.
"""

from fractions import Fraction
from itertools import combinations
from math import ceil, floor, isfinite, log10, pi, sqrt

import numpy as np
from scipy.optimize import brentq


def _matrix(values, *, distances):
    x = np.asarray(values, dtype=float)
    if (x.shape != (8, 8) or not np.isfinite(x).all()
            or not np.allclose(x, x.T, atol=1e-13, rtol=0)
            or np.any(np.diag(x) != 0)):
        raise ValueError("Expected eight symmetric station readings with zero diagonal")
    if distances and np.any(x[np.triu_indices(8, 1)] <= 0):
        raise ValueError("Station distances must be positive")
    return (x + x.T) / 2


def all_radii_obstruction(distances):
    """Bound the distance alteration needed for any round S³, or flat R³.

    Six points on S³ have chord distances in R⁴. For geodesic extrema l,u
    the chord ratio is at least l/u by concavity of sin on [0, pi/2].
    If squared chords lie in [L²,U²], centering their matrix yields five
    Gram eigenvalues on 1-perp at least (6L²-4U²)/4: compare with equal
    squared distances (L²+U²)/2 and bound each error row sum. Five positive
    directions contradict affine dimension four. Test the resulting necessary
    ratio sqrt(2/3) in exact arithmetic on the supplied binary64 distances.
    """
    d = _matrix(distances, distances=True)
    best = None
    ratio = sqrt(2 / 3)
    for indices in combinations(range(8), 6):
        edges = [d[i, j] for i, j in combinations(indices, 2)]
        low, high = min(edges), max(edges)
        bound = (low - ratio * high) / (1 + ratio)
        if best is None or bound > best[0]:
            best = bound, indices, low, high
    bound, indices, low, high = best
    low_q, high_q = Fraction.from_float(low), Fraction.from_float(high)
    excluded = 3 * low_q**2 > 2 * high_q**2
    units = max(0, floor(bound * 10**10))
    lower = Fraction(units, 10**10)
    if excluded and 3 * (low_q - lower)**2 <= 2 * (high_q + lower)**2:
        lower -= Fraction(1, 10**10)
    return dict(excluded=excluded, witness_indices=list(indices),
                minimum_distance=float(low), maximum_distance=float(high),
                maximum_distance_change_lower=float(lower),
                exact_lower_rational=str(lower))


def _curvature(d):
    # Schur complement of the first unit vector in the spherical Gram.
    # Dividing by z² removes the trivial rank collapse at zero curvature.
    length = d.max()
    scaled = d / length

    def smallest(angle):
        v = -0.5 * scaled**2 * np.sinc(angle * scaled / (2 * pi))**2
        row = v[0, 1:]
        schur = v[1:, 1:] - row[:, None] - row - angle**2 * np.outer(row, row)
        return np.linalg.eigvalsh(schur)[0]

    # Find a simple PSD rank-four boundary on the principal geodesic branch.
    # This search reports one calibration; it makes no uniqueness claim.
    nodes = np.linspace(0, pi, 129)
    old = smallest(0)
    for left, right in zip(nodes[:-1], nodes[1:], strict=True):
        new = smallest(right)
        if (old < 0 < new) or (new <= 0 < old):
            angle = brentq(smallest, left, right, xtol=2e-14)
            if angle > 1e-8:
                return (angle / length)**2
        old = new
    raise ValueError("No simple positive spherical calibration found on the grid")


def sphere_fit(distances):
    """Report a five-station calibration and its independent held residuals."""
    d = _matrix(distances, distances=True)
    curvature = _curvature(d[:5, :5])
    z = sqrt(curvature)
    gram = np.cos(z * d)
    spectrum = np.linalg.eigvalsh(gram[:5, :5])
    anchors = gram[:4, :4]
    if abs(spectrum[0]) > 1e-10 or spectrum[1] < 1e-8 or np.linalg.eigvalsh(anchors)[0] < 1e-8:
        raise ValueError("Five calibration points need a simple null and four independent anchors")
    coefficients = np.linalg.solve(anchors, gram[:4, :])
    predicted = gram[:, :4] @ coefficients
    norm_squared = predicted.diagonal()
    if np.min(norm_squared) < 1e-14:
        raise ValueError("A held station has no directional prediction from these anchors")
    norms = np.sqrt(norm_squared)
    angles = np.arccos(np.clip(predicted / np.outer(norms, norms), -1, 1))
    comparisons = [
        dict(indices=[i, j], measured_distance=float(d[i, j]),
             normalized_prediction=float(angles[i, j] / z),
             gram_residual=float(predicted[i, j] - gram[i, j]))
        for i, j in combinations(range(4, 8), 2)
    ]
    return dict(curvature=curvature, inverse_radius=z, calibration_indices=list(range(5)),
                anchor_indices=list(range(4)), calibration_gram_eigenvalues=spectrum.tolist(),
                all_distances_on_shortest_branch=bool(z * d.max() <= pi),
                held_norm_squared=norm_squared[5:].tolist(),
                maximum_held_diagonal_residual=float(max(abs(norm_squared[5:] - 1))),
                held_comparisons=comparisons,
                maximum_held_distance_error=float(max(abs(row["measured_distance"] - row["normalized_prediction"])
                                                        for row in comparisons)),
                all_radii=all_radii_obstruction(d))


def sphere_curvature_slope(distances, change, curvature):
    """Differentiate the five-point PSD boundary, with a fixed set of anchors."""
    d = _matrix(distances, distances=True)[:5, :5]
    delta = _matrix(change, distances=False)[:5, :5]
    z = sqrt(curvature)
    values, basis = np.linalg.eigh(np.cos(z * d))
    if not (z > 0 and abs(values[0]) < 1e-9 and values[1] > 1e-8):
        raise ValueError("Expected a positive calibration with a simple Gram null")
    v = basis[:, 0]
    weights = np.outer(v, v) * np.sin(z * d)
    denominator = np.sum(weights * d)
    if abs(denominator) < 1e-10 * d.max():
        raise ValueError("The spherical calibration is tangent in radius")
    return float(-2 * curvature * np.sum(weights * delta) / denominator)


def sphere_fixed_volume(curvature, flat_radius):
    """Area and dA/dK at fixed V=4π ell³/3 on a supplied round three-sphere."""
    from mpmath import mp

    if not (isfinite(curvature) and isfinite(flat_radius) and curvature > 0 and flat_radius > 0):
        raise ValueError("Expected positive curvature and flat reference radius")
    # Both the enclosed volume and its area derivative have cancellations
    # at a small aperture. Preserve their precision in the flat limit.
    log_angle = log10(curvature)/2 + log10(flat_radius)
    with mp.workdps(50 + max(0, ceil(-4 * log_angle))):
        k, ell = mp.mpf(curvature), mp.mpf(flat_radius)
        flat_angle = mp.sqrt(k) * ell
        target = 2 * flat_angle**3 / 3
        if target >= mp.pi:
            raise ValueError("Requested volume exceeds the available spherical ball")
        lo, hi = mp.mpf(0), 2*flat_angle if flat_angle < mp.mpf("0.25") else mp.pi
        for _ in range(160):
            x = (lo + hi) / 2
            if x - mp.sin(x) * mp.cos(x) < target:
                lo = x
            else:
                hi = x
        x = (lo + hi) / 2
        area = 4 * mp.pi * mp.sin(x)**2 / k
        volume_factor = x - mp.sin(x) * mp.cos(x)
        slope = 4 * mp.pi / k**2 * (mp.mpf(3)/2 * volume_factor * mp.cot(x) - mp.sin(x)**2)
        return dict(proper_radius=float(x / mp.sqrt(k)), area=float(area),
                    area_per_curvature_at_fixed_volume=float(slope))
