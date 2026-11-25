"""Finite, rotation-covariant proper-volume and null screen-area decoders."""

import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erf, erfi

from matched_einstein_limit.rate_geometry.static_directions import design

A, B = 0.5, 0.25
ALPHA, BETA, FIRST = 5 / 24, 3 / 8, 3 / 4


def frame():
    directions, poles, _, _ = design()
    return directions[poles]


def tensor(directions, samples):
    """Inverse of the six icosahedral quadratic evaluations."""
    return 1.25 * np.einsum(
        "i,ij,ik->jk", samples, directions, directions
    ) - 0.25 * sum(samples) * np.eye(3)


def spatial(directions, center, pairs, radius):
    reference = 4 * np.pi * (B * radius) ** 3 / 3
    reads = -6 * (pairs.mean(axis=1) - center) / (reference * (A * radius) ** 2)
    return tensor(directions, reads)


def null(directions, area, affine):
    # Each opening-angle pi/3 cap subtends pi steradians.
    reads = -6 * (area / (np.pi * affine**2) - 1) / affine**2
    smoothed = tensor(directions, reads.mean(axis=1))
    even = (smoothed - ALPHA * np.trace(smoothed) * np.eye(3)) / BETA
    odd = directions.T @ (reads[:, 0] - reads[:, 1]) / (8 * FIRST)
    return even, odd


def solid_directions(axis, lower, order):
    """Unit directions and normalized solid-angle quadrature on a spherical cap."""
    t, weights = leggauss(order)
    z = (1 + lower + (1 - lower) * t) / 2
    phi = 2 * np.pi * np.arange(3 * order) / (3 * order)
    helper = np.eye(3)[np.argmin(abs(axis))]
    e = np.cross(axis, helper)
    e /= np.linalg.norm(e)
    f = np.cross(axis, e)
    vectors = z[:, None, None] * axis + np.sqrt(1 - z**2)[:, None, None] * (
        np.cos(phi)[:, None] * e + np.sin(phi)[:, None] * f
    )
    return vectors.reshape(-1, 3), np.repeat(weights / (6 * order), 3 * order)


def spatial_volumes(directions, shape, radius, order=8):
    # Euclidean normal coordinates for g=I+[x]_cross^T B [x]_cross.
    # Independently integrate sqrt(det(g)) over each actual tagged region.
    n, w = solid_directions(np.array([0.0, 0.0, 1.0]), -1, order)
    nodes, q = leggauss(order)
    t = (nodes + 1) / 2
    grid = (t[:, None, None] * n).reshape(-1, 3)
    weights = (1.5 * q[:, None] * t[:, None] ** 2 * w).ravel()
    volume = 4 * np.pi * (B * radius) ** 3 / 3

    def read(offset):
        x = radius * (B * grid + offset)
        cross = np.cross(x[:, None, :], np.eye(3))
        h = np.eye(3) + cross @ shape @ np.swapaxes(cross, 1, 2)
        return volume * np.dot(weights, np.sqrt(np.linalg.det(h)))

    center = read(0.0)
    pairs = np.array([[read(sign * A * n) for sign in (1, -1)] for n in directions])
    return center, pairs


def screen_areas(directions, even, odd, affine, order=8):
    # In exp(2 phi) Minkowski space with phi(u,un) = -q u²/4, the
    # affine parameter is integral exp(-q u²/2) du, the area factor u²exp(-q u²/2).
    def read(axis):
        n, weights = solid_directions(axis, 0.5, order)
        q = np.einsum("ij,jk,ik->i", n, even, n) + 2 * n @ odd
        u = np.full_like(q, affine)
        nonzero = q != 0
        scale = np.sqrt(abs(q[nonzero]) / 2)
        for _ in range(6):
            z = u[nonzero] * scale
            integral = (
                np.where(q[nonzero] > 0, erf(z), erfi(z)) * np.sqrt(np.pi) / (2 * scale)
            )
            u[nonzero] -= (integral - affine) * np.exp(q[nonzero] * u[nonzero] ** 2 / 2)
        return np.pi * np.dot(weights, u**2 * np.exp(-q * u**2 / 2))

    return np.array([[read(sign * n) for sign in (1, -1)] for n in directions])


def report():
    n = frame()
    shape = np.array([[0.6, 0.2, -0.1], [0.2, 0.5, 0.12], [-0.1, 0.12, 0.8]])
    expected = 3 * (shape - np.trace(shape) * np.eye(3))
    even = np.array([[0.8, -0.3, 0.15], [-0.3, -0.6, 0.24], [0.15, 0.24, 0.5]])
    odd = np.array([0.25, -0.12, 0.4])
    rows = []
    for r in (0.5, 0.25, 0.125, 0.0625):
        center, pairs = spatial_volumes(n, shape, r)
        metric = spatial(n, center, pairs, r)
        areas = screen_areas(n, even, odd, r)
        even_out, odd_out = null(n, areas, r)
        rows.append(
            dict(
                radius=r,
                center_volume=center,
                paired_volumes=pairs.tolist(),
                paired_screen_areas=areas.tolist(),
                spatial_error=float(np.linalg.norm(metric - expected)),
                null_even_error=float(np.linalg.norm(even_out - even)),
                null_odd_error=float(np.linalg.norm(odd_out - odd)),
            )
        )
    return dict(
        axes=n.tolist(),
        spatial_target=expected.tolist(),
        null_even_target=even.tolist(),
        null_odd_target=odd.tolist(),
        prescribed_smooth_comparisons=rows,
    )


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
