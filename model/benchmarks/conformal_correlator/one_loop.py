"""Local half-BPS four-point function in a supplied four-dimensional N=4 theory."""

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import spence

# The first fixture has four equal pairwise separations. All points avoid the inversion pole.
REGULAR = np.eye(4)
SKEW = np.array([[1, 0, 0, 0], [1, 2, 0, 0], [0, 1, 1, 0], [2, 0, 1, 1]], dtype=float)


def distances(points):
    points = np.asarray(points, dtype=float)
    if points.shape != (4, 4) or not np.isfinite(points).all():
        raise ValueError("Four finite Euclidean four-vectors are required")
    differences = points[:, None] - points[None, :]
    squared = np.sum(differences**2, axis=-1)
    if np.min(squared[np.triu_indices(4, 1)]) <= 0:
        raise ValueError("Local insertions must be separated")
    return squared


def box(points, order):
    """Euclidean integral divided by pi squared, in projective Feynman parameters."""
    if type(order) is not int or order < 1:
        raise ValueError("Quadrature order must be a positive integer")
    d = distances(points)
    nodes, weights = leggauss(order)
    u = (1 + nodes[:, None]) / 2
    v = (1 + nodes[None, :]) / 2
    weights2 = weights[:, None] * weights[None, :] / 4
    value = 0.0
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            k, l = [m for m in range(4) if m != i and m != j]
            # Set alpha_i=1 and (alpha_j, alpha_k, alpha_l)=r(1,u,v).
            # The radial Jacobian cancels r^-2; integrate 1/(A+rB)^2 exactly.
            a = d[i, j] + u * d[i, k] + v * d[i, l]
            b = u * d[j, k] + v * d[j, l] + u * v * d[k, l]
            value += np.sum(weights2 / (a * (a + b)))
    return float(value)


def cross_ratio_box(points):
    """Independent dilogarithm for fixtures with nonreal cross-ratio z."""
    d = distances(points)
    scale = d[0, 2] * d[1, 3]
    u, v = d[0, 1] * d[2, 3] / scale, d[0, 3] * d[1, 2] / scale
    real = (1 + u - v) / 2
    if u <= real**2:
        raise ValueError(
            "The complex cross-ratio check requires a nondegenerate fixture"
        )
    z = complex(real, np.sqrt(u - real**2))
    phi = 2 * (spence(1 - z).imag + np.log(abs(z)) * np.angle(1 - z)) / z.imag
    return float(phi / scale)


def correlator(points, coupling, order=24, pairs=((0, 1), (2, 3))):
    """Return the protected two-point factors and the order-a connected correction for SU(2)."""
    if not np.isfinite(coupling) or coupling < 0:
        raise ValueError("The real squared coupling must be finite and nonnegative")
    if (
        len(pairs) != 2
        or any(len(pair) != 2 for pair in pairs)
        or sorted(i for pair in pairs for i in pair) != list(range(4))
    ):
        raise ValueError("The polarizations must form two disjoint pairs")
    d = distances(points)
    d12, d34 = [d[i, j] for i, j in pairs]
    # c_I is the protected two-point normalization. Cross-pair polarizations are orthogonal.
    color_dimension = 3
    c_i = color_dimension / (2 * (4 * np.pi**2) ** 2)
    two_point = [c_i / d12**2, c_i / d34**2]
    disconnected = np.prod(two_point)
    coefficient = -2 * d12 * d34 * box(points, order) / color_dimension
    return dict(
        two_point=two_point,
        disconnected=disconnected,
        normalized_connected_per_unit_coupling=coefficient,
        connected_one_loop=coupling * disconnected * coefficient,
        normalized_four_point_through_one_loop=1 + coupling * coefficient,
    )


def run():
    rows = []
    coupling = 0.01
    for name, points in (("regular", REGULAR), ("skew", SKEW)):
        base = correlator(points, coupling)
        inverted = points / np.sum(points**2, axis=1)[:, None]
        perm = [2, 0, 3, 1]
        moved_pairs = ((perm.index(0), perm.index(1)), (perm.index(2), perm.index(3)))
        correction = base["connected_one_loop"]
        rows.append(
            dict(
                fixture=name,
                points=points.tolist(),
                free=correlator(points, 0.0),
                interacting=base,
                box=box(points, 24),
                relative_quadrature_change=abs(box(points, 12) / box(points, 24) - 1),
                relative_dilog_error=abs(cross_ratio_box(points) / box(points, 24) - 1),
                relative_permutation_error=abs(
                    correlator(points[perm], coupling, pairs=moved_pairs)[
                        "connected_one_loop"
                    ]
                    / correction
                    - 1
                ),
                relative_dilation_error=abs(
                    correlator(2 * points, coupling)["connected_one_loop"]
                    * 2**8
                    / correction
                    - 1
                ),
                relative_inversion_error=abs(
                    correlator(inverted, coupling)["connected_one_loop"]
                    / (correction * np.prod(np.sum(points**2, axis=1) ** 2))
                    - 1
                ),
            )
        )
    return dict(
        theory="Euclidean SU(2) N=4 supersymmetric Yang-Mills",
        coupling_a=coupling,
        perturbative_order=1,
        fixtures=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
