"""Read full and half diamonds independently of stress in conformal geometry controls."""

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq


HERE = Path(__file__).resolve().parent
WEIGHTS = np.array([-9 / 16, 81 / 160, 9 / 160])
QUAD_NODES, QUAD_WEIGHTS = leggauss(5)


def flat_volume(duration):
    return np.pi * duration**4 / 24


def endpoints(s, boost, curvature_parameter):
    # U=t+x, W=t-x in a null-rotated chart with constant transverse position.
    # A^2 dU/ds = boost, dW/ds = 1/boost are the proper geodesic equations.
    def elapsed(u):
        return u + curvature_parameter * u**3 / 3 - boost * s

    u = brentq(elapsed, -1.5, 1.5, xtol=5e-15)
    return u, s / boost


def marked_volume(start, stop, boost, curvature_parameter, mark):
    u0, w0 = endpoints(start, boost, curvature_parameter)
    u1, w1 = endpoints(stop, boost, curvature_parameter)
    u = (u1 + u0) / 2 + (u1 - u0) / 2 * QUAD_NODES
    # Integrate transverse disks and W analytically. The U integrand has degree eight.
    integrand = (
        (u - u0)
        * (u1 - u)
        * (1 + curvature_parameter * u * u) ** 2
        * (1 + mark * u * u)
    )
    return np.pi * (w1 - w0) ** 2 * np.dot(QUAD_WEIGHTS, integrand) / 8


def direction_counts(duration, boost, curvature_parameter, mark):
    pairs = ((-duration / 2, duration / 2), (-duration / 2, 0.0), (0.0, duration / 2))
    counts = np.array(
        [marked_volume(a, b, boost, curvature_parameter, mark) for a, b in pairs]
    ) / flat_volume(duration)
    return np.array([counts[0] - counts[1:].sum(), *counts[1:]])


def control_panel(duration, nx, curvature_parameter, mark, rng):
    boosts = (1.0, (5 + 4 * nx) / 3, (5 - 4 * nx) / 3)
    means = np.stack(
        [direction_counts(duration, v, curvature_parameter, mark) for v in boosts]
    )
    # Disjoint pieces of one Poisson census per direction. Each half has
    # 1/16 the flat volume of the full diamond. Their mean has a factor eight.
    assert np.all(means > 0)
    correction = np.array([1.0, -13 / 3, -13 / 3])
    exposure = int(2**32 * (0.25 / duration) ** 6)
    counts = rng.poisson(means * exposure)
    scale = 36 / duration**2

    def inferred(rows):
        return float(scale * (WEIGHTS @ (rows @ correction)))

    variance = (WEIGHTS**2) @ (means @ correction**2) / exposure
    return dict(
        duration=duration,
        null_x=nx,
        exposure=exposure,
        counts=counts.tolist(),
        exact_center_ricci=-2 * curvature_parameter * (1 + nx) ** 2,
        mean_read=inferred(means),
        poisson_read=inferred(counts / exposure),
        poisson_standard_error=float(scale * np.sqrt(variance)),
        full_only_mean=float(30 / duration**2 * (WEIGHTS @ means.sum(axis=1))),
    )


def controls():
    seed = 190271
    rng = np.random.default_rng(seed)
    return dict(
        seed=seed,
        census="synthetic Poisson marked volume; central intensity one",
        panels=[
            dict(
                conformal_quadratic=c,
                marker_quadratic=mark,
                reads=[
                    control_panel(h, nx, c, mark, rng)
                    for nx in (1.0, -1.0, 0.0, 2**-0.5)
                    for h in (0.25, 0.125, 0.0625, 0.03125)
                ],
            )
            for c, mark in ((0.0, 0.5), (-0.25, 0.0), (-0.25, 0.5), (0.25, 0.5))
        ],
    )


def sharp_support(board, route):
    """Check every off-site block used by the shortest-path amplitude proof."""
    from scipy.sparse.csgraph import connected_components
    from ..common_front.measure import hop

    adjacency = board.adjacency().tocsr()
    assert connected_components(adjacency, directed=False)[0] == 1
    identity = np.eye(2, dtype=int)
    rows = []
    for charge in (1, -1):
        for odd in (False, True):
            blocks = hop(board, route, charge, odd).tobsr(blocksize=(4, 4))
            blocks.sort_indices()
            adjacency.sort_indices()
            np.testing.assert_array_equal(blocks.indices, adjacency.indices)
            np.testing.assert_array_equal(blocks.indptr, adjacency.indptr)
            marked = np.flatnonzero(np.any(blocks.data != np.eye(4), axis=(1, 2)))
            assert len(marked) == 4
            seen = set()
            for i in marked:
                origin = np.searchsorted(blocks.indptr, i, side="right") - 1
                edge = frozenset((int(origin), int(blocks.indices[i])))
                block = blocks.data[i]
                if edge == frozenset(route[:2]):
                    factor = block[::2, ::2]
                    np.testing.assert_array_equal(block, np.kron(factor, identity))
                else:
                    assert edge == frozenset(route[2:])
                    factor = block[:2, :2]
                    np.testing.assert_array_equal(block, np.kron(identity, factor))
                assert np.min(factor) >= 0 and np.max(factor) > 0
                seen.add(edge)
            assert len(seen) == 2
            rows.append(
                dict(
                    charge=charge,
                    odd=odd,
                    directed_blocks=len(blocks.data),
                    nonidentity_blocks=len(marked),
                )
            )
    return rows


def native_relay():
    from scipy.sparse.csgraph import shortest_path
    from ..common_front import certificate as front
    from ..common_front.measure import Ruler, graph, optics

    board = optics.pairing.formation.federation(1)
    route = optics.pairing.formation.region(board)
    ruler = Ruler(board, route, 3)
    low, high, _, scale = front.weights(ruler, front.inputs())
    distances = []
    for costs in (low, high):
        _, _, network = graph(board.n, [route[0]], costs[ruler.index], ruler.edges)
        assert board.n * int(costs.max()) < 2**52
        d = shortest_path(network, directed=False, indices=route[0])
        assert np.all(np.isfinite(d)) and np.all(d == d.astype(np.int64))
        distances.append(d.astype(np.int64))

    def volume(time):
        # Each port contributes the allowed waiting time between the two tips.
        # A smaller travel distance gives the upper bound on the census.
        return [
            sum(max(0, time * scale - 2 * int(a)) for a in d)
            for d in reversed(distances)
        ]

    panels = []
    for time in (2, 4, 8, 16, 32, 64):
        lo, hi = volume(time)
        lo_half, hi_half = volume(time // 2)
        panels.append(
            dict(
                duration=time,
                full_ticks=[lo, hi],
                half_ticks=[lo_half, hi_half],
                ratio=[lo / hi_half, hi / lo_half],
            )
        )
    return dict(
        station=int(route[0]),
        ports=board.n,
        order=ruler.order,
        clock="reciprocal blank-link unit",
        ticks_per_unit=scale,
        census="one marker per port per clock unit",
        panels=panels,
        sharp_support=sharp_support(board, route),
    )


def run():
    return dict(geometric_comparison=controls(), native_rest_relay=native_relay())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(result)
    else:
        print(result, end="")
