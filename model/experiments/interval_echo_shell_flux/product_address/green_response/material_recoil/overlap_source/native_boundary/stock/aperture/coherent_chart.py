"""Echo-defined marker masks on a stored coherent one-defect chart."""

from itertools import product
import numpy as np
from scipy.sparse.csgraph import dijkstra

from . import occupied_release as flight
from . import phase_receipts as phase
from . import thermal_radar as radar

HERE = flight.HERE
VOLUMES = (165, 192, 219, 240)


def geometry(board, edges, clamp, center=None):
    live, hops, central = radar.chart(board, edges)
    if center is None:
        center = central
    x, y = np.where(np.triu(live.toarray(), 1))
    # Keep flags in the actual source order, even if sorting the recorded edges would permute them.
    near = [
        (a, b) if (hops[center, a], a) < (hops[center, b], b) else (b, a)
        for a, b in edges
    ]
    survey = {}
    for c in (0, 1):
        g = live.copy()
        g.data[:] = 1 + c * (g.indices == clamp)
        direct = (65 * dijkstra(g, indices=center, directed=True)).astype(int)
        assisted = [
            direct[a] + 67 + (65 * dijkstra(g, indices=b, directed=True)).astype(int)
            for a, b in near
        ]
        for f in range(4):
            arrival = direct.copy()
            for k, route in enumerate(assisted):
                if f & (2 >> k):
                    arrival = np.minimum(arrival, route)
            for reverse in (False, True):
                order = np.argsort(arrival[::-1] if reverse else arrival, kind="stable")
                if reverse:
                    order = board.n - 1 - order
                for volume in VOLUMES:
                    mask = np.zeros(board.n, bool)
                    mask[order[:volume]] = True
                    assert sum(mask) == volume
                    cut = mask[x] ^ mask[y]
                    survey[(volume, reverse, c, f)] = (mask, int(sum(cut)))
    return center, live, survey


def ordinary_site(board, edges, center, first):
    # Select the clamp with the largest minimum ordinary-only change over
    # the pre-enrolled volumes and tie orders. Zero flags = no written shortcut.
    _, hops, _ = radar.chart(board, edges)
    candidates = []
    for x in range(board.n):
        if x in (center, first, *np.ravel(edges)) or hops[center, x] > 8:
            continue
        _, _, survey = geometry(board, edges, x, center)
        if not all(
            mask[x] and mask[np.ravel(edges)].all() for mask, _ in survey.values()
        ):
            continue
        changes = [
            abs(survey[n, b, 1, 0][1] - survey[n, b, 0, 0][1])
            for n in VOLUMES[:-1]
            for b in (False, True)
        ]
        candidates.append(
            (sum(v != 0 for v in changes), min(changes), sum(changes), -x)
        )
    best = max(candidates)
    return -best[-1], dict(
        eligible=len(candidates),
        nonzero_panels=best[0],
        minimum_links=best[1],
        sum_absolute_links=best[2],
    )


def strata(native):
    p, q = flight.sector_prior(native)
    records = []
    probs = []
    for c, h, s, f in product(range(2), range(4), range(9), range(4)):
        a, b = flight.SECTORS[s]
        clicks = np.sin(phase.dwell.DWELL * np.array([a, b]) / 2) ** 2
        v = flight.CONDITIONAL[h, s] * np.prod(
            np.where([f // 2, f % 2], clicks, 1 - clicks)
        )
        if v == 0:
            continue
        records.append((c, h, s, f))
        probs.append([0.5 * p[h] * v, (0.75 if c else 0.25) * p[h] * v, 0.5 * q[h] * v])
    weights = np.array(probs).T
    np.testing.assert_allclose(weights.sum(1), 1, atol=1e-14, rtol=0)
    return np.array(records, int), weights
