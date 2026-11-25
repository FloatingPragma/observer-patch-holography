"""Read every compact factor of the same 143-defect parent against a physical store."""

from collections import defaultdict
from math import asin, log, pi, sqrt
import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path

from ...row_04_compact import material_reference as store

WAITS = 2.0 ** np.arange(-6, 13)
SHOTS = 2**18
CLOCK_SHOTS = 2**28


def read(bank):
    """Unwrap X,Y one doubling at a time, independent of the line identity."""
    xy = bank[:, :, :2]
    z = (xy[:, :, 0] - xy[:, :, 1]) / bank.sum(axis=2)
    wrapped = np.arctan2(z[:, 1], z[:, 0])
    phase = wrapped[0]
    for angle in wrapped[1:]:
        phase = angle + 2 * pi * round((2 * phase - angle) / (2 * pi))
    return phase / WAITS[-1]


def clock(rng):
    law = store.LAW
    ids = [law.SRC.index(x) for x in ((0, 1), (1, 0))]
    dimer = store.anchor.free_seam(store.K)[np.ix_(ids, ids)]
    result = []
    for wait in WAITS:
        xy = []
        for port, extra in ((0, 0), (1, pi / 2)):
            bright = abs(expm(-1j * (wait + extra) * dimer)[port, 0]) ** 2
            b = int(rng.binomial(CLOCK_SHOTS, np.clip(round(bright, 13), 0, 1)))
            xy.append([b, CLOCK_SHOTS - b, 0])
        result.append(xy)
    return np.asarray(result)


def phase_radius(starts, banks):
    eps = sqrt(2 * log(2 * banks / 0.01) / starts)
    # Include the numerical probability-grid resolution of the simulated banks.
    return asin(sqrt(2) * (eps + 2e-9))


def measure(board, charts, records, store_port, native_home, wave):
    rng = np.random.default_rng(350671)
    store_owner = store_port // 12
    assert board.degree[store_port] == 6 and all(
        r["owner"] != store_owner for r in records
    )
    pred = shortest_path(
        board.adjacency(),
        indices=store_port,
        directed=False,
        unweighted=True,
        return_predecessors=True,
    )[1]
    pooled = defaultdict(lambda: np.zeros((len(WAITS), 2, 3), dtype=np.int64))
    sizes = defaultdict(int)
    max_route = 0
    for r in records:
        owner, sign, tag = r["owner"], r["signed_load"], r["compact_branch"]
        gate = r["sharp_port"] % 12
        # A chain of full swaps, followed by the single endpoint splitter
        # and then the inverse swaps, works through arbitrary spectator words.
        route = [r["sharp_port"]]
        while route[-1] != store_port:
            route.append(int(pred[route[-1]]))
        assert all(b in board.neighbours[a] for a, b in zip(route, route[1:]))
        max_route = max(max_route, len(route) - 1)
        x = wave(charts[owner], sign, tag, gate)[:, gate]
        h = native_home(charts[owner], sign)
        for neighbour in board.neighbours[store_port]:
            if neighbour // 12 == owner:
                h[neighbour % 12, neighbour % 12] -= store.K + 0.5
                assert abs(x[neighbour % 12]) < 1e-12
        energy, frame = np.linalg.eigh(h)
        amp = np.exp(-1j * WAITS[:, None] * energy) @ abs(frame.conj().T @ x) ** 2
        split, combiner = store.optics(sign)
        kind = f"{sign},{tag}"
        sizes[kind] += 1
        for j, a in enumerate(amp):
            for q, phase in enumerate((1, -1j)):
                # A routed private phase on the store arm sets the second quadrature.
                arms = combiner @ np.array([a * split[0], phase * split[1]])
                probs = np.maximum(
                    np.r_[abs(arms) ** 2, abs(split[0]) ** 2 * (1 - abs(a) ** 2)], 0
                )
                np.testing.assert_allclose(probs.sum(), 1, atol=3e-11)
                pooled[kind][j, q] += rng.multinomial(
                    SHOTS, np.round(probs / probs.sum(), 12)
                )
    c = clock(rng)
    unit = read(c)
    numbanks = (len(pooled) + 1) * 2 * len(WAITS)
    dc = phase_radius(CLOCK_SHOTS, numbanks) / WAITS[-1]
    observed = 0
    error = 0
    rows = {}
    for kind, counts in sorted(pooled.items()):
        n = sizes[kind]
        phase = phase_radius(n * SHOTS, numbanks)
        assert 3 * phase < pi and 69 * WAITS[0] + phase < pi
        f = read(counts)
        observed += n * f
        error += n * phase / WAITS[-1]
        rows[kind] = dict(
            owners=n, read_angular_frequency=f, pooled_counts=counts.tolist()
        )
    assert unit - dc > 0
    interval = [(observed - error) / (unit + dc), (observed + error) / (unit - dc)]
    return dict(
        store_port=store_port,
        maximum_routed_seam_distance=max_route,
        waits=WAITS.tolist(),
        shots_per_owner=SHOTS,
        clock_shots=CLOCK_SHOTS,
        owner_groups=rows,
        dimer_counts=c.tolist(),
        dimer_frequency=unit,
        formation_in_dimer_units=observed / unit,
        simultaneous_interval_99=interval,
    )
