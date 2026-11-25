"""Physical routes, local pulses and invariant material sectors."""

from functools import cache
from itertools import combinations
from math import atan2, pi

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from ..conjugate_flight.measure import closed_walk, seam_data


@cache
def fixture():
    """Choose two disjoint arms from one owner to two different neighboring owners."""
    b = federation(0)
    bridges = sorted(e for e in b.external if e[0] // 12 == 0)

    def paths(g, n, forbidden):
        def extend(path):
            if len(path) == n:
                yield path
                return
            for p in sorted(b.neighbours[path[-1]]):
                if p // 12 == g // 12 and p not in path and p not in forbidden:
                    yield from extend([*path, p])

        yield from extend([g])

    g, h = bridges[0]
    j, k = bridges[1]
    for left in paths(g, 5, {j}):
        for second in paths(j, 3, set(left)):
            first = left[::-1] + next(paths(h, 3, set()))
            other = second[::-1] + next(paths(k, 3, set()))
            if len(set(first + other)) != 14:
                continue
            active = {tuple(sorted(t)) for r in (first, other) for t in zip(r, r[1:])}
            touched = {e for e in b.seams if set(e).intersection(first + other)}
            # Each phase source can open beside a filled buffer on its own owner's side.
            shutters = {}
            for u in [first[2], first[7], other[0], other[5]]:
                available = [
                    p
                    for p in b.neighbours[u]
                    if p not in first + other
                    and b.degree[p] == 5
                    and p // 12 == u // 12
                    and p not in shutters.values()
                ]
                if available:
                    shutters[u] = min(available)
            if len(shutters) == 4:
                return first, other, sorted(touched - active), shutters
    raise AssertionError("no two-owner receiver layout")


def geometry_record():
    import networkx as nx

    b = federation(0)
    a, d, cuts, posts = fixture()
    g = nx.Graph(b.seams)
    # Park the source inventory off the bench during the writes on unit pairs.
    forbidden = set(a + d) | set(posts.values()) | {p for e in cuts for p in e}
    stores = [p for p in range(b.n) if p not in forbidden][:5]
    w = np.ones(b.n, dtype=int)
    for p, n in zip(stores, [0, 0, 0, 2, 2], strict=True):
        w[p] = n

    def transpose(u, v, graph):
        path = nx.shortest_path(graph, u, v)
        for x, y in closed_walk(path):
            w[x], w[y] = w[y], w[x]

    private = [p for p in a + d if b.degree[p] == 5]
    for p in private:
        transpose(stores[0], p, g)
        assert w[p] == 0
        transpose(stores[0], p, g)
    shutter_edges = {tuple(sorted(e)) for e in posts.items()}
    for e in cuts:
        if e not in shutter_edges:
            assert all(w[p] == 1 for p in e)
            g.remove_edge(*e)
    # Native swaps on the open graph put each defect into its laboratory source.
    for s, p in zip(stores, [a[0], a[3], d[1], a[4], d[2]], strict=True):
        transpose(s, p, g)
    assert all(w[p] == 1 for e in shutter_edges for p in e)
    graph_routes = []
    for r, buffer, sites in (
        (a, a[2], [a[0], a[3]]),
        (a, a[7], [a[5]]),
        (d, d[0], [d[1]]),
        (d, d[5], [d[3]]),
    ):
        for p in sites:
            bench = g.subgraph(r + [posts[buffer]])
            graph_routes.append(nx.shortest_path(bench, p, posts[buffer]))
    local = candidate.native_law
    echo = local.echo()
    for a0 in (0, 1):
        col = np.eye(99)[:, local.SRC.index((a0, a0))]
        np.testing.assert_allclose(echo @ (echo @ col), col, atol=1e-13)
    return dict(
        charged_arm=a,
        second_arm=d,
        receiver_owners=[a[-1] // 12, d[-1] // 12],
        private_zero_writes=private,
        silent_seams=[list(e) for e in cuts],
        shutter_posts=sorted([list(e) for e in posts.items()]),
        initial_stores=stores,
        address_phase_routes=graph_routes,
    )


def sector(n, zeros, doubles, strength, cut=None):
    """Compress the full blank-seam columns after equal-load writes on exterior edges."""
    blank, _ = seam_data(strength)
    # The blank columns are integral halves at the declared contrast settings.
    op = np.rint(blank * 2).astype(int)
    np.testing.assert_allclose(op, blank * 2, atol=1e-12)
    words = []
    for vacancies in combinations(range(n), zeros):
        for positive in combinations(sorted(set(range(n)) - set(vacancies)), doubles):
            w = np.ones(n, dtype=int)
            w[list(vacancies)], w[list(positive)] = 0, 2
            words.append(tuple(map(int, w)))
    index = {w: i for i, w in enumerate(words)}
    h = np.zeros((len(words), len(words)), dtype=int)
    for j, w in enumerate(words):
        for a in range(n - 1):
            if a == cut:
                continue
            b = a + 1
            column = op[:, 3 * w[a] + w[b]]
            for z in np.flatnonzero(column):
                q = list(w)
                q[a], q[b] = divmod(int(z), 3)
                h[index[tuple(q)], j] += column[z]
    return words, h / 2


def rotate_frame(block, strength, loads):
    """One physical partial swap between opposite private phase dwells."""
    # Fix the first eigenvector with positive entries in both order words.
    angle = atan2(-2 * block[0, 1], block[1, 1] - block[0, 0]) / 2
    blank, private = seam_data(strength)
    p, q = loads
    rows = [3 * p + q, 3 * q + p]
    swap = expm(-2j * angle * blank[np.ix_(rows, rows)])
    rate = strength * (private[3 * p] - private[3 * q])

    def phase(theta):
        dwell = (np.sign(rate) * theta % (2 * pi)) / abs(rate)
        return np.diag([np.exp(-1j * strength * private[3 * i] * dwell) for i in loads])

    frame = phase(-pi / 2) @ swap @ phase(pi / 2)
    frame *= np.exp(-1j * np.angle(frame[0, 0]))
    frame = np.real_if_close(frame)
    np.testing.assert_allclose(
        block @ frame, frame @ np.diag(np.diag(frame.T @ block @ frame)), atol=1e-12
    )
    return frame


def local_frames(strength):
    # Relative to filled material; separate one-vacancy and neutral-pair calibrations.
    wc, hc = sector(8, 1, 0, strength)
    wn, hn = sector(6, 1, 1, strength)
    cc = [wc.index(w) for w in [(0, 1, 1, 1, 1, 1, 1, 1), (1, 0, 1, 1, 1, 1, 1, 1)]]
    nn = [wn.index(w) for w in [(1, 0, 2, 1, 1, 1), (1, 2, 0, 1, 1, 1)]]
    blocks = [hc[np.ix_(cc, cc)], hn[np.ix_(nn, nn)]]
    frames = [
        rotate_frame(x, strength, loads)
        for x, loads in zip(blocks, [(0, 1), (0, 2)], strict=True)
    ]
    return blocks, frames


def frame_columns(words, arms, frames):
    """Ordered product-energy channels, with filled unaddressed ports."""
    from itertools import product

    orders = [(0, 1)] + [(0, 2)] if len(arms) == 2 else [(0, 2)]
    slots = []
    for bits in product((0, 1), repeat=len(arms)):
        w = [1] * len(words[0])
        for (left, r), b, order in zip(arms, bits, orders, strict=True):
            w[left], w[r] = order if not b else order[::-1]
        slots.append(words.index(tuple(w)))
    local = frames[0] if len(arms) == 1 else np.kron(*frames)
    return np.eye(len(words))[:, slots] @ local, slots


def read_local(y, words, settings):
    """Each observer rotates and decodes only its own orbitals after its buffer reads one."""
    out = y.copy()
    labels = np.full((len(words), len(settings)), 2, dtype=int)
    index = {w: i for i, w in enumerate(words)}
    for j, (left, r, buffer, order, frame) in enumerate(settings):
        for i, w in enumerate(words):
            if w[buffer] != 1 or (w[left], w[r]) != order:
                continue
            z = list(w)
            z[left], z[r] = z[r], z[left]
            pair = [i, index[tuple(z)]]
            labels[pair, j] = [0, 1]
            out[pair] = frame.conj().T @ out[pair]
    p = np.zeros((3,) * len(settings))
    np.add.at(p, tuple(labels.T), abs(out) ** 2)
    return p
