"""Physical ring, ordered wedge control and moving private-reference exposures."""

from functools import lru_cache
from itertools import combinations
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

RING = (0, 5, 4, 9, 8, 6, 2, 11)
PAIR = ((0, 1, 4), (4, 5, 1), (1, 2, 5), (5, 6, 2))
ONE = tuple((j, j + 1, None) for j in range(4))
LAW = candidate.native_law


def geometry():
    board = federation(0)
    private = [i for i in range(12) if board.degree[i] == 5]

    def first(prefix):
        if len(prefix) == 8:
            return tuple(prefix) if prefix[0] in board.neighbours[prefix[-1]] else None
        for p in private:
            if p not in prefix and (not prefix or p in board.neighbours[prefix[-1]]):
                found = first([*prefix, p])
                if found:
                    return found
        return None

    if first([]) != RING:
        raise ArithmeticError("Guide differs from the enrolled private-port cycle")
    # Address read can use the seven blank links of this path. Starting
    # from differences the maximum is the occupied level because at most two holes exist.
    for count in (1, 2):
        for sites in combinations(range(8), count):
            loads = np.ones(8, dtype=int)
            loads[list(sites)] = 0
            rel = np.r_[0, np.cumsum(np.diff(loads))]
            if not np.array_equal(rel - rel.max() + 1, loads):
                raise ArithmeticError("Terminal differences do not recover the holes")
    return dict(
        ports=list(RING), reader_edges=[list(x) for x in zip(RING[:-1], RING[1:])]
    )


@lru_cache(None)
def generator(count, k, law):
    """Unordered sharp words. Only a tagged crossing in an ordered wedge gains its parity."""
    words = list(combinations(range(8), count))
    idx = {w: i for i, w in enumerate(words)}
    holes = np.array([[int(p in w) for p in range(8)] for w in words])
    h = np.zeros((len(words), len(words)))
    edge = LAW.seam_term(k)
    for col, word in enumerate(words):
        for a in range(8):
            b = (a + 1) % 8
            pair = (int(a not in word), int(b not in word))
            loc = LAW.SRC.index(pair)
            h[col, col] += edge[loc, loc]
            if pair[0] == pair[1]:
                continue
            old, new = (a, b) if a in word else (b, a)
            moved = [new if v == old else v for v in word]
            parity = sum(u > v for j, u in enumerate(moved) for v in moved[j + 1 :])
            sign = (-1) ** parity if law == "exchange_minus" else 1
            if law == "flux_pi" and a == 7:
                sign = -sign
            h[idx[tuple(sorted(moved))], col] += (
                sign * edge[LAW.SRC.index(pair[::-1]), loc]
            )
    return words, holes, h


def source(count):
    return (0,) if count == 1 else (0, 4)


def target(count):
    return (4,) if count == 1 else (2, 6)


def initial_stocks():
    import networkx as nx

    board = federation(0)
    graph = nx.Graph(board.seams)
    used = set(RING)
    paths = []
    for origin in (RING[0], RING[4]):
        candidates = nx.single_source_shortest_path(graph, origin)
        path = min(
            (v for i, v in candidates.items() if i not in used),
            key=lambda v: (len(v), v),
        )
        paths.append(path)
        used.add(path[-1])
    return paths


def trace(count, kind, duration, h=0.25):
    k = 0 if kind == "zero_contrast" else 8
    words, holes, H = generator(count, k, kind)
    pulse = expm(-1j * h / 8 * H)
    if kind == "parked":
        pulse = np.diag(np.exp(-1j * h / 8 * H.diagonal()))
    src = words.index(source(count))
    wave = np.zeros((len(words), 2), dtype=complex)
    wave[src, :] = 1
    steps = round(duration / h)
    if abs(duration / h - steps) > 1e-12 or steps % 2:
        raise ValueError("Need an even number of complete appointments")
    snapshots = [(0, 0, np.abs(wave.T) ** 2)]
    s = (np.arange(steps) + 0.5) / steps
    for leg, (a, b, park) in enumerate(ONE if count == 1 else PAIR):
        start = np.column_stack([holes[:, a], holes[:, -a % 8]])
        delta = np.column_stack([holes[:, b], holes[:, -b % 8]]) - start
        if park is not None:
            start += np.column_stack([holes[:, park], holes[:, -park % 8]])
        # Relative to the unit host, positive private contrast k n² contributes -k per vacancy.
        bank = np.exp((0.5j * k * h) * (start[None] + s[:, None, None] * delta[None]))
        for half, chunk in enumerate(np.array_split(bank, 2), 1):
            for phase in chunk:
                wave = phase * (pulse @ (phase * wave))
            snapshots.append((leg + 1, half, np.abs(wave.T) ** 2))
    norms = np.sum(abs(wave) ** 2, axis=0)
    if np.max(abs(norms - 1)) > 5e-9:
        raise ArithmeticError("Norm drift")
    # Remove only numerical norm error; both arms exhaust the Hilbert space.
    wave /= np.sqrt(norms)
    exits = np.array([abs(wave[:, 0] + sign * wave[:, 1]) ** 2 / 4 for sign in (1, -1)])
    incoherent = np.repeat(np.mean(abs(wave.T) ** 2, axis=0)[None] / 2, 2, axis=0)
    return words, snapshots, exits, incoherent
