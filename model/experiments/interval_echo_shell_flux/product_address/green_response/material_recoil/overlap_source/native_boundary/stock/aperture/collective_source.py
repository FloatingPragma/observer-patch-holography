"""Stationary communication on several native four-owner overlap preparations."""

from fractions import Fraction as Q
from functools import lru_cache
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, shortest_path

from . import tidal_release as tidal

PARTIAL = (2, 4, 7, 9)


@lru_cache(maxsize=1)
def records():
    _, words, cells, *_ = tidal.lab.energy.pair.setup()
    original = np.full(len(words), Q(1, len(words)), object)
    repaired = tidal.lab.energy.least_change(cells) / len(words)
    # Only the last external contact can carry a bit on the chosen ruler.
    assert not np.isin(words[:, 4:6], (0, *PARTIAL)).any()
    n = (words[:, :4] == 1).sum(1)
    assert np.all(words[:, :4].sum(1) == 4)
    blank = words[:, -1] == 0
    partial = np.isin(words[:, -1], PARTIAL)
    assert all(n[partial] >= 2)
    out = {}
    for name, weights in (("original", original), ("source", repaired)):
        out[name] = tuple(weights @ f for f in (blank, partial, n, partial * n))
    return words, {"original": original, "source": repaired}, out


def chances(volume, arms, singles=0, *, moments=None):
    """Single global pool, independent returned records before mixing.
    Each partial head sees its own initial population through a coincidence.
    """
    stats = records()[2] if moments is None else moments
    means = [stats[k][2] for k in arms]
    all_singles = sum(means) + singles
    return [
        a + b / 4 + (c + b * (all_singles - n)) / (4 * volume)
        for name in arms
        for a, b, n, c in [stats[name]]
    ]


def ordinary_coefficients():
    # Match the repaired success rate to the original at fixed transmitting
    # heads plus uniformly pooled ordinary singles. Applies when m>0.
    _, _, stats = records()
    a, b, n, c = stats["original"]
    A, B, N, C = stats["source"]
    return ((4 * (A - a) + B - b) / b, (B * N - b * n) / b, (C - B * N - c + b * n) / b)


def check_gate():
    # All four transmitting labels communicate on the occupied support.
    # Blank transpositions on the connected graph permute every load.
    # Symmetric half-pulse steps preserve the uniform stationary sector.
    gate = tidal.fastest_relay().reshape(2, 11, 3, 3)
    for present in product(range(3), repeat=2):
        n1, n2 = present
        pop = (240 - n1 - n2, n1, n2)
        groups = tidal.lab.stationary_heads(pop)
        for m in range(11):
            if m not in (0, *PARTIAL):
                assert not gate[:, groups[m]].any()
        assert np.all(gate[:, 0] == 32)
        if n1 >= 2:
            for m in PARTIAL:
                assert set(groups[m]) == set(PARTIAL)
            # Uniform stationary weight over the four heads. Sum tickets
            # first, retaining the receiver's target in each orientation.
            by_target = np.array([1, 2, 1]) * 32
            sums = gate[:, list(PARTIAL)].sum(1)
            assert np.all(sums[0] == by_target[None, :])
            assert np.all(sums[1] == by_target[:, None])


def fractional_donor(rng, fraction, size):
    if fraction.denominator <= np.iinfo(np.int64).max:
        return rng.integers(fraction.denominator, size=size) < fraction.numerator
    # Compare a uniform binary expansion with the exact rational in base 2^32.
    # Only a tied digit requires more random bits.
    chosen = np.zeros(size, bool)
    pending = np.arange(size)
    num = fraction.numerator
    while len(pending) and num:
        digit, num = divmod(num << 32, fraction.denominator)
        draw = rng.integers(2**32, size=len(pending))
        chosen[pending[draw < digit]] = True
        pending = pending[draw == digit]
    return chosen


def pooled_draw(volume, arms, emitter, singles, seed, *, ensemble=None):
    """Direct sharp-word, material-position and gate sampling, with erasures."""
    words, probs, moments = records() if ensemble is None else ensemble
    rng = np.random.default_rng(seed)
    N = 2**20
    gate = tidal.fastest_relay().reshape(2, 11, 3, 3)
    total = 0
    p = chances(volume, arms, singles, moments=moments)[emitter]
    for _ in range(N // 8192):
        selected = [
            words[rng.choice(len(words), 8192, p=np.array(probs[k], float))]
            for k in arms
        ]
        material = np.stack([w[:, :4] for w in selected])
        donors = Q(singles)
        base = donors.numerator // donors.denominator
        n1 = (material == 1).sum((0, 2)) + base
        remainder = donors - base
        if remainder:
            n1 += fractional_donor(rng, remainder, len(n1))
        assert np.max(n1 + (material == 2).sum((0, 2))) <= volume
        n2 = (material == 2).sum((0, 2))
        # Draw the receiving site's load from the population. The sender
        # replaces its resident load by the registered codeword in each slot.
        ticket = rng.integers(volume, size=len(n1))
        target = (ticket < n1).astype(int) + 2 * (ticket >= volume - n2)
        m = selected[emitter][:, -1].copy()
        mask = np.isin(m, PARTIAL)
        m[mask] = rng.choice(PARTIAL, mask.sum())
        # All closed-head classes have identically silent receiver tickets.
        orient = rng.integers(2, size=len(m))
        success = np.where(orient == 0, gate[0, m, 0, target], gate[1, m, target, 0])
        total += int((rng.integers(32, size=len(m)) < success).sum())
    eps = float(np.sqrt(np.log(2 / 0.0001) / (2 * N)))
    assert abs(total / N - float(p)) <= eps
    return dict(
        trials=N, successes=total, p=str(p), simultaneous_per_bank_error=eps, seed=seed
    )


def preparation_routes(board, layout):
    """Lift each addressed compact appointment to the physical external seam.
    Blank full pulses transpose the loads; reversing the two routing strings
    returns the outgoing loads to their own homes and restores the other ports.
    """
    home = tidal.lab.energy.pair.stock.observations.MATTER
    g = tidal.lab.gates()[0]
    for left, right in product(range(3), repeat=2):
        assert (
            g[3 * right + left, 3 * left + right] == 4
            and g[:, 3 * left + right].sum() == 4
        )
    found = []
    for owners, marks in layout:
        for j, (left, right) in enumerate(marks):
            paths = []
            assert (left, right) in board.external and (
                left // 12,
                right // 12,
            ) == owners[j : j + 2]
            for endpoint in (left, right):
                owner = endpoint // 12
                start = 12 * owner + home
                trail = {start: [start]}
                waiting = [start]
                for node in waiting:
                    if node == endpoint:
                        break
                    for q in board.neighbours[node]:
                        if q // 12 == owner and q not in trail:
                            trail[q] = [*trail[node], q]
                            waiting.append(q)
                path = trail[endpoint]
                # The transpose test holds for arbitrary spectator loads.
                # Verify as a permutation of all twelve address labels.
                perm = list(range(12))
                for u, v in zip(path[:-1], path[1:], strict=True):
                    x, y = u % 12, v % 12
                    perm[x], perm[y] = perm[y], perm[x]
                assert perm[endpoint % 12] == home
                for u, v in reversed(list(zip(path[:-1], path[1:], strict=True))):
                    x, y = u % 12, v % 12
                    perm[x], perm[y] = perm[y], perm[x]
                assert perm == list(range(12))
                paths.append(path)
            found.append(paths)
    return found


def packing(board, spacing, limit):
    """Enrol ascending disjoint paths, preserving a connected blank carrier."""
    owner_edge = {(a // 12, b // 12): (a, b) for a, b in board.external}
    adj = np.zeros((board.cells, board.cells), bool)
    for a, b in owner_edge:
        adj[a, b] = adj[b, a] = True
    owner_distance = shortest_path(csr_matrix(adj), unweighted=True).astype(int)
    upward = [
        np.flatnonzero(adj[i] & (np.arange(board.cells) > i)).tolist()
        for i in range(board.cells)
    ]
    paths = [
        (a, b, c, d)
        for a in range(board.cells)
        for b in upward[a]
        for c in upward[b]
        for d in upward[c]
    ]
    used = set()
    chosen = []
    blocked = adj.copy()
    for _ in range(limit):
        available = [p for p in paths if not used.intersection(p)]
        if used and spacing == "dispersed":
            nearest = owner_distance[:, sorted(used)].min(1)
            available.sort(key=lambda p: (-nearest[list(p)].min(), p))
        for p in available:
            trial = blocked.copy()
            seams = []
            for a, b in zip(p[:-1], p[1:], strict=True):
                trial[a, b] = trial[b, a] = False
                seams.append(owner_edge[a, b])
            if connected_components(csr_matrix(trial), directed=False)[0] == 1:
                chosen.append((p, seams))
                used.update(p)
                blocked = trial
                break
        else:
            break
    return chosen
