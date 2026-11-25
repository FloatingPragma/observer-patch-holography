"""Return-time geometry of the sharp two-load words at four exterior portions of a patch."""

import json
from dataclasses import dataclass
from itertools import combinations, permutations
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq

from .. import measure as source

HERE = Path(__file__).resolve().parent
STATIONS = 4
EDGES = tuple(combinations(range(STATIONS), 2))
RAYS = tuple(permutations(range(STATIONS), 2))
TICKS = 24
MEAN_TIMEOUT = 4096


def clock_scores(transition, occupied, cap=TICKS):
    """Capped A-only -> B-only -> A-only first-return means, by initial word."""
    starts = np.column_stack([occupied[:, a] & ~occupied[:, b] for a, b in RAYS])
    opposite = np.column_stack([occupied[:, b] & ~occupied[:, a] for a, b in RAYS])
    outbound = np.zeros(starts.shape)
    homeward = np.zeros(starts.shape)
    for _ in range(cap):
        # Clock counts each interval up to the cap, including the return interval.
        both = np.column_stack(
            [np.where(opposite, homeward, outbound), np.where(starts, 0, homeward)]
        )
        predicted = 1 + transition @ both
        outbound, homeward = np.hsplit(predicted, 2)
    return starts, outbound * starts * (source.EXPOSURE / 2)


def discounted_pair(transition, occupied_a, occupied_b, mean_wait):
    """Exact linear systems for a return before an independent geometric timeout."""
    z = 1 - 1 / mean_wait
    targets = (occupied_a & ~occupied_b, occupied_b & ~occupied_a)
    matrices = []
    for target in targets:
        rem = ~target
        matrix = np.eye(sum(rem)) - z * transition[np.ix_(rem, rem)]
        matrices.append(cho_factor(matrix, check_finite=False))

    def hit(i, boundary):
        mask = targets[i]
        x = np.empty(len(mask))
        x[mask] = boundary
        b = z * transition[np.ix_(~mask, mask)] @ boundary
        x[~mask] = cho_solve(matrices[i], b, check_finite=False)
        return x

    ret = []
    for a, b in ((0, 1), (1, 0)):
        return_a = hit(a, np.ones(sum(targets[a])))
        success = hit(b, return_a[targets[b]])
        cost = mean_wait * (1 - success) * source.EXPOSURE / 2
        ret.append((targets[a], targets[a] * cost, success, return_a))
    return ret


def discounted_scores(transition, occupied, mean_wait):
    result = {}
    for a, b in EDGES:
        runs = discounted_pair(transition, occupied[:, a], occupied[:, b], mean_wait)
        result[a, b], result[b, a] = runs
    return tuple(np.column_stack([result[r][i] for r in RAYS]) for i in (0, 1))


def distances(pop, starts, scores):
    directed = (pop @ scores) / (pop @ starts)
    return np.array(
        [
            (directed[RAYS.index((a, b))] + directed[RAYS.index((b, a))]) / 2
            for a, b in EDGES
        ]
    )


def geometry(edges):
    # Cayley/Gram reconstruction with the first actual rim port as origin.
    squared = np.zeros((4, 4), dtype=edges.dtype)
    for (a, b), d in zip(EDGES, edges, strict=True):
        squared[a, b] = squared[b, a] = d * d
    gram = (squared[0, 1:, None] + squared[0, None, 1:] - squared[1:, 1:]) / 2
    volume = np.sqrt(np.linalg.det(gram)) / 6
    areas = []
    for i, j, k in combinations(range(4), 3):
        a, b, c = squared[i, j], squared[i, k], squared[j, k]
        areas.append(np.sqrt(4 * a * b - (a + b - c) ** 2) / 4)
    return np.array([sum(areas), volume]), gram


def differentials(pop, starts, scores):
    enrollment = pop @ starts
    means = (pop @ scores) / enrollment
    partials = (scores - starts * means) / enrollment
    de = np.column_stack(
        [
            (partials[:, RAYS.index((a, b))] + partials[:, RAYS.index((b, a))]) / 2
            for a, b in EDGES
        ]
    )
    edges = distances(pop, starts, scores)
    jac = np.array(
        [geometry(edges + 1e-20j * col)[0].imag / 1e-20 for col in np.eye(6)]
    )
    return de @ jac, means


@dataclass
class Survey:
    board_cells: int
    cells: list
    words: list
    rim: list
    h: np.ndarray
    transition: np.ndarray
    occupied: np.ndarray
    cut_count: np.ndarray
    entropy_scores: np.ndarray
    uniform: np.ndarray
    starts: np.ndarray
    scores: np.ndarray


def survey(level):
    board = source.federation(level)
    adjacent = sorted(
        {
            b // 12 if a // 12 == 0 else a // 12
            for a, b in board.external
            if 0 in (a // 12, b // 12)
        }
    )
    cells = [0, adjacent[0]]
    ports, contacts, border, words = source.apparatus(board, cells)
    assert len(ports) == 24 and len(border) == 4
    h = source.generator(board, words, contacts)
    prob, _ = source.sampling_kernel(h)
    occupied = np.array([[a in word for a in sorted(border)] for word in words])
    k = occupied.sum(1)
    # The common log(total words) cancels from all trace-zero directions.
    entropy = -np.log([comb(len(border), int(i)) for i in k])
    pop = np.full(len(words), 1 / len(words))
    starts, scores = discounted_scores(prob, occupied, MEAN_TIMEOUT)
    return Survey(
        board.n // 12,
        cells,
        words,
        sorted(border),
        h,
        prob,
        occupied,
        k,
        entropy,
        pop,
        starts,
        scores,
    )


def survey_directions(s, grad):
    E = s.h.diagonal()
    k = s.cut_count
    uniform = np.array([(k == j) / sum(k == j) for j in range(3)])
    directions = dict(
        zero_to_one=uniform[1] - uniform[0],
        zero_to_two=uniform[2] - uniform[0],
        one_to_two=uniform[2] - uniform[1],
    )
    sectors = uniform @ E
    e = np.cross(np.ones(3), sectors)
    e *= np.sign(e[np.flatnonzero(e)[0]])
    directions["matched_energy"] = (e / e[e > 0].sum()) @ uniform
    idx = np.flatnonzero(k == 0)
    lo = idx[np.argmin(E[idx])]
    hi = idx[np.argmax(E[idx])]
    joint = np.array([0.5, -1.0, 0.5]) @ uniform
    c = -(joint @ E) / (E[hi] - E[lo])
    joint[hi] += c
    joint[lo] -= c
    directions["matched_energy_and_cut_load"] = joint / joint[joint > 0].sum()
    visible = k > 0
    hidden = visible / visible.sum() - (~visible) / (~visible).sum()
    c = -(hidden @ E) / (E[hi] - E[lo])
    hidden[hi] += c
    hidden[lo] -= c
    directions["hidden_interior_frequency"] = hidden
    pairs = []
    for cls, energy in sorted(set(zip(k, E, strict=True))):
        ids = np.flatnonzero((k == cls) & (E == energy))
        v = grad[ids, 1]
        pairs.append((v.max() - v.min(), ids[v.argmin()], ids[v.argmax()]))
    _, a, b = max(pairs)
    comp = np.zeros(len(k))
    comp[a] = -1
    comp[b] = 1
    return directions, comp, [s.words[a], s.words[b]]


def responses(s, starts, scores):
    grad, means = differentials(s.uniform, starts, scores)
    directions, comp, controller = survey_directions(s, grad)
    ctrl = comp @ grad
    readings = {}
    for key, raw in directions.items():
        coeff = -(raw @ grad[:, 1]) / ctrl[1]
        t = raw + coeff * comp
        readings[key] = dict(
            area=float(t @ grad[:, 0]),
            volume=float(t @ grad[:, 1]),
            entropy=float(t @ s.entropy_scores),
            energy=float(t @ s.h.diagonal()),
            cut=float(t @ s.cut_count),
            control=float(coeff),
        )
    area = readings["zero_to_one"]["area"]
    eta = -readings["zero_to_one"]["entropy"] / area
    for v in readings.values():
        v["residual"] = v["entropy"] + eta * v["area"]
    edges = distances(s.uniform, starts, scores)
    base, gram = geometry(edges)
    return dict(
        gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
        edges=edges.tolist(),
        means=means.tolist(),
        baseline=base.tolist(),
        volume_control=controller,
        controller_slopes=ctrl.tolist(),
        directions=readings,
        calibrated_eta=eta,
    )


def fifth(s, edges):
    port = min(v for word in s.words for v in word if v not in s.rim)
    a = np.array([port in w for w in s.words])
    ds = []
    for j in range(4):
        rays = discounted_pair(s.transition, a, s.occupied[:, j], MEAN_TIMEOUT)
        ds.append(sum(s.uniform @ c / (s.uniform @ h) for h, c, _, _ in rays) / 2)
    gram = geometry(np.asarray(edges))[1]
    y = (np.diag(gram) + ds[0] ** 2 - np.square(ds[1:])) / 2
    err = ds[0] ** 2 - y @ np.linalg.solve(gram, y)
    return dict(
        port=port,
        ranges=ds,
        extra_squared_length=float(err),
        relative_extra_squared_length=float(err / max(ds) ** 2),
    )


def run():
    rows = []
    for level in range(3):
        s = survey(level)
        timed = responses(s, s.starts, s.scores)
        all_attempts = responses(
            s,
            np.ones(s.starts.shape),
            s.scores + (~s.starts) * MEAN_TIMEOUT * source.EXPOSURE / 2,
        )
        short_starts, short_scores = clock_scores(s.transition, s.occupied, TICKS)
        short = responses(s, short_starts, short_scores)
        herald = s.uniform @ s.starts
        rows.append(
            dict(
                host=s.board_cells,
                cells=s.cells,
                stations=s.rim,
                herald_probability=herald.tolist(),
                timed=timed,
                all_attempts=all_attempts,
                short_clock=short,
                held_fifth=fifth(s, timed["edges"]),
            )
        )
    return dict(
        short_clock_cap=TICKS,
        mean_timeout_ticks=MEAN_TIMEOUT,
        physical_tick=source.EXPOSURE,
        rows=rows,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
