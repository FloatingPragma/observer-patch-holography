"""The same native scalar triangle at two owners, with complete ternary echoes."""

from functools import cache
from itertools import combinations, pairwise
from math import atan2, pi

import numpy as np

from model.executable_model import federation
from model.experiments.interval_echo_gauge.exterior_register import material
from .. import measure as scalar
from .sampling import multinomial


@cache
def transport():
    graph, previous = material.shared_report.layout()
    board = federation(0)
    owners = {"source": previous["alice"], "receiver": previous["bob"]}
    local = {
        name: [12 * owner + i for i in scalar.PORTS] for name, owner in owners.items()
    }
    read_paths = {}
    reference_loads = {}
    for name, owner in owners.items():
        a = graph.subgraph(range(12 * owner, 12 * (owner + 1)))
        triangle = local[name]
        assert all(a.has_edge(i, j) for i, j in combinations(triangle, 2))
        assert all(board.degree[i] == 5 for i in triangle)
        standards = previous["readers"]["alice" if name == "source" else "bob"]
        for val, (_, ref) in zip((0, 2), standards, strict=True):
            reference_loads[ref] = val
        read_paths[name] = [
            [(material.shared_report.path(a, p, g), r) for g, r in standards]
            for p in triangle
        ]
    transit = graph.copy()
    transit.remove_nodes_from(reference_loads)
    paths = [
        material.shared_report.path(transit, a, b)
        for a, b in zip(local["source"], local["receiver"], strict=True)
    ]
    perm, cut = list(range(board.n)), list(range(board.n))
    for path in paths:
        moves = list(pairwise(path))
        for u, v in moves + moves[:-1][::-1]:
            perm[u], perm[v] = perm[v], perm[u]
            if u // 12 == v // 12:
                cut[u], cut[v] = cut[v], cut[u]
    expected = list(range(board.n))
    for trail in paths:
        u, v = trail[0], trail[-1]
        expected[u], expected[v] = expected[v], expected[u]
    assert perm == expected
    lookup = material.shared_report.certify()
    raw = []
    for j, p in enumerate(local["source"]):
        stock = np.zeros(board.n, dtype=int)
        stock[p] = 1
        for r, v in reference_loads.items():
            stock[r] = v
        row = [
            material.read(stock, read_paths["source"], lookup),
            material.read(stock[perm], read_paths["receiver"], lookup),
            material.read(stock[cut], read_paths["receiver"], lookup),
        ]
        expected_bits = "".join("00" if i == j else "01" for i in range(3))
        assert row == [expected_bits, expected_bits, "010101"]
        raw.append(row)
    labels = [r[0] for r in raw]
    mapped = [labels.index(r[1]) for r in raw]
    # With no load at the remote triangle any analyzer pulse fixes its blank
    # vacuum. A cut link therefore keeps its malformed word under all settings.
    law = scalar.candidate.native_law
    ex, con = law.local_term_parts()
    np.testing.assert_allclose(
        (ex + scalar.candidate.KAPPA * con)[:, law.SRC.index((0, 0))], 0, atol=2e-13
    )
    return dict(
        owners=owners,
        triangles=local,
        transfer_paths=paths,
        connected_columns=mapped,
        ternary_read_paths=read_paths,
        address_clicks=raw,
        path_contacts=sum(2 * len(p) - 3 for p in paths),
    )


@cache
def state(doublet):
    u, d = doublet
    norm = np.hypot(abs(u), abs(d))
    word = scalar.ry(0, 2, -2 * atan2(1, norm))
    if norm:
        word += scalar.ry(0, 1, -2 * atan2(abs(d), abs(u)))
    word += [scalar.phase(i, np.angle(z)) for i, z in enumerate(doublet)]
    prepared = scalar.execute(word)[:, 0]
    expected = np.asarray([u, d, 1], complex) / np.sqrt(1 + norm**2)
    np.testing.assert_allclose(prepared, expected, atol=3e-13)
    return prepared


@cache
def read(doublet, target=None, axis=None, site="source"):
    wave = state(doublet)
    if site == "receiver":
        transfer = np.eye(3)[:, transport()["connected_columns"]]
        wave = transfer @ wave
    elif site == "cut":
        return np.array([0.0, 0.0, 0.0, 1.0])
    elif site != "source":
        raise ValueError("Unknown triangle read location")
    if target is not None:
        phi = 0 if axis == "X" else -pi / 2
        read_word = scalar.reader(("pair", target, 2, phi))
        rotated = scalar.execute(read_word) @ wave
        check = scalar.reading_vector(("pair", target, 2, phi)) @ wave
        np.testing.assert_allclose(abs(rotated) ** 2, abs(check) ** 2, atol=2e-13)
        wave = rotated
    return abs(wave) ** 2


@cache
def coin(fair):
    dwell = pi / 2 if fair else pi / 3
    wave = scalar.execute(scalar.rx(0, 1, dwell))[:, 0]
    expected = [0.5, 0.5, 0] if fair else [0.75, 0.25, 0]
    np.testing.assert_allclose(abs(wave) ** 2, expected, atol=2e-13)
    return abs(wave) ** 2


SHOTS = 2**18
SEED = 686900
GRID = 2**36
# Every individual variable used for inference is bounded. Each joint raw
# bank has at most 81 cells. The maximum TV rounding error is 81 / GRID.
UNIT_ERROR = np.sqrt(np.log(2 * 1024 / 0.01) / (2 * SHOTS)) + 81 / GRID


def acquire(rng, probability):
    p = np.asarray(probability)
    np.testing.assert_allclose(p.sum(), 1, atol=1e-12)
    weights = np.floor(p.ravel() * GRID + 0.5).astype(np.int64)
    weights[np.argmax(weights)] += GRID - weights.sum()
    assert min(weights) >= 0
    assert np.abs(weights / GRID - p.ravel()).sum() / 2 <= 81 / GRID
    counts = multinomial(rng, SHOTS, weights).reshape(p.shape)
    raw = [
        [*(int(i) for i in indices), int(counts[indices])]
        for indices in zip(*np.nonzero(counts))
    ]
    assert sum(row[-1] for row in raw) == SHOTS
    return raw


def joint(*probabilities):
    result = 1.0
    for p in probabilities:
        result = np.multiply.outer(result, p)
    return result


def sample_mean(raw, score):
    return sum(score(*r[:-1]) * r[-1] for r in raw) / SHOTS


def ratio_interval(numerator, denominator, width):
    low, high = denominator - UNIT_ERROR, denominator + UNIT_ERROR
    if low <= 0:
        return None
    endpoints = [
        x / y for x in (numerator - width, numerator + width) for y in (low, high)
    ]
    return [min(endpoints), max(endpoints)]
