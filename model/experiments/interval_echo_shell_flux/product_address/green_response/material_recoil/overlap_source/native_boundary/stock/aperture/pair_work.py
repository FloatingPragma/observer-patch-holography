"""Compression work and a mobile ordinary reader of two transported home records."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt

import networkx as nx
import numpy as np
from scipy.linalg import expm

from . import closed_energy as energy
from .closed_insertion import local_spectra
from .measure import check

HERE = energy.HERE
SHOTS = 2**18
DURATION = 128
LEVELS = np.array([0, 1, 8, 9, 32, 33])


def mount(board, home, marks):
    graph = nx.Graph()
    graph.add_edges_from(sorted(board.seams))
    graph.remove_edges_from(marks)
    blocked = set(np.ravel(marks)) | set(home)
    clean = sorted(v for v in graph if v not in blocked and not set(graph[v]) & blocked)
    g = graph.subgraph(clean)
    for left_end in clean:
        for a in sorted(g[left_end]):
            for b in sorted(set(g[a]) - {left_end}):
                if a // 12 != b // 12 or b in g[left_end]:
                    continue
                for r in sorted(set(g[b]) - {a}):
                    path = [left_end, a, b, r]
                    if (
                        g.subgraph(path).number_of_edges() == 3
                        and len(set(board.degree[path])) == 1
                    ):
                        return graph, path
    raise ValueError("No clean equal-degree path")


def route(graph, flat, initial, final, spectators, swap):
    placed = list(initial)
    routes = []
    for j, dest in enumerate(final):
        avoid = set(spectators) | {placed[1 - j]}
        path = nx.shortest_path(graph.subgraph(set(graph) - avoid), placed[j], dest)
        for x, y in zip(path, path[1:]):
            code = 3 * flat[:, x] + flat[:, y]
            new = swap[code]
            flat[:, x], flat[:, y] = new // 3, new % 3
        routes.append(path)
        placed[j] = dest
    return routes


def whole(board, seams, w, flat, h):
    left, right = np.array(board.seams).T
    labels = np.zeros((len(w), len(left)), int)
    for j, seam in enumerate(seams):
        labels[:, board.seams.index(seam)] = w[:, 4 + j]
    code = 9 * labels + 3 * flat[:, left] + flat[:, right]
    return h.diagonal()[code].sum(1) + 32 * ((flat**2) @ (6 - board.degree))


def read_compression(w, weights, rng, spectra):
    pair = w[:, :2]
    settings = [
        pair,
        np.column_stack((pair[:, 0], np.zeros(len(w), int))),
        np.column_stack((np.zeros(len(w), int), pair[:, 1])),
    ]
    raw = []
    for ports in settings:
        draw = rng.choice(len(w), SHOTS, p=weights)
        level = spectra[ports[draw, 0], ports[draw, 1], rng.integers(2, size=SHOTS)]
        raw.append(np.bincount(level, minlength=34)[LEVELS])
    means = np.asarray(raw) @ LEVELS / SHOTS
    return dict(
        counts=np.asarray(raw).tolist(), work=float(means[0] - means[1] - means[2])
    )


def root_band(square):
    # Exact enclosure of sqrt(square) for a nonnegative integer.
    scale = 10**12
    x = isqrt(square * scale**2)
    return F(x, scale), F(x + 1, scale)


def reader(t, h, degree, n, m):
    # Calculate the two native diagonal insertion energies relative to the
    # unchanged sources. The diagonal of 4H on a blank bond is tabulated.
    def insert(x):
        return int(
            h[3 * t + x, 3 * t + x]
            - h[x, x]
            + (degree - 1) * h[3 * t, 3 * t]
            + 32 * t * t * (6 - degree)
        )

    left, right = insert(n), insert(m)
    hop = int(h[3 * t, t])
    assert hop == -2
    gap_sq = (left - right) ** 2 + 4 * hop**2
    low, high = root_band(gap_sq)
    level = ((left + right - high) / 8, (left + right - low) / 8)
    gap = np.sqrt(gap_sq) / 4
    reach = (4 * hop**2 / gap_sq) * (1 - np.sinc(DURATION * gap / np.pi)) / 2
    matrix = np.array([[left, hop], [hop, right]]) / 4
    np.testing.assert_allclose(
        np.linalg.eigvalsh(matrix)[0], float(level[0]), atol=1e-12
    )
    # Check both the propagation and the uniform-time average. The sampler
    # below draws a fresh exposure and a binary port record on each run.
    for time in (0.25, 4.0):
        exact = (4 * hop**2 / gap_sq) * np.sin(time * gap / 2) ** 2
        np.testing.assert_allclose(
            abs(expm(-1j * time * matrix)[1, 0]) ** 2, exact, atol=1e-12
        )
    return level, gap, 4 * hop**2 / gap_sq, reach


def scaleband(band, p):
    return tuple(p * x for x in band)


def sub(a, b):
    return a[0] - b[1], a[1] - b[0]


def ratio(a, b):
    assert b[1] < 0
    ends = [x / y for x in a for y in b]
    return [float(min(ends)), float(max(ends))]


def run():
    _, w, cells, *_ = energy.pair.setup()
    board, home, marks, *_ = energy.returned(w)
    h, h0 = energy.pair.stock.kernels()[2:]
    # Full blank pi pulse under the actual stiffness law is an exact SWAP.
    swap = np.arange(9).reshape(3, 3).T.ravel()
    np.testing.assert_allclose(
        expm(-1j * np.pi * h[:9, :9] / 4), np.eye(9)[swap], atol=1e-13
    )
    g, path = mount(board, home, marks)
    left_end, a, b, r = path
    initial = np.zeros((len(w), board.n), int)
    initial[:, home] = w[:, :4]
    flat = initial.copy()
    routes = route(g, flat, home[:2], [a, r], home[2:], swap)
    far = whole(board, marks, w, flat, h)
    flat[:, [r, b]] = flat[:, [b, r]]
    near = whole(board, marks, w, flat, h)
    x, y = w[:, :2].T
    np.testing.assert_array_equal(flat[:, [a, b]], w[:, :2])
    np.testing.assert_array_equal((x - y) ** 2, (x - 1) ** 2 + 3 * (y - 1) ** 2)
    np.testing.assert_array_equal(x != y, x != 1)

    def contact(k):
        d = k.diagonal()
        return d[3 * x + y] - d[3 * x] - d[y]

    np.testing.assert_array_equal(near - far, contact(h))
    flat[:, [r, b]] = flat[:, [b, r]]
    for p in reversed(routes):
        for i, j in zip(p[:0:-1], p[-2::-1]):
            flat[:, [i, j]] = flat[:, [j, i]]
    np.testing.assert_array_equal(initial, flat)
    mobile_routes = route(g, flat, home[:2], [left_end, r], home[2:], swap)
    np.testing.assert_array_equal(flat[:, [left_end, r]], w[:, :2])
    # Check the projected mobile matrix diagonals against the entire mounted
    # board on every accepted word, for both probes and both placements.
    before = whole(board, marks, w, flat, h)
    degree = int(board.degree[a])
    for t in (1, 2):
        for site, neighbor in ((a, x), (b, y)):
            flat[:, site] = t
            diag = h.diagonal()
            expected = (
                diag[3 * t + neighbor]
                - diag[neighbor]
                + (degree - 1) * diag[3 * t]
                + 32 * t * t * (6 - degree)
            )
            np.testing.assert_array_equal(
                whole(board, marks, w, flat, h) - before, expected
            )
            flat[:, site] = 0

    original = np.full(len(w), F(1, len(w)), object)
    weights = dict(original=original, repaired=energy.least_change(cells) / len(w))
    levels = local_spectra()
    spectra = np.array([[levels[i, j] for j in range(3)] for i in range(3)])
    rng = np.random.default_rng(1338627)
    arms, joints, work = {}, {}, {}
    for key, p in weights.items():
        assert p @ x == p @ y == 1 and sum(p) == 1
        v, f = p @ ((x - 1) ** 2), p @ ((y - 1) ** 2)
        cost = (p @ contact(h)) / 4
        assert cost == -17 + F(65, 4) * f + F(3, 4) * v
        work[key] = cost
        arms[key] = dict(
            partner_variance=str(v),
            middle_variance=str(f),
            work=str(cost),
            exchange_only_work=str((p @ contact(h0)) / 4),
            read=read_compression(w, np.array(p, float), rng, spectra),
        )
        joint = np.array(
            [[sum(p[(x == i) & (y == j)]) for j in range(3)] for i in range(3)], object
        )
        joints[key] = joint
        joints[key + "_unpaired"] = np.outer(joint.sum(1), joint.sum(0))
    bound = F(65, 4) * F(1, 14) + F(3, 4) * (F(5, 28) - F(11, 28))
    assert bound == 1 and work["repaired"] - work["original"] >= bound
    delta = arms["repaired"]["read"]["work"] - arms["original"]["read"]["work"]
    # Allocate 0.5% error to six spectral means, 0.5% to 12 arrival banks.
    width = 6 * 33 * np.sqrt(np.log(2 * 6 / 0.005) / (2 * SHOTS))
    assert (
        delta - width > 0
        and abs(delta - float(work["repaired"] - work["original"])) <= width
    )
    for key, pair in (("vacuum", (0, 0)), ("ordinary", (1, 0))):
        joints[key] = np.full((3, 3), F(), object)
        joints[key][pair] = F(1)
    mobile = {}
    arrival_width = np.sqrt(np.log(2 * 12 / 0.005) / (2 * SHOTS))
    bands = []
    for t in (1, 2):
        mx = [reader(t, h, degree, i, j) for i in range(3) for j in range(3)]
        gaps, amplitudes, means = np.array([s[1:] for s in mx]).T
        row, levels = {}, {}
        for key, joint in joints.items():
            p = joint.ravel()
            lo, hi = zip(*(scaleband(v[0], q) for v, q in zip(mx, p, strict=True)))
            levels[key] = sum(lo), sum(hi)
            rngidx = rng.choice(9, SHOTS, p=np.asarray(p, float))
            time = rng.uniform(0, DURATION, SHOTS)
            k = int(
                np.sum(
                    rng.random(SHOTS)
                    < amplitudes[rngidx] * np.sin(gaps[rngidx] * time / 2) ** 2
                )
            )
            expected = float(np.asarray(p, float) @ means)
            assert abs(k / SHOTS - expected) <= arrival_width
            row[key] = dict(
                lower_level=[float(v) for v in levels[key]],
                expected_arrival=expected,
                arrivals=k,
            )
        ordinary = sub(levels["ordinary"], levels["vacuum"])
        row["shifts_in_ordinary_units"] = {
            tag: ratio(sub(levels[right], levels[left]), ordinary)
            for tag, left, right in (
                ("repair", "original", "repaired"),
                ("original_pairing", "original_unpaired", "original"),
                ("repaired_pairing", "repaired_unpaired", "repaired"),
            )
        }
        bands.append(row["shifts_in_ordinary_units"]["repair"])
        mobile[str(t)] = row
    recorded = json.loads((HERE / "fluctuation_balls.json").read_bytes())["arms"]
    u, v = (recorded[k]["measured_variance_bands"] for k in ("original", "repaired"))
    # The cited source banks assign simultaneous 99% bounds to the port
    # variances. Propagate both endpoints before consulting compression data.
    source_band = [
        (65 * (v[2][i] - u[2][1 - i]) + 3 * (v[1][i] - u[1][1 - i])) / 4 for i in (0, 1)
    ]
    assert 0 < source_band[0] < work["repaired"] - work["original"] < source_band[1]
    return dict(
        source_count_prediction_99=source_band,
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "PAIR_WORK.md", HERE / "fluctuation_balls.json")
        },
        home=home,
        path=path,
        degree=degree,
        compression_routes=routes,
        mobile_routes=mobile_routes,
        shots=SHOTS,
        levels=LEVELS.tolist(),
        arms=arms,
        required_increment_lower=str(bound),
        canonical_increment=str(work["repaired"] - work["original"]),
        measured_increment_99=[delta - width, delta + width],
        arrival_simultaneous_halfwidth=arrival_width,
        mobile=mobile,
        common_repair_scalar_excluded=max(v[0] for v in bands)
        > min(v[1] for v in bands),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "pair_work.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        check(json.loads(dst.read_bytes()), result)
    print(
        result["path"], result["canonical_increment"], result["measured_increment_99"]
    )
    for t, row in result["mobile"].items():
        print(t, row["shifts_in_ordinary_units"])
