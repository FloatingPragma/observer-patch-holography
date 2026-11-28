"""Test a midpoint potential and a stationary ordinary load on an exterior counter."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import factorial
import json

import networkx as nx
import numpy as np
from scipy import sparse
import sympy as sp

from model.executable_model import federation
from model.experiments.interval_echo_einstein_station import shared_variance, volume
from model.experiments.mobile_memory.calibrated import geometry
from . import initial_acceleration as meter
from . import midpoint_bulk as bulk

routes = meter.routes
iv, HERE = routes.iv, routes.HERE
SEED = 28122812


def compact_mode(lap):
    local = sp.Matrix(lap[:12, :12].toarray().astype(int))
    gateways = sorted(set(lap[12:, :12].nonzero()[1].tolist()))
    g = sp.eye(12)[gateways, :]
    rows = []
    for _ in range(12):
        rows.append(g)
        g = g * local
    basis = sp.Matrix.hstack(*sp.Matrix.vstack(*rows).nullspace())
    restricted = (basis.T * basis).inv() * basis.T * local * basis
    assert local * basis == basis * restricted
    values = sorted(restricted.eigenvals().items())
    modes = (restricted - values[0][0] * sp.eye(restricted.rows)).nullspace()
    assert len(modes) == 1
    vector = basis * modes[0]
    integer = np.array(vector * sp.ilcm(*[v.q for v in vector]), dtype=int).ravel()
    if integer[np.flatnonzero(integer)[0]] < 0:
        integer *= -1
    full = np.zeros(lap.shape[0], int)
    full[:12] = integer
    np.testing.assert_array_equal(lap @ full, int(values[0][0]) * full)
    return full / np.linalg.norm(full), {
        "gateways": gateways,
        "home_amplitudes_before_normalization": integer.tolist(),
        "compact_laplacian_spectrum": [[int(z), n] for z, n in values],
    }


def stationarity(h, psi, rng):
    site = next(i for i in np.flatnonzero(psi) if h[i, i] == 2)
    control = np.eye(1, len(psi), site).ravel()
    starts = np.column_stack([psi, control])
    e, u = np.linalg.eigh(h.toarray())
    coeff = u.T @ starts
    state = {}
    for time in (0, 1, 4, 16, 64):
        path = u @ (np.exp(-1j * e * time)[:, None] * coeff)
        probs = (abs(path) ** 2).reshape(-1, 12, 2).sum(axis=1)
        state[str(time)] = [routes.lab.counts(p / p.sum(), rng) for p in probs.T]
    delay = 1 / 32768
    phase = (abs(coeff) ** 2).T @ np.sin((e + 40) * delay)
    err = (40 + max(e) + 1e-10) ** 3 * delay**2 / 6
    energy = []
    for signal in phase:
        tape = routes.lab.counts([(1 + signal) / 2, (1 - signal) / 2], rng)
        f = tape[0] / routes.lab.SHOTS
        eps = 2 * routes.reg.src.probability_error(np.asarray(f), routes.lab.SHOTS)
        band = routes.reg.band(2 * f - 1 - eps, 2 * f - 1 + eps) / delay + iv.mpf(
            [-err, err]
        )
        energy.append({"counts": tape, "energy": routes.bounds(band)})
    return {
        "times": state,
        "energy": energy,
        "single_port_control": int(site),
        "bypass_time": delay,
    }


def potential_table(board, ball, graph, points):
    hop = np.asarray(
        nx.floyd_warshall_numpy(graph, nodelist=list(range(board.cells))), int
    )
    dots = [
        iv.mpf(1),
        iv.sqrt(5) / 3,
        iv.mpf(1) / 3,
        -iv.mpf(1) / 3,
        -iv.sqrt(5) / 3,
        iv.mpf(-1),
    ]
    np.testing.assert_allclose(
        points @ points.T, np.array([float(x.mid) for x in dots])[hop], atol=1e-12
    )

    def dot(a, b):
        return dots[hop[a, b]]

    def word(v, i, j):
        return -1 / iv.sqrt(1 + (1 + dot(i, j)) / 2 - dot(v, i) - dot(v, j))

    pairs = list(combinations(ball, 2))
    edges = sorted(
        (a, b)
        for edge in board.external
        for a, b in (edge, edge[::-1])
        if a // 12 not in ball
        and b // 12 not in ball
        and hop[0, a // 12] + 1 == hop[0, b // 12]
    )
    result = []
    for a, b in edges:
        u, v = a // 12, b // 12
        inc = [word(v, i, j) - word(u, i, j) for i, j in pairs]
        ordinary = 1 / iv.sqrt(2 - 2 * dot(0, u)) - 1 / iv.sqrt(2 - 2 * dot(0, v))
        result.append(
            {
                "ports": [int(a), int(b)],
                "owners": [u, v],
                "ordinary": routes.bounds(ordinary),
                "pair_increments": list(map(routes.bounds, inc)),
                "mean": routes.bounds(sum(inc) / len(inc)),
            }
        )
    return result, pairs


def superpose(psi, x, y, a, b, q):
    def part(v):
        p = psi[x] * (y == v)
        return (p + psi[y] * (x == v)) / np.sqrt(2) if q == 1 else p

    left, right = part(a), part(b)
    z = np.column_stack([left + right, left - right]) / np.sqrt(2)
    np.testing.assert_allclose(z.T @ z, np.eye(2), atol=1e-14)
    return z


def species_controls(board, one, adj, psi, a, b, common, rng):
    score = (np.arange(board.n) // 12 == a // 12).astype(int)
    start = np.zeros((board.n, 2))
    start[a] = 1 / np.sqrt(2)
    start[b] = np.array([1, -1]) / np.sqrt(2)
    local = 12 * (a // 12) + int(np.flatnonzero(common[a % 12])[0])
    assert adj[a, local] and not adj[b, local]
    records = {}
    terms = meter.contact_terms()
    for q in (1, 2):
        raw = meter.read(one, start, score, rng)
        h, x, y = meter.two_loads(one, adj, terms[q], q)
        out = {}
        for name, source in (
            ("settled", psi),
            ("neighbor", np.eye(1, board.n, local).ravel()),
        ):
            obs = meter.read(
                h, superpose(source, x, y, a, b, q), score[x] + score[y], rng
            )

            def split(r):
                p, m = r["phases"]
                return (
                    iv.mpf(p["held_acceleration"]) - iv.mpf(m["held_acceleration"])
                ) / 2, (p["exact_acceleration"] - m["exact_acceleration"]) / 2

            z, exact = split(obs)
            ref, base = split(raw)
            expected = 0 if name == "settled" else -terms[q] / 2
            assert abs(exact - base - expected) < 1e-11 and expected in (z - ref)
            out[name] = {
                "raw": obs,
                "held": routes.bounds(z - ref),
                "predicted_native": expected,
            }
        records[str(q)] = {"empty": raw, **out}
    return records, start, score


def source_fields(board, cal, pairs, arms, starts, score, rng):
    active = routes.graph(board, cal)
    step = 1 / 48
    times = np.arange(5) * step
    table = np.zeros((2, 5, 2))
    for i, j in pairs:
        for k, w in np.ndindex(9, 9):
            adj, degree = active(i, j, k + 1, w + 1)
            h = (sparse.diags(degree, dtype=float) - adj) * 0.5
            norm = max(degree)
            power = starts
            path = np.broadcast_to(starts, (5, *starts.shape)).astype(complex).copy()
            for n in range(1, 20):
                power = h @ power
                path += ((-1j * times) ** n / factorial(n))[:, None, None] * power
            assert (
                np.exp(4 * step * norm) * (4 * step * norm) ** 20 / factorial(20)
                < 1e-23
            )
            prob = (abs(path[:, score == 1]) ** 2).sum(axis=1)
            table += arms[:, k, w, None, None] * prob / len(pairs)
    result = {}
    weights = list(map(float, meter.WEIGHTS))
    err = (
        (2 * norm) ** 10
        * step**8
        * sum(abs(w) * j**10 for j, w in enumerate(weights))
        / factorial(10)
    )
    for q in (1, 2):
        held = iv.mpf(0)
        counts = []
        for arm, sgn in zip(table, [-1, 1], strict=True):
            tapes = []
            for j in range(1, 5):
                phases = []
                for sign, prob in zip((1, -1), arm[j], strict=True):
                    tape = routes.lab.counts([prob, 1 - prob], rng)
                    phases.append(tape)
                    p = tape[0] / routes.lab.SHOTS
                    eps = routes.reg.src.probability_error(
                        np.asarray(p), routes.lab.SHOTS
                    )
                    held += (
                        sgn
                        * sign
                        * weights[j]
                        * routes.reg.band(p - eps, p + eps)
                        / (2 * step**2)
                    )
                tapes.append(phases)
            counts.append(tapes)
        # All four equal launch-time probabilities cancel; two units of
        # stencil error cover the half-difference of four accelerations.
        held += iv.mpf([-2 * err - 1e-8, 2 * err + 1e-8])
        result[str(q)] = {"counts": counts, "held_conditional": routes.bounds(held)}
    return result


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    board = federation(0)
    graph, points = geometry._atlas([board])
    prior = routes.shells.dist.read("word_routes")
    edges, pairs = potential_table(board, prior["source_ball"], graph[0], points[0])
    a, b = edges[0]["ports"]
    code, common = routes.codebook()
    cal, _ = routes.calibrate(code)
    adj, deg = routes.graph(board, cal)(-1, -1, 0, 0)
    lap = sparse.diags(deg, dtype=float) - adj
    psi, mode = compact_mode(lap)
    one = lap * 0.5
    mode["observer_reads"] = stationarity(one, psi, rng)
    controls, starts, score = species_controls(board, one, adj, psi, a, b, common, rng)
    src = volume.radar_source_and_hand(board)[-1]
    laws, _ = shared_variance.source_words(src)
    _, fresh, shared = np.rint(2592 * laws).astype(int)
    energy = np.array(src["source_word_energies"])
    cost = (energy[:, None] + energy).ravel()
    total = cost[:, None] + cost
    for e in np.unique(total):
        assert fresh[total == e].sum() == shared[total == e].sum()
    overlap = source_fields(
        board, cal, pairs, np.stack([fresh, shared]) / 2592, starts, score, rng
    )
    base = routes.shells.dist.read("resolved_shells")
    word = [F(base["word_event_probabilities"][k]) for k in ("fresh", "shared")]
    price = iv.mpf(base["conditional_entropy_slope"])
    local = [
        bulk.native_read(F(2, board.cells * (board.cells - 1)), price, word, rng)
        for _ in pairs
    ]
    qsum = sum(iv.mpf(v["charge"]) for v in local)
    measured = (
        sum(
            iv.mpf(row["charge"]) * iv.mpf(d)
            for row, d in zip(local, edges[0]["pair_increments"], strict=True)
        )
        / qsum
    )
    for d in overlap.values():
        d["held_per_all_starts"] = routes.bounds(
            iv.mpf(d["held_conditional"])
            * routes.shells.interval(F(prior["pair_chance"]))
        )
        d["allowed_conditional_potential_coupling"] = routes.bounds(
            2 * iv.mpf(d["held_conditional"]) / measured
        )
    return {
        "seed": SEED,
        "shots": routes.lab.SHOTS,
        "source_starts": routes.shells.STARTS,
        "potential_edges": edges,
        "pairs": [list(p) for p in pairs],
        "individual_pair_reads": local,
        "potential_difference_from_records": routes.bounds(measured),
        "ordinary_mode": mode,
        "detector_controls": controls,
        "overlap": overlap,
        "overlap_step": "1/48",
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "SETTLED_MIDPOINT_FORCE.md",
                "midpoint_bulk.json",
                "word_routes.json",
                "initial_acceleration.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    path = HERE / "settled_midpoint_force.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    print(data["potential_edges"][0], data["potential_difference_from_records"])
    print({q: a["held_conditional"] for q, a in data["overlap"].items()})
