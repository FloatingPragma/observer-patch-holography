"""Measure one clock graph and predict native matter transport for both load species."""

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import factorial
import json

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply
from scipy.special import jv

from model.executable_model import cell, federation
from model.experiments.interval_echo_einstein_station import shared_variance, volume
from . import resolved_shells as shells

HERE, iv, lab, reg = shells.HERE, shells.iv, shells.dist.lab, shells.reg
bounds = shells.bounds
SEED = 28115731
TIMES = (1, 2, 4)
FIELD_NAMES = ("fresh", "shared", "cal_positive", "cal_negative")


def perfect_matchings(edges, vertices):
    if not vertices:
        yield ()
        return
    first = min(vertices)
    for a, b in edges:
        if a == first and b in vertices:
            for tail in perfect_matchings(edges, vertices - {a, b}):
                yield ((a, b), *tail)


def codebook():
    edges = sorted(cell().seams)
    common = cell().adjacency().toarray().astype(bool)
    selected = []
    for matching in perfect_matchings(edges, set(range(12))):
        proposal = common.copy()
        for a, b in matching:
            proposal[a, b] = proposal[b, a] = False
        if max(shortest_path(proposal, indices=0, unweighted=True)) <= 3:
            common = proposal
            selected.append(matching)
        if len(selected) == 10:
            break
    assert len(selected) == 10
    return selected, common


def calibrate(code):
    # A blank exchange moves an isolated q by one edge in π. A sharp
    # record gives only a phase, independent of q and its mobile endpoint.
    law = volume.candidate.native_law
    ex, co = law.local_term_parts()
    h = ex + volume.candidate.KAPPA * co
    values = (0, 32, 33)
    labels = [volume.label(0, 1), volume.label(2, 1), volume.label(2, -1)]
    changed = h.copy()
    sector = np.array([sum(v) == 2 for v in law.SRC])
    changed[:9, :9] += ex[:9, :9] * np.outer(sector, sector)
    raw, changed_arrivals = {}, {}
    for q in (1, 2):
        indices = [law.SRC.index(z) for z in ((q, 0), (0, q))]
        blank = h[np.ix_(indices, indices)]
        columns = np.eye(law.DIM)[:, indices]
        np.testing.assert_allclose(h @ columns, columns @ blank, atol=2e-13, rtol=0)
        vacuum = np.eye(law.DIM)[:, law.SRC.index((0, 0))]
        np.testing.assert_allclose(h @ vacuum, 0, atol=2e-13, rtol=0)
        np.testing.assert_allclose(
            blank, np.eye(2) * (0.5 + 8 * q * q) - 0.5 * (1 - np.eye(2)), atol=2e-13
        )
        changed_block = changed[np.ix_(indices, indices)]
        eye = np.eye(2)
        np.testing.assert_allclose(
            changed_block - 8 * q * q * eye, (blank - 8 * q * q * eye) * q, atol=2e-13
        )
        changed_arrivals[str(q)] = float(
            abs(expm(-1j * np.pi * changed_block)[1, 0]) ** 2
        )
        assert abs(abs(expm(-1j * np.pi * blank)[1, 0]) ** 2 - 1) < 1e-12
        for E, m in zip(values, labels, strict=True):
            cols = np.kron(
                m[:, None], np.eye(law.N_SRC)[:, [law.SRC.index((0, 0)), *indices]]
            )
            np.testing.assert_allclose(h @ cols, E * cols, atol=2e-13, rtol=0)
        # Archive actual sharp arrival outcomes for every internal contact.
        raw[str(q)] = [
            [int(e not in matching) for e in sorted(cell().seams)] for matching in code
        ]
    assert raw["1"] == raw["2"]
    return raw, changed_arrivals


def graph(board, calibration):
    edges = sorted(cell().seams)
    internal = [
        np.array([e for e, b in zip(edges, bits, strict=True) if b], int)
        for bits in calibration["1"]
    ]
    common = np.array(board.external)

    def active(i, j, k, w):
        mats = [
            internal[k if h == i else w if h == j else 0] + 12 * h
            for h in range(board.cells)
        ]
        a, b = np.vstack([common, *mats]).T
        u = np.r_[a, b]
        v = np.r_[b, a]
        adjacency = sparse.csr_matrix(
            (np.ones(len(u)), (u, v)), shape=(board.n, board.n)
        )
        deg = np.bincount(u, minlength=board.n)
        assert np.array_equal(deg, board.degree - 1)
        return adjacency, deg

    return active


def flights(adjacency, deg, launches, controls=False):
    # The common position generator L/2 follows from one blank hop per
    # adjacent pair and one written internal seam at every port.
    h = (sparse.diags(deg, dtype=float) - adjacency) * 0.5
    scale = float(max(deg)) / 2
    x = h / scale - sparse.eye(len(deg), format="csr")
    start = np.zeros((len(deg), len(launches)))
    start[launches, np.arange(len(launches))] = 1
    times = np.array([*TIMES, 8])
    order = 72
    z = scale * max(times)
    tail = (
        2
        * (z / 2) ** (order + 1)
        / factorial(order + 1)
        * np.exp(z * z / (4 * (order + 2)))
        / (1 - z / (2 * (order + 2)))
    )
    assert tail < 1e-25
    prev = start
    curr = x @ start
    coeff = (
        jv(np.arange(order + 1)[:, None], scale * times)
        * ((-1j) ** np.arange(order + 1))[:, None]
    )
    coeff[1:] *= 2
    Y = coeff[0, :, None, None] * prev + coeff[1, :, None, None] * curr
    for k in range(2, order + 1):
        prev, curr = curr, 2 * (x @ curr) - prev
        Y += coeff[k, :, None, None] * curr
    probs = abs(Y) ** 2
    np.testing.assert_allclose(probs.sum(axis=1), 1, atol=2e-12, rtol=0)
    if controls:
        ref = expm_multiply(-4j * h, start, traceA=-4j * h.diagonal().sum())
        np.testing.assert_allclose(abs(ref) ** 2, probs[2], atol=1e-13, rtol=0)
    return probs.reshape(4, -1, 12, len(launches)).sum(axis=2).transpose(0, 2, 1)


def native_sharp_check(board, code, sites, words, launches, predicted):
    # Independently assemble the one-load block from all physical seams
    # (including written ones), without using the surveyed adjacency.
    law = volume.candidate.native_law
    ex, cost = law.local_term_parts()
    h = ex + 8 * cost
    memories = [
        np.eye(11)[0],
        volume.label(0, 1),
        volume.label(2, 1),
        volume.label(2, -1),
    ]
    marks = {}
    for owner in range(board.cells):
        idx = words[sites.index(owner)] + 1 if owner in sites else 0
        letters = list(divmod(idx - 1, 3)) if idx else [0, 0]
        for rank, (a, b) in enumerate(code[idx]):
            marks[(12 * owner + a, 12 * owner + b)] = 1 + (
                letters[rank] if rank < 2 else 0
            )
    for q in (1, 2):
        local = []
        for m in memories:
            cols = np.kron(
                m[:, None],
                np.eye(9)[:, [law.SRC.index(v) for v in ((0, 0), (q, 0), (0, q))]],
            )
            block = cols.T @ h @ cols
            np.testing.assert_allclose(h @ cols, cols @ block, atol=2e-13, rtol=0)
            local.append(block)
        full = np.diag((8 * q * q * (6 - board.degree)).astype(float))
        for a, b in board.seams:
            block = local[marks.get((min(a, b), max(a, b)), 0)]
            full[np.ix_([a, b], [a, b])] += block[1:, 1:] - np.eye(2) * block[0, 0]
        base = np.zeros((board.n, len(launches)))
        base[launches, np.arange(len(launches))] = 1
        gen = -4j * sparse.csr_matrix(full)
        end = expm_multiply(gen, base, traceA=gen.diagonal().sum())
        observed = (abs(end) ** 2).reshape(board.cells, 12, -1).sum(axis=1).T
        np.testing.assert_allclose(observed, predicted, rtol=0, atol=5e-12)


def conditional(board, code, calibration, ball, launches, modes):
    active = graph(board, calibration)
    size = len(ball) * (len(ball) - 1)
    # Source standards are symmetric under interchanging the two occupied cells.
    # The contribution from either address order is exactly the same after mixing.
    f = np.zeros((len(modes), 4, len(launches), board.cells))
    froms, tos = np.triu_indices(board.cells, 1)
    radar = np.zeros((len(modes), len(froms), board.n))
    shape = np.zeros((board.cells, board.cells), dtype=np.int64)
    signed_law = np.rint((modes[1] - modes[0]) * 2592).astype(int)
    for n, (i, j) in enumerate(combinations(ball, 2)):
        for k in range(9):
            for w in range(9):
                weights = modes[:, k, w] * 2 / size
                adjacency, deg = active(i, j, k + 1, w + 1)
                dist = shortest_path(
                    adjacency, indices=np.arange(0, board.n, 12), unweighted=True
                )[:, ::12]
                route = dist[froms, tos].astype(int)
                radar[:, np.arange(len(froms)), route] += weights[:, None]
                shape += 2 * signed_law[k, w] * dist.astype(int)
                prob = flights(
                    adjacency, deg, launches, controls=(n == 0 and k == w == 0)
                )
                if n == 0 and (k, w) == (5, 8):
                    native_sharp_check(board, code, (i, j), (k, w), launches, prob[2])
                f += weights[:, None, None, None] * prob
    end = int(np.nonzero(radar.sum(axis=(0, 1)) > 0)[0][-1]) + 1
    return (
        f,
        radar[:, :, :end],
        shape,
        list(zip(froms.tolist(), tos.tolist(), strict=True)),
    )


def contrast(tapes, scale=1):
    # Binary/multinomial means are decoded on their original per-start scale.
    # For signed differences constant offsets in a score may be dropped.
    q = np.asarray(tapes, dtype=float) / lab.SHOTS
    bands = reg.src.probability_error(q, lab.SHOTS)
    return (q[1] - q[0]) * scale, (bands[1] + bands[0]) * scale


def read_score(difference, score):
    value, error = difference
    score = np.asarray(score, dtype=float)
    mid = value @ score
    pad = error @ abs(score) + 3e-14 * (1 + abs(score).max())
    return reg.band(mid - pad, mid + pad)


def capped(p, rng, chance):
    # The final bin records nonadmission; all conditional terminal bins remain.
    prob = p * float(chance)
    return lab.counts(np.r_[prob, 1 - prob.sum()], rng)


def report(data, score, prob, q, kappa):
    # The TV of shared-fresh for a doubly admitted word pair is 1/72.
    direct = read_score(contrast([data["fresh"], data["shared"]]), [*score, 0])
    standard = (
        read_score(contrast([data["cal_negative"], data["cal_positive"]]), score) / 72
    )
    return {
        "held": bounds(direct),
        "from_source": bounds(q * standard / kappa),
        "conditional_standard": bounds(standard),
        "population": float(prob),
    }


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    code, common = codebook()
    calibration, changed_arrivals = calibrate(code)
    board = federation(0)
    assert set(board.seams).issuperset(board.external)
    external = np.array(board.external) // 12
    owners = np.zeros((board.cells, board.cells))
    for a, b in external:
        owners[a, b] = owners[b, a] = 1
    ring = shortest_path(owners, indices=0, unweighted=True).astype(int)
    ball = np.flatnonzero(ring <= 1).tolist()
    rays = [int(i) for j in (1, 2) for i in np.flatnonzero(ring == j)[[0, -1]]]
    starts = [12 * i for i in rays]
    P = Fraction(len(ball) * (len(ball) - 1), board.cells * (board.cells - 1))
    src = volume.radar_source_and_hand(board)[-1]
    laws, _ = shared_variance.source_words(src)
    reverse, fresh, shared = np.rint(2592 * laws).astype(int)
    delta = shared - fresh
    assert np.array_equal(delta, fresh - reverse) and np.abs(delta).sum() == 72
    assert np.array_equal(
        delta.reshape(9, 3, 3).transpose(0, 2, 1).reshape(9, 9), -delta
    )
    modes = np.array(
        [
            fresh / 2592,
            shared / 2592,
            np.maximum(delta, 0) / 36,
            np.maximum(-delta, 0) / 36,
        ]
    )
    assert all(np.array_equal(z, z.T) for z in modes)
    kappa = iv.mpf(shells.dist.read("resolved_shells")["conditional_entropy_slope"])
    pos = delta > 0
    pf = Fraction(int(fresh[pos].sum()), 2592)
    ps = Fraction(int(shared[pos].sum()), 2592)
    a, x = shells.event_counts(pf, P, rng)
    b, y = shells.event_counts(ps, P, rng)
    q = (y - x) * kappa / shells.interval(ps - pf)
    # An independently reversed routing code for the second physical word.
    # Average the two placements because admissions have ordered uniform addresses.
    turned = modes[:2].reshape(2, 9, 3, 3).transpose(0, 1, 3, 2).reshape(2, 9, 9)
    turned = (turned + turned.transpose(0, 2, 1)) / 2
    modes = np.concatenate([modes, turned])
    flight, radar, shape, station_pairs = conditional(
        board, code, calibration, ball, starts, modes
    )
    for table in (flight, radar):
        np.testing.assert_allclose(
            table[5] - table[4], table[0] - table[1], atol=2e-14, rtol=0
        )
    radar_rng = np.random.default_rng(SEED + 6)
    reversal_rng = np.random.default_rng(SEED + 7)
    np.testing.assert_allclose(
        flight[1] - flight[0], (flight[2] - flight[3]) / 72, rtol=0, atol=1e-13
    )
    tables = {"radar": {}, "flight": {}}
    for name, arrays in (("radar", radar), ("flight", flight)):
        for index in np.ndindex(arrays.shape[1:-1]):
            if name == "radar":
                key = "/".join(map(str, station_pairs[index[0]]))
                steps = np.arange(arrays.shape[-1])
                metrics = {"travel_time": steps}
                tags = ("1", "2")
            else:
                t, j = index
                key = f"{[*TIMES, 8][t]}/{rays[j]}"
                metrics = {"squared_radius": ring**2, "entry": (ring <= 1).astype(int)}
                tags = ("1", "2") if t < 3 else ("changed_double",)
            row = {}
            for tag in tags:
                arms = arrays[(slice(None), *index)]
                stream = radar_rng if name == "radar" else rng
                counts = {
                    n: (capped(p, stream, P) if h < 2 else lab.counts(p, stream))
                    for h, (n, p) in enumerate(zip(FIELD_NAMES, arms[:4], strict=True))
                }
                reversal = [capped(v, reversal_rng, P) for v in arms[4:]]
                row[tag] = {
                    "counts": counts,
                    "reversed_code_counts": reversal,
                    "reversed_code_held": {
                        label: bounds(read_score(contrast(reversal), [*score, 0]))
                        for label, score in metrics.items()
                    },
                    "metrics": {
                        label: report(
                            counts,
                            score,
                            float(P) * ((arms[1] - arms[0]) @ score),
                            q,
                            kappa,
                        )
                        for label, score in metrics.items()
                    },
                }
            tables[name][key] = row
    return {
        "seed": SEED,
        "all_pairs_seed": SEED + 6,
        "reversed_code_seed": SEED + 7,
        "changed_blank_swap_arrivals": changed_arrivals,
        "starts": lab.SHOTS,
        "source_starts": shells.STARTS,
        "cells": board.cells,
        "source_ball": ball,
        "launch_cells": rays,
        "clock_radii": ring.tolist(),
        "clock_eccentricity": int(max(ring)),
        "pair_chance": str(P),
        "source_counts": {"fresh": a, "shared": b},
        "source_nats": bounds(q),
        "code": [[list(e) for e in m] for m in code],
        "internal_edges": list(map(list, sorted(cell().seams))),
        "common_clock_ranges": shortest_path(common, indices=0, unweighted=True)
        .astype(int)
        .tolist(),
        "sharp_swap_arrivals": calibration,
        "ranging_matrix_exact_numerator": shape.tolist(),
        "ranging_matrix_denominator": 2592 * board.cells * (board.cells - 1),
        "tables": tables,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in ("WORD_ROUTES.md", "resolved_shells.json")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    path = HERE / "word_routes.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    print("ball", data["source_ball"], "q", data["source_nats"])
    print(
        "clock",
        [
            (k, z["1"]["metrics"]["travel_time"]["held"])
            for k, z in data["tables"]["radar"].items()
            if abs(z["1"]["metrics"]["travel_time"]["population"]) > 1e-10
        ],
    )
    for k, z in data["tables"]["flight"].items():
        print(k, next(iter(z.values()))["metrics"])
