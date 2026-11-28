"""Read entropy, size and boundary from the same clock-defined pair regions."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

from mpmath.ctx_iv import MPIntervalContext
import networkx as nx
import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import paired_ranging as survey

HERE = survey.HERE
IV = MPIntervalContext()
IV.dps = 40
SEED = 28122862
N = 2**52


def incidence(pools, graph):
    ids = {pair: i for i, g in enumerate(pools) for pair in g}
    links = set()
    for (a, b), i in ids.items():
        for first, other in ((a, b), (b, a)):
            for u in graph[first]:
                j = ids[tuple(sorted((u, other)))]
                if i != j:
                    links.add(tuple(sorted((i, j))))
    return ids, np.array(sorted(links)).T


def sample(values, weights, scale, rng):
    # Draw clock or word readouts with all unsuccessful admissions included.
    hist = np.bincount(values.ravel(), weights=weights.ravel()).astype(np.int64)
    support = np.flatnonzero(hist)
    hist = hist[support]
    assert int(hist.sum()) == scale
    counts = rng.multinomial(N, hist / scale)
    return support, counts


def interval(q):
    x = F(q)
    return IV.mpf(x.numerator) / x.denominator


def held(values, weights, rng, entropy=False):
    total = int(weights.sum())
    support, counts = sample(values, weights, total, rng)
    score = [
        -IV.log(interval(F(int(k), 2592))) if entropy else IV.mpf(int(k))
        for k in support
    ]
    mean = sum(z * int(n) for z, n in zip(score, counts, strict=True)) / N
    # Fewer than 10000 scalar estimates in the complete experiment.
    e = (
        max(score, key=lambda z: float(z.b)) - min(score, key=lambda z: float(z.a))
    ) * IV.sqrt(IV.log(2000000) / (2 * N))
    bound = mean + IV.mpf([-1, 1]) * e
    return {
        "histogram": [[int(k), int(n)] for k, n in zip(support, counts, strict=True)]
    }, bound


def bounds(x):
    return [
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    ]


def geometry(values, weights, rng):
    # Three independent spectra: reference, positive and negative signed controls.
    reads = []
    estimate = []
    for w in weights:
        tape, band = held(values, w, rng)
        reads.append(tape)
        estimate.append(band)
    return reads, estimate[0], (estimate[1] - estimate[2]) / 72


def word_cost(gate, stock):
    nat = np.broadcast_to(stock, gate.shape)
    omitted = np.sum(nat * (~gate), axis=-1, keepdims=True)
    # Every hidden label maps to the same recorded miss. Costs depend on
    # the reference probability of that combined output given its address.
    cost = np.where(gate, nat, omitted)
    assert cost.min() > 0
    return cost


def source_spectra(gate, stock, delta, rng):
    nat = np.broadcast_to(stock, gate.shape)
    cost = word_cost(gate, stock)
    tapes = []
    bands = []
    for read in (nat, 72 * np.maximum(delta, 0), 72 * np.maximum(-delta, 0)):
        tape, band = held(cost, np.broadcast_to(read, gate.shape), rng, entropy=True)
        tapes.append(tape)
        bands.append(band)
    return tapes, bands[0], (bands[1] - bands[2]) / 72


def stationarity(src, base):
    change = {k: (v[0] + v[1]) / 2 for k, v in src.items()}
    radial = {k: v[1] - v[0] for k, v in base.items()}
    if 0 in radial["V"]:
        return {"volume_compensator_unresolved": True}
    adjust = -change["V"] / radial["V"]
    dS = change["S"] + adjust * radial["S"]
    dA = change["A"] + adjust * radial["A"]
    return {
        "coin_slope": bounds(adjust),
        "fixed_volume_entropy": bounds(dS),
        "fixed_volume_boundary": bounds(dA),
        "fixed_volume_screen_energy": bounds(
            change["screen"] + adjust * radial["screen"]
        ),
        "fixed_volume_pair_energy": bounds(change["pair"] + adjust * radial["pair"]),
        "required_coefficient": None if 0 in dA else bounds(-dS / dA),
    }


def calibration_obstruction(rows, tag, finite):
    # Eliminate all three coefficients exactly, before using the observed
    # entropy enclosures. A nonnegative witness tests maximum inequalities too.
    from itertools import combinations
    import sympy as sp

    keys = ("A", "screen", "pair")
    design = []
    entropy = []
    for row in rows:
        obs = row["readings"][tag]
        if finite:
            trial = obs["finite_dose"]
            design.append([F(trial["exact_changes"][k]) for k in keys])
            entropy.append(trial["measured_entropy_change"])
        else:
            lo, hi = obs["population"]

            def average_source(k):
                return (F(lo[k][1]) + F(hi[k][1])) / 2

            def radial(k):
                return F(hi[k][0]) - F(lo[k][0])

            coin = -average_source("V") / radial("V")
            design.append([average_source(k) + coin * radial(k) for k in keys])
            entropy.append(obs["balance"]["fixed_volume_entropy"])
    best = None
    positive = None
    for group in combinations(range(len(rows)), 4):
        matrix = sp.Matrix(
            [[sp.Rational(design[i][k]) for i in group] for k in range(3)]
        )
        for vector in matrix.nullspace():
            for sign in (1, -1):
                norm = sum(abs(x) for x in vector)
                weights = [F(sign * x / norm) for x in vector]
                assert all(
                    sum(w * design[i][k] for i, w in zip(group, weights, strict=True))
                    == 0
                    for k in range(3)
                )
                minimum = sum(
                    w * F(entropy[i][0 if w >= 0 else 1])
                    for i, w in zip(group, weights, strict=True)
                )
                if minimum <= 0:
                    continue
                candidate = (minimum, group, weights)
                if best is None or minimum > best[0]:
                    best = candidate
                if min(weights) >= 0 and (positive is None or minimum > positive[0]):
                    positive = candidate

    def packed(z):
        if z is None:
            return None
        minimum, group, weights = z
        return {
            "center_weights": [
                [rows[i]["center"], str(w)]
                for i, w in zip(group, weights, strict=True)
                if w
            ],
            "certified_minimum_max_error_nats": str(minimum),
        }

    return {"equality": packed(best), "one_sided": packed(positive)}


def run():
    board = survey.bulk.federation(0)
    graphs, points = survey.bulk.geometry._atlas([board])
    graph = graphs[0]
    hop = nx.floyd_warshall_numpy(graph, nodelist=range(board.cells)).astype(int)
    pools, means = survey.midpoint_classes(points[0], hop)
    ids, links = incidence(pools, graph)
    _, _, plans = survey.shells(means)
    centers = [plans[1][0][0], *[t[1] for t in plans[1][:6]], ids[(0, 0)]]
    query = survey.queries(pools, [(c, i) for c in centers for i in range(len(pools))])
    active = survey.routes.graph(
        board, survey.routes.calibrate(survey.routes.codebook()[0])[0]
    )

    def distances(i, j, k, w):
        adj = active(i, j, k, w)[0]
        d = shortest_path(adj, indices=np.arange(0, board.n, 12), unweighted=True)[
            :, ::12
        ].astype(int)
        return survey.deliver(d, query).reshape(2, len(centers), len(pools))

    absence = distances(-1, -1, 0, 0)
    deadlines = (
        np.quantile(absence, [1 / 3, 2 / 3], axis=-1, method="inverted_cdf")
        .astype(int)
        .transpose(1, 2, 0)
    )
    spots = [pair for pair in ids if pair[0] != pair[1]]
    shape = (len(spots), 81, 2, len(centers), 2)
    sizes = np.zeros(shape, dtype=np.int16)
    areas = np.zeros_like(sizes)
    gates = np.zeros(shape, bool)
    boosts = {k: np.zeros(shape, dtype=np.int32) for k in ("screen", "pair")}
    energies = np.asarray(
        survey.routes.volume.radar_source_and_hand(board)[-1]["source_word_energies"],
        int,
    )
    word_energy = (energies[:, None] + energies[None, :]).ravel()
    for si, (i, j) in enumerate(spots):
        for k, w in np.ndindex(9, 9):
            d = distances(i, j, k + 1, w + 1)
            inside = d[:, :, None, :] <= deadlines[..., None]
            at = (si, 9 * k + w)
            sizes[at] = inside.sum(axis=-1)
            areas[at] = (inside[..., links[0]] != inside[..., links[1]]).sum(axis=-1)
            gates[at] = inside[..., ids[i, j]]

            def kernel(g):
                return np.maximum(deadlines**2 - d[..., g, None] ** 2, 0)

            boosts["screen"][at] = word_energy[k] * kernel(ids[i, i]) + word_energy[
                w
            ] * kernel(ids[j, j])
            boosts["pair"][at] = (word_energy[k] + word_energy[w]) * kernel(ids[i, j])
    stocks, turned = survey.standards(board)
    rng = np.random.default_rng(SEED)
    out = {}
    for mode, choice in enumerate(("fixed", "best")):
        out[choice] = []
        for ci, c in enumerate(centers):
            alts = {}
            for tag, ref, sign in (("native", 0, 1), ("reversed_second", 5, -1)):
                base = {k: [] for k in ("S", "A", "V", "screen", "pair")}
                src = {k: [] for k in base}
                finite_entropy = []
                population = []
                tapes = []
                for di in range(2):
                    ix = (slice(None), slice(None), mode, ci, di)
                    geom_weights = [
                        np.broadcast_to(stocks[k].ravel(), sizes[ix].shape)
                        for k in (ref, 3, 4)
                    ]
                    panel = {}
                    exact = {}
                    for key, array in (("V", sizes), ("A", areas), *boosts.items()):
                        tape, avg, diff = geometry(array[ix], geom_weights, rng)
                        factor = (
                            int(2 * deadlines[mode, ci, di]) if key in boosts else 1
                        )
                        avg /= factor
                        diff /= factor
                        scale = int(geom_weights[0].sum()) * factor
                        exact[key] = [
                            F(int((array[ix] * geom_weights[0]).sum()), scale),
                            sign
                            * F(
                                int(
                                    (
                                        array[ix] * (geom_weights[1] - geom_weights[2])
                                    ).sum()
                                ),
                                72 * scale,
                            ),
                        ]
                        panel[key] = tape
                        base[key].append(avg)
                        src[key].append(sign * diff)
                    population.append(exact)
                    order = (
                        [(stocks[0], stocks[1])]
                        if tag == "native"
                        else [(turned[0], turned[1]), (turned[0].T, turned[1].T)]
                    )
                    panel["S"] = []
                    avg = diff = h_end = 0
                    panel["endpoint_entropy"] = []
                    for referee, endpoint in order:
                        tape, v, d = source_spectra(
                            gates[ix],
                            referee.ravel(),
                            sign * (stocks[2] - stocks[1]).ravel(),
                            rng,
                        )
                        panel["S"].append(tape)
                        end_tape, he = held(
                            word_cost(gates[ix], endpoint.ravel()),
                            np.broadcast_to(endpoint.ravel(), gates[ix].shape),
                            rng,
                            entropy=True,
                        )
                        panel["endpoint_entropy"].append(end_tape)
                        h_end += he / len(order)
                        avg += v / len(order)
                        diff += d / len(order)
                    finite_entropy.append(h_end)
                    tapes.append(panel)
                    base["S"].append(avg)
                    src["S"].append(diff)

                def matched(key, coin, t):
                    return sum(
                        frac * (population[ii][key][0] + t * population[ii][key][1])
                        for ii, frac in enumerate((1 - coin, coin))
                    )

                half = F(1, 2)
                v0 = matched("V", half, 0)
                vl, vh = [matched("V", a, 1) for a in (0, 1)]
                coin = (v0 - vl) / (vh - vl)
                assert 0 < coin < 1 and matched("V", coin, 1) == v0
                alpha = interval(coin)
                coin_entropy = -alpha * IV.log(alpha) - (1 - alpha) * IV.log(1 - alpha)
                ds = (
                    sum(
                        h * w
                        for h, w in zip(finite_entropy, [1 - alpha, alpha], strict=True)
                    )
                    - sum(base["S"]) / 2
                    + coin_entropy
                    - IV.log(2)
                )
                finite = {
                    "upper_deadline_fraction": str(coin),
                    "matched_volume": str(v0),
                    "measured_entropy_change": bounds(ds),
                    "exact_changes": {
                        k: str(matched(k, coin, 1) - matched(k, half, 0)) for k in exact
                    },
                }
                alts[tag] = {
                    "finite_dose": finite,
                    "population": [
                        {k: list(map(str, v)) for k, v in d.items()} for d in population
                    ],
                    "balance": stationarity(src, base),
                    "source": {k: list(map(bounds, v)) for k, v in src.items()},
                    "reference": {k: list(map(bounds, v)) for k, v in base.items()},
                    "tapes": tapes,
                }
            out[choice].append(
                {
                    "center": c,
                    "depth": float(np.linalg.norm(means[c])),
                    "deadlines": deadlines[mode, ci].tolist(),
                    "home_admitted_without_words": (
                        absence[mode, ci, ids[0, 0]] <= deadlines[mode, ci]
                    ).tolist(),
                    "readings": alts,
                }
            )
    return {
        "seed": SEED,
        "shots": N,
        "markers": len(pools),
        "boundary_links": links.shape[1],
        "metric_results": out,
        "calibration_obstructions": {
            k: {
                tag: {
                    test: calibration_obstruction(rows, tag, finite)
                    for test, finite in (("slope", False), ("finite", True))
                }
                for tag in ("native", "reversed_second")
            }
            for k, rows in out.items()
        },
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("PAIRED_BALLS.md", "paired_ranging.json", "resolved_shells.json")
        },
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    z = run()
    path = HERE / "paired_balls.json"
    if args.write:
        path.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    for k, rows in z["metric_results"].items():
        print(
            k,
            [
                (
                    r["center"],
                    r["readings"]["native"]["balance"].get("required_coefficient"),
                )
                for r in rows
            ],
        )
