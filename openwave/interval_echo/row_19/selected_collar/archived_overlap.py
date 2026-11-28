"""Measure the entropy of locally archived write histories and ordinary material."""

from fractions import Fraction as F
from itertools import combinations
from math import log, sqrt

from mpmath import iv
import numpy as np
import sympy as sp

from . import coherent_balls as coherent

steady = coherent.steady
local = coherent.local
bounds = coherent.balls.bounds
FLAGS = steady.FLAGS
SEED = 28122920
N = coherent.N


def probabilities(counts, err):
    return [iv.mpf([max(0, k / N - err), min(1, k / N + err)]) for k in counts]


def entropy(p):
    return -sum(z * iv.log(z) for z in p if z != 0)


def local_edges(sites, edges):
    index = {v: i for i, v in enumerate(sites)}
    return [(index[a], index[b]) for a, b in edges if a in index and b in index]


def quiet_span(h):
    operator = sp.Matrix(h.tolist())
    columns = [sp.ones(len(h), 1)]
    while True:
        nxt = operator * columns[-1]
        basis = sp.Matrix.hstack(*columns)
        if basis.row_join(nxt).rank() == len(columns):
            return basis
        columns.append(nxt)


def null_read(sites, deposit, span, edges, blank, quiet_weight, rng):
    # Isolate the component of a localized ordinary load orthogonal to
    # every possible quiet history, including every finite waiting time.
    k = sp.Matrix.hstack(*span[sites, :].columnspace())
    e = sp.eye(len(sites))[:, sites.index(deposit)]
    v = e - k * (k.T * k).inv() * k.T * e
    length = F((v.T * v)[0])
    if length == 0:
        return None
    assert v.T * k == sp.zeros(1, k.cols)
    unit = np.array(v, float).ravel() / sqrt(float(length))
    _, vecs = np.linalg.eigh(np.outer(unit, unit))
    rotation, circuit = coherent.trap.local_rotation(vecs.T, local_edges(sites, edges))
    meter = rotation[-1]
    np.testing.assert_allclose(abs(meter @ unit), 1, atol=5e-13)
    assert abs(meter @ blank[np.ix_(sites, sites)] @ meter) < 1e-12
    assert abs(meter[sites.index(deposit)] ** 2 - float(length)) < 1e-12
    w = float(quiet_weight)
    reads = {}
    for dose in (F(0), F(1, 128), F(1, 32), F(1)):
        hit = float(dose * quiet_weight * length)
        reads[str(dose)] = coherent.draws([hit, w - hit, 1 - w], rng)
    return {
        "sites": sites,
        "quiet_rank_in_ball": k.cols,
        "exact_null_vector": [str(x) for x in v],
        "ordinary_hit_given_quiet": str(length),
        "ordinary_joint_hit": str(length * quiet_weight),
        "readout": circuit,
        "hit_miss_other_counts": reads,
    }


def conditional_read(sites, edges, layouts, rng, err):
    # The far inner ball is a single port and its source read is classical.
    # The surrounding six-port ball reaches both history registers.
    per_flag, conditional_entropy, vac = {}, [], []
    for f in FLAGS:
        m = layouts[f][np.ix_(sites, sites)]
        _, basis = np.linalg.eigh(m)
        u, circuit = coherent.trap.local_rotation(basis.T, local_edges(sites, edges))
        x = u @ m @ u.T
        count = coherent.draws([*x.diagonal(), 1 - np.trace(x)], rng)
        sample = [iv.mpf(k) / N for k in count]
        pair_counts = []
        rows = np.zeros(len(sites))
        for i, j in combinations(range(len(sites)), 2):
            avg = (x[i, i] + x[j, j]) / 2
            pair = [
                coherent.draws([avg + y, avg - y, 1 - 2 * avg], rng)
                for y in (x[i, j].real, x[i, j].imag)
            ]
            off = sum(abs(k[0] - k[1]) / (2 * N) + err for k in pair)
            rows[[i, j]] += off
            pair_counts.append(pair)
        # The reported diagonal histogram and the true conditional matrix
        # have trace distance <= d/2 times the max absolute row error.
        d = len(count)
        eps = iv.mpf(d) * (max(rows) + err + 1e-14) / 2
        assert float(eps.b) < 1 - 1 / d
        radius = entropy([eps, 1 - eps]) + eps * iv.log(d - 1)
        ent = entropy(sample) + iv.mpf([-float(radius.b), float(radius.b)])
        conditional_entropy.append(ent)
        vac.append(probabilities(count, err)[-1])
        per_flag[str(f)] = {
            "mode_counts": count,
            "phase_counts": pair_counts,
            "circuit": circuit,
            "conditional_entropy": bounds(ent),
            "trace_distance_upper": float(eps.b),
        }
    return per_flag, conditional_entropy, vac


def finite_pair(row, setup, layouts, weights, calibrated, rng, err):
    near, far = setup[3]
    inner, outer = row["ball_sites"]
    assert inner == [far] and near in outer and setup[4][0] not in outer
    reads = {}
    single = {}
    for tag in ("fresh", "shared", "ordinary"):
        ws = weights["fresh" if tag == "ordinary" else tag]
        cells = np.zeros((2, 2))
        for f, w in zip(FLAGS, ws, strict=True):
            seen = 0 if tag == "ordinary" else layouts[f][far, far]
            cells[f[1]] += float(w) * np.array([seen, 1 - seen])
        reads[tag] = coherent.draws(cells.ravel(), rng)
        single[tag] = probabilities(reads[tag], err)
    source = single["fresh"]
    inner_s = entropy(source)
    slopes = {
        t: -sum((x - y) * iv.log(y) for x, y in zip(single[t], source, strict=True))
        for t in ("shared", "ordinary")
    }
    flag_counts = {
        k: coherent.draws(list(map(float, v)), rng) for k, v in weights.items()
    }
    fresh, shared = [probabilities(flag_counts[t], err) for t in ("fresh", "shared")]
    _, ent, vacuum = calibrated
    outer_s = entropy(fresh) + sum(w * h for w, h in zip(fresh, ent, strict=True))
    outer_ds = sum(
        (x - y) * (h - iv.log(y)) for x, y, h in zip(shared, fresh, ent, strict=True)
    )
    outer_os = sum(
        w * (-iv.log(v) - h) for w, v, h in zip(fresh, vacuum, ent, strict=True)
    )
    gains = {}
    for tag, outer_slope, area in (
        ("shared", outer_ds, "record_source"),
        ("ordinary", outer_os, "ordinary_standard"),
    ):
        a = iv.mpf(row[area]["area_tangent"])
        coin = iv.mpf(row[area]["public_probability_tangent"])
        slope = (slopes[tag] + outer_slope) / 2 + coin * (outer_s - inner_s)
        gains[tag] = (a, slope)
    a, s = gains["shared"]
    ao, so = gains["ordinary"]
    return {
        "inner_joint_counts": reads,
        "flag_counts": flag_counts,
        "overlap_gain": bounds(-s / a),
        "ordinary_gain": bounds(-so / ao),
        "determinant": bounds(a * so - ao * s),
        "additional_entropy_target": bounds(a * so / ao - s),
    }


def run():
    iv.dps = 45
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    setup, _ = steady.symmetric_instrument(board)
    layouts, h = coherent.stationary_layouts(setup)
    span = quiet_span(h)
    saved = local.inputs("restored_balls")["surveys"]
    wait = F(
        local.inputs("stationary_polarization")["mirror_matched"][
            "combined_trace_error_upper"
        ]
    )
    err = sqrt(log(2 * 100000 / 0.01) / (2 * N)) + float(wait) + 5e-13
    forward, reverse, _ = steady.parity(board, setup[1])
    deposit = setup[4][0]
    far_outer = next(
        z["ball_sites"][1]
        for z in saved
        if z["center"] == setup[2][setup[3][1]] and z["hop_radii"] == [0, 1]
    )
    calibration = conditional_read(far_outer, setup[6], layouts, rng, err)
    surveys = []
    distinct = {}
    for name, law in (("forward", forward), ("reversed", reverse)):
        weights = {
            name: [(1 - p) * int(f == (0, 0)) + p * law[f] for f in FLAGS]
            for name, p in (("fresh", F(25, 36)), ("shared", F(101, 144)))
        }
        for row in saved:
            if row["direction"] != name:
                continue
            keys = []
            for sites in row["ball_sites"]:
                if deposit not in sites or not set(setup[3]) <= set(sites):
                    continue
                key = name + ":" + str(sites)
                if key not in distinct:
                    distinct[key] = null_read(
                        sites,
                        deposit,
                        span,
                        setup[6],
                        layouts[0, 0],
                        weights["fresh"][0],
                        rng,
                    )
                if distinct[key]:
                    keys.append(key)
            out = {
                "direction": name,
                "center": row["center"],
                "radii": row["hop_radii"],
            }
            if keys:
                out["null_measurements"] = keys
                # The public ball selector is fair at the reference and its
                # bounded derivative cannot change the leading logarithm.
                out["ordinary_log_coefficient_lower"] = str(
                    sum(F(distinct[k]["ordinary_joint_hit"]) for k in keys) / 2
                )
            else:
                out["finite"] = finite_pair(
                    row, setup, layouts, weights, calibration, rng, err
                )
            surveys.append(out)
    return {
        "seed": SEED,
        "starts": N,
        "quiet_global_rank": span.cols,
        "sites": setup[2],
        "far_outer_conditionals": calibration[0],
        "null_reads": distinct,
        "surveys": surveys,
        "inputs": steady.exp.save_inputs(
            "ARCHIVED_OVERLAP.md",
            "restored_balls.json",
            "coherent_balls.json",
            "stationary_polarization.json",
        ),
    }


if __name__ == "__main__":
    result = run()
    steady.exp.emit("archived_overlap", result)
    for row in result["surveys"]:
        print(
            row["direction"],
            row["center"],
            row["radii"],
            row.get("ordinary_log_coefficient_lower")
            or (row["finite"]["overlap_gain"], row["finite"]["ordinary_gain"]),
        )
