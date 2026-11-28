"""Test regional entropy stationarity on counts alone and enumerate coarse time records."""

import argparse
import json
from functools import lru_cache
from hashlib import sha256
from itertools import combinations

import numpy as np
from mpmath import iv

from . import returned_head_lab as lab
from . import returned_sources as src
from . import six_station_read as six

SCORE = np.arange(2, 13) / 4
FINE = 10
MASK = ((np.arange(1 << FINE)[:, None] >> np.arange(FINE)) & 1).astype(float)


@lru_cache(maxsize=1)
def partitions():
    """Unlabeled, disjoint bins on ten successful ticks, with zero padding."""
    parts = [[]]
    for bit in range(FINE):
        v = 1 << bit
        parts = [
            row[:i] + [row[i] | v] + row[i + 1 :]
            for row in parts
            for i in range(len(row))
        ] + [row + [v] for row in parts]
    padded = np.array([x + [0] * (FINE - len(x)) for x in parts], dtype=np.int16)
    assert len(parts) == 115975 and np.all(padded.sum(1) == (1 << FINE) - 1)
    return padded


def band(a, b):
    return iv.mpf([float(a), float(b)])


def timing_slopes(tape, n, contrasts=src.CONTRASTS):
    means = {}
    for key, rows in tape.items():
        prob = np.array(rows, dtype=float) / n
        # Invert Bernstein for a capped mean. Bound its variance by the
        # second moment of the deficit below timeout using the same event.
        high = prob + src.probability_error(prob, n)
        var = high @ (3 - SCORE) ** 2
        a = np.log(2 * src.FAMILY / 0.01) / n
        e = np.sqrt(2 * a * var) + 5 * a / 3 + 1e-13
        mid = 3 + prob @ (SCORE - 3)
        means[key] = [band(m - u, m + u) for m, u in zip(mid, e, strict=True)]
    slopes = np.column_stack(
        [
            [s * (x - y) for x, y in zip(means[hi], means[lo], strict=True)]
            for hi, lo, s in contrasts
        ]
    )
    return means, slopes


def geometry(tape, n, contrasts=src.CONTRASTS):
    means, slopes = timing_slopes(tape, n, contrasts)
    return six.geometry_slopes(means["reference"], slopes)


def cross(a, v):
    return [
        a[1] * v[2] - a[2] * v[1],
        a[2] * v[0] - a[0] * v[2],
        a[0] * v[1] - a[1] * v[0],
    ]


def interval_product(a, b):
    return np.stack(
        (
            a[..., 0] * b[..., 0],
            a[..., 1] * b[..., 0],
            a[..., 0] * b[..., 1],
            a[..., 1] * b[..., 1],
        ),
        axis=-1,
    )


def score_balance(q, weight):
    return sum(w * band(*r) for w, r in zip(weight, q, strict=True))


def separated(limits):
    return bool(limits[0] > 0 or limits[1] < 0)


def subset_probabilities(tape, n):
    out = {}
    for key, rows in tape.items():
        raw = np.array(rows)
        assert raw.dtype.kind in "iu" and np.all(raw >= 0) and np.all(raw.sum(1) == n)
        prob = raw / n
        err = src.probability_error(prob, n)
        # All unions share the single event enclosing the original 11 bins.
        mid = prob[:, :-1] @ MASK.T
        width = err[:, :-1] @ MASK.T
        # Complements may be tighter. Subtraction in the tail has 1e-15 slack.
        width = np.minimum(width, err.sum(1)[:, None] - width + 1e-15)
        mid = np.column_stack((mid, prob[:, -1]))
        width = np.column_stack((width, err[:, -1]))
        out[key] = np.stack(
            (np.maximum(0, mid - width), np.minimum(1, mid + width)), axis=-1
        )
    return out


def time_scores(probs, indices):
    p = probs["reference"][indices, 1:, :]
    assert np.min(p) > 0
    pooled = p.mean(axis=0)
    prices = {
        "connected": np.stack(
            (np.log(pooled[..., 0] / p[..., 1]), np.log(pooled[..., 1] / p[..., 0])),
            axis=-1,
        ),
        "conditional_entropy": -np.log(p[..., ::-1]),
    }
    delta = np.array(
        [
            scale * (probs[hi][indices, 1:, :] - probs[lo][indices, 1:, ::-1])
            for hi, lo, scale in src.CONTRASTS
        ]
    )
    results = {}
    for label, price in prices.items():
        z = interval_product(delta, price)
        # For each possible group record its independent contribution in
        # this region. Mask zero contributes exactly zero. Timeout is last.
        vals = np.stack((z.min(-1).mean(1), z.max(-1).mean(1)), axis=-1)
        results[label] = np.pad(vals, ((0, 0), (1, 0), (0, 0)))
    return results


def evaluate_partitions(table, weights):
    # All partitions and their publicly labeled lotteries share the same
    # geometric cofactors. Multiply and sum only after fixing the contrast.
    W = np.array([six.bounds(w) for w in weights])
    prod = interval_product(table, W[:, None, :])
    slots = np.stack((prod.min(-1).sum(0), prod.max(-1).sum(0)), axis=-1)
    ids = partitions()
    vals = slots[ids].sum(1) + slots[-1]
    # Cover binary64 accumulation and transcendentals; the timer bounds
    # are far wider than the roundoff at the absolute scales involved.
    pad = 3e-14 * (np.max(np.abs(W)) + np.max(abs(vals)) + 1e-25)
    vals += np.array([-pad, pad])
    lo, hi = vals[:, 0], vals[:, 1]

    def qsum(bins):
        return table[:, [*bins, -1], :].sum(1)

    return {
        "total": len(ids),
        "excluded": int(np.count_nonzero((lo > 0) | (hi < 0))),
        "all_lotteries": [float(lo.min()), float(hi.max())],
        "full_time": qsum(1 << np.arange(FINE)).tolist(),
        "echo_only": qsum([2**FINE - 1]).tolist(),
    }, vals


def measure_one(q, A, V, w):
    Q = [band(*p) for p in q]
    res = {
        "charge": six.bounds(Q[-1]),
        "stationarity_residual": six.bounds(score_balance(q, w)),
    }
    if not w[2].a <= 0 <= w[2].b:
        res["ordinary_area_coefficient"] = six.bounds(
            (-Q[0] * V[1] + Q[1] * V[0]) / w[2]
        )
        res["held_prediction"] = six.bounds(-(w[0] * Q[0] + w[1] * Q[1]) / w[2])
    return res


def rank_witness(matrix):
    # A square minor bounded away from zero excludes every nonzero common
    # linear combination, even allowing an arbitrary sign on each component.
    candidates = []
    d = len(matrix[0])
    for ids in combinations(range(len(matrix)), d):
        center = np.array([[float(matrix[i][j].mid) for j in range(d)] for i in ids])
        scales = np.linalg.norm(center, axis=0)
        if min(scales) > 0:
            candidates.append((abs(np.linalg.det(center / scales)), ids))
    for _, indices in sorted(candidates, reverse=True):
        minor = np.array([matrix[i] for i in indices], dtype=object)
        det = six.survey.geometry.curved.det_interval(minor)
        if not det.a <= 0 <= det.b:
            return {"quartet_indices": list(indices), "determinant": six.bounds(det)}
    return None


def run(record):
    iv.dps = 35
    n = record["shots_per_setting"]
    tapes = record["tapes"]
    inc, returned, unresolved = src.responses(tapes, n)
    timing = {k: t["time"] for k, t in tapes.items()}
    probs = subset_probabilities(timing, n)
    summary = []
    rank = []
    family_values = {}
    time_exclusion = {}
    for sites, A, V in geometry(timing, n):
        ids = [i for i, r in enumerate(six.survey.ROUTES) if all(x in sites for x in r)]
        w = cross(A, V)
        # Normalize by a positive known scale for a better conditioned
        # numerical table. The determinant-zero condition is homogeneous.
        factor = max(float(abs(x).b) for x in w)
        assert factor > 0
        w = [x / factor for x in w]
        reg = returned[ids].mean(0)
        menus = {
            "incoming_heads": inc,
            "returned_heads": reg,
            "incoming_less_returned": inc - reg[:, ::-1],
        }
        row = {
            "stations": list(sites),
            "area_slopes": [six.bounds(x) for x in A],
            "volume_slopes": [six.bounds(x) for x in V],
            "cofactor_scale": factor,
            "time_partitions": {},
        }
        for label, table in time_scores(probs, ids).items():
            stats, vals = evaluate_partitions(table, w)
            # Intersect surviving deterministic choices across ALL regions.
            ok = (vals[:, 0] <= 0) & (vals[:, 1] >= 0)
            if label not in time_exclusion:
                time_exclusion[label] = np.ones(len(ok), bool)
            time_exclusion[label] &= ok
            family_values.setdefault(label, []).append(vals)
            for mode in ("full_time", "echo_only"):
                menus[label + "_" + mode] = np.array(stats.pop(mode))
            row["time_partitions"][label] = stats
        row["meters"] = {
            label: measure_one(q, A, V, [x * factor for x in w])
            for label, q in menus.items()
        }
        labels = [
            "incoming_heads",
            "returned_heads",
            "connected_full_time",
            "conditional_entropy_full_time",
        ]
        rank.append([score_balance(menus[label], w) for label in labels])
        summary.append(row)
    return {
        "input_sha256": sha256(
            (lab.HERE / "returned_head_lab.json").read_bytes()
        ).hexdigest(),
        "quartets": summary,
        "unresolved_returned_reference_blocks_out_of_60": unresolved,
        "surviving_time_conventions": {
            k: int(v.sum()) for k, v in time_exclusion.items()
        },
        "excluded_all_public_lotteries_in_quartets": {
            k: [
                i
                for i, v in enumerate(vals)
                if separated([v[:, 0].min(), v[:, 1].max()])
            ]
            for k, vals in family_values.items()
        },
        "joint_meter_order": labels,
        "nonzero_joint_minor": rank_witness(rank),
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    srcfile = lab.HERE / "returned_head_lab.json"
    dest = lab.HERE / "regional_balance.json"
    data = run(json.loads(srcfile.read_bytes()))
    if cli.parse_args().write:
        dest.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dest.read_bytes()), data)
    print(
        data["surviving_time_conventions"],
        data["excluded_all_public_lotteries_in_quartets"],
        data["nonzero_joint_minor"],
    )
