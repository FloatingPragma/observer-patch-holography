"""Measure the modular response lost when two different references cross the same channel."""

import argparse
import json
from hashlib import sha256
from itertools import combinations

import numpy as np
from scipy.linalg import eigvalsh

from . import regional_balance as reg
from . import returned_sources as src

HERE = reg.lab.HERE
CONTRASTS = (*src.CONTRASTS, ("product", "reference", 1.0))


def fitted_order(X, e, P, f):
    I = np.eye(4)
    lower = eigvalsh(P - f * I, X + e * I)[0]
    upper = eigvalsh(P + f * I, X - e * I)[-1]
    assert lower > 0
    return float(np.log(lower) - 1e-12), float(np.log(upper) + 1e-12)


def contrast_blocks(stock):
    for hi, lo, s in CONTRASTS:
        Hp, ep = stock[hi]
        Hm, em = stock[lo]
        yield s * (Hp - Hm), s * (ep + em)
    yield stock["reference"]  # Absolute deficit uses minus the X expectation.


def instrument(tapes, n, returned, limits=None):
    obs = {name: src.observed(t, n, returned) for name, t in tapes.items()}
    x, ex = obs["reference"]
    p, ep = obs["product"]
    pairs = list(contrast_blocks(obs))
    estimate = np.zeros((5, 2))
    local = np.zeros((5, 2))
    true_limits = []
    failures = 0
    for b in range(2):
        X, P = x[b], p[b]
        exb, epb = 2 * ex[b], 2 * ep[b]
        ev, U = np.linalg.eigh(X)
        fv, V = np.linalg.eigh(P)
        if not returned:
            order = fitted_order(X, exb, P, epb)
            true_limits.append(order)
        else:
            order = limits[b]
        good = ev[0] > exb and fv[0] > epb
        if good:
            lp = (V * np.log(fv)) @ V.conj().T
            K = lp - (U * np.log(ev)) @ U.conj().T
            kp = epb / (fv[0] - epb) + 1e-13
            err = exb / (ev[0] - exb) + kp + 1e-13

            def coefficients(op):
                return (
                    abs(np.einsum("ijkl,lk->ij", reg.lab.heads.PRODUCTS, op)).sum() / 4
                )

            coeff = coefficients(K)
            base_coeff = coefficients(lp)
        else:
            failures += 1
        for k, (D, e) in enumerate(pairs):
            norm = np.linalg.norm(D[b], ord="nuc") + 4 * e[b]
            limit = norm * max(map(abs, order)) + 1e-15
            lo, hi = -limit, limit
            if good:
                val = np.trace(D[b] @ K).real
                width = coeff * e[b] + err * norm + 1e-15
                lo, hi = max(lo, val - width), min(hi, val + width)
                if not returned:
                    val = -np.trace(D[b] @ lp).real
                    width = base_coeff * e[b] + kp * norm + 1e-15
                    local[k] += [val - width, val + width]
            assert lo <= hi
            estimate[k] += [lo, hi]
    if returned:
        # Timeout is a single orthogonal recorded flag. Read its frequency
        # on the same Pauli appointment in each arm, retaining all starts.
        t = {k: v[0][0][-1] / n for k, v in tapes.items()}
        timeout = {
            k: reg.band(
                f - src.probability_error(f, n), f + src.probability_error(f, n)
            )
            for k, f in t.items()
        }
        log = reg.iv.ln(timeout["product"] / timeout["reference"])
        moments = [s * (timeout[hi] - timeout[lo]) for hi, lo, s in CONTRASTS]
        for i, v in enumerate([*moments, timeout["reference"]]):
            estimate[i] += reg.six.bounds(v * log)
    return estimate, local, true_limits, failures, pairs


def classical(probs, ids):
    x = probs["reference"][ids, 1:, :]
    p = probs["product"][ids, 1:, :]
    assert np.min(p) > 0 and np.min(x) > 0
    price = np.stack(
        (np.log(p[..., 0] / x[..., 1]), np.log(p[..., 1] / x[..., 0])), axis=-1
    )
    d = np.array(
        [
            s * (probs[hi][ids, 1:, :] - probs[lo][ids, 1:, ::-1])
            for hi, lo, s in CONTRASTS
        ]
    )
    products = reg.interval_product(d, price)
    table = np.stack((products.min(-1).mean(1), products.max(-1).mean(1)), axis=-1)
    return np.pad(table, ((0, 0), (1, 0), (0, 0)))


def necessary_gain(A, V, full, orders, pairs):
    """Enclose admissible eta for ANY regional quantum channel on the same two inputs."""
    lower, upper = 0.0, np.inf
    rows = []
    # Route, parity and public random flags can all be included in one
    # CPTP map. Use the widest order interval across the two input blocks.
    radius = (max(o[1] for o in orders) - min(o[0] for o in orders)) / 2
    for i, j in combinations(range(4), 2):
        # W = V_j D_i - V_i D_j has zero volume and zero total trace.
        # The identity part of the pulled-back logarithm drops out.
        vs = [V[j], -V[i]]
        d = np.array([pairs[k][0] for k in (i, j)])
        e = np.array([pairs[k][1] for k in (i, j)])
        mid = np.array([float(v.mid) for v in vs])
        widths = np.array([float((v - v.mid).__abs__().b) for v in vs])
        norm = 0.0
        for b in range(2):
            s = np.array([np.linalg.norm(block, ord="nuc") for block in d[:, b]])
            norm += (
                np.linalg.norm(np.tensordot(mid, d[:, b], axes=1), ord="nuc")
                + widths @ s
                + (abs(mid) + widths) @ (4 * e[:, b])
            )
        c = V[j] * A[i] - V[i] * A[j]
        if c.a <= 0 <= c.b:
            continue
        entropy = V[j] * reg.band(*full[i]) - V[i] * reg.band(*full[j])
        limit = (norm + 1e-15) * radius + 1e-16
        interval = reg.six.bounds((reg.band(-limit, limit) - entropy) / c)
        lower = max(lower, interval[0])
        upper = min(upper, interval[1])
        rows.append({"directions": [i, j], "gain": interval})
    return {
        "lower": lower,
        "upper": upper if np.isfinite(upper) else None,
        "pair_tests": rows,
    }


def run():
    reg.iv.dps = 35
    paths = [HERE / "returned_head_lab.json", HERE / "product_pulse.json"]
    original, product = [json.loads(p.read_bytes()) for p in paths]
    assert (
        original["shots_per_setting"] == product["shots_per_setting"]
        and original["ports"] == product["ports"]
    )
    assert product["source_bank_sha256"] == sha256(paths[0].read_bytes()).hexdigest()
    n = product["shots_per_setting"]
    tapes = {**original["tapes"], "product": product["tape"]}
    incoming, base, orders, _, pairs = instrument(
        {k: r["incoming"] for k, r in tapes.items()}, n, False
    )
    # The manufacturer prepares P by sharing parity alone. Check the two
    # individual marginals against the independently measured reference.
    ref, ex = src.observed(tapes["reference"]["incoming"], n, False)
    prod, ep = src.observed(tapes["product"]["incoming"], n, False)
    # Use the four measured one-head and parity coefficient families.
    coords = np.einsum("ijkl,blk->bij", reg.lab.heads.PRODUCTS, ref - prod).real
    assert np.all(abs(coords[:, :, 0]) <= (ex + ep)[:, None]) and np.all(
        abs(coords[:, 0, :]) <= (ex + ep)[:, None]
    )
    corrected = incoming + base
    # The incoming one-head/parity price has exactly zero source, rotation,
    # and erasure slope on the supplied matched preparations.
    assert np.all((base[1:4, 0] <= 0) & (base[1:4, 1] >= 0))
    base[1:4] = 0
    all_out = []
    unresolved = 0
    for i in range(30):
        r, _, _, bad, _ = instrument(
            {k: t["returned"][i] for k, t in tapes.items()}, n, True, orders
        )
        all_out.append(r)
        unresolved += bad
    all_out = np.array(all_out)
    times = {k: r["time"] for k, r in tapes.items()}
    probs = reg.subset_probabilities(times, n)
    quartets = []
    minors = {"pair": [], "with_ordinary": []}
    admitted = np.ones(len(reg.partitions()), bool)
    for sites, A, V in reg.geometry(times, n, CONTRASTS):
        ids = [
            i for i, r in enumerate(reg.six.survey.ROUTES) if all(x in sites for x in r)
        ]
        transmitted = all_out[ids].mean(0)
        lost = incoming - transmitted[:, ::-1]
        meters = {
            "transmitted": transmitted,
            "lost": lost,
            "ordinary_plus_lost": base + lost,
        }
        q = classical(probs, ids)
        # Put the incoming and ordinary scores into the mandatory timeout
        # slot as a constant. Every time partition includes this term once.
        total_table = -q[:, :, ::-1].copy()
        total_table[:, -1] += incoming[:4] + base[:4]
        result = {
            "stations": list(sites),
            "absolute_data_processing_loss": (-lost[4, ::-1]).tolist(),
            "symmetric_data_processing_loss": lost[3].tolist(),
            "erasure_area_slope": reg.six.bounds(A[3]),
            "erasure_volume_slope": reg.six.bounds(V[3]),
        }
        assert (
            result["absolute_data_processing_loss"][1] >= 0
            and result["symmetric_data_processing_loss"][1] >= 0
        )
        for held in (2, 3):
            sub = [0, 1, held]
            a, v = A[sub], V[sub]
            w = reg.cross(a, v)
            factor = max(float(abs(x).b) for x in w)
            wn = [x / factor for x in w]
            label = ("matched_joint", "decorrelation")[held - 2]
            result[label] = {
                k: reg.measure_one(score[sub], a, v, w) for k, score in meters.items()
            }
            stats, vals = reg.evaluate_partitions(total_table[sub], wn)
            admitted &= (vals[:, 0] <= 0) & (vals[:, 1] >= 0)
            stats.pop("full_time")
            stats.pop("echo_only")
            result[label]["ordinary_plus_lost_time"] = stats
            row = [reg.score_balance(z[sub], wn) for z in (base, incoming, transmitted)]
            minors["pair"].append(row[1:])
            minors["with_ordinary"].append(row)
        result["necessary_any_channel_gain"] = necessary_gain(
            A, V, corrected, orders, pairs
        )
        quartets.append(result)
    witnesses = {}
    for name, rows in minors.items():
        witness = reg.rank_witness(rows)
        if witness is not None:
            indices = witness.pop("quartet_indices")
            witness["rows"] = [
                {
                    "stations": quartets[i // 2]["stations"],
                    "intervention": ("matched_joint", "decorrelation")[i % 2],
                }
                for i in indices
            ]
        witnesses[name] = witness
    return {
        "input_sha256": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in [*paths, HERE / "PRODUCT_INTERVENTION.md"]
        },
        "quartets": quartets,
        "input_log_order": [list(pair) for pair in orders],
        "incoming_modular_scores": incoming.tolist(),
        "ordinary_one_head_scores": base.tolist(),
        "unresolved_return_reference_blocks": unresolved,
        "surviving_classical_time_conventions": int(admitted.sum()),
        "common_meter_rank": witnesses,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    data = run()
    dest = HERE / "transmission_balance.json"
    if cli.parse_args().write:
        dest.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dest.read_bytes()), data)
    print(data["surviving_classical_time_conventions"], data["common_meter_rank"])
