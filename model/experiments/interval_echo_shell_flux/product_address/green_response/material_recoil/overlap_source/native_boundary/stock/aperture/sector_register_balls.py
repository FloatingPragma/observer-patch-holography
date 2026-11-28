"""Measure entropy tangents and equal-volume boundaries on a common two-member chart."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from . import sector_register_pairs as pairs

HERE = pairs.collective.HERE


def timed_charts(adj, uv, cfg, periods):
    u, v = uv.T
    types = np.array([periods[tuple(c)] for c in cfg])
    times, group = np.unique(
        (types[:, None] + types).reshape(-1, len(u)), axis=0, return_inverse=True
    )
    ranks = []
    for t in times:
        graph = adj.copy()
        graph[u, v] = graph[v, u] = t
        d = shortest_path(csr_matrix(graph)).astype(int)
        assert np.all(d[:, v] <= d[:, u] + t) and np.all(d[:, u] <= d[:, v] + t)
        ranks.append(np.argsort(np.argsort(d, axis=1, kind="stable"), axis=1))
    return np.array(ranks), group, times


def fixed_for_every_register(ports, adj, uv, center, volume):
    """Certify the same ranked ball for any mean of the measured link periods."""
    c = np.searchsorted(ports, center)
    homes = np.searchsorted(ports, [8, 20, 32, 44])
    u, v = uv.T
    maximal = adj.copy()
    maximal[u, v] = maximal[v, u] = (
        1 + np.isin(u, homes).astype(int) + np.isin(v, homes).astype(int)
    )
    low, high = [
        shortest_path(csr_matrix(g), indices=c).astype(int) for g in (adj, maximal)
    ]
    order = np.argsort(low, kind="stable")
    inside, outside = order[:volume], order[volume:]
    last = max(zip(high[inside], ports[inside], strict=True))
    first = min(zip(low[outside], ports[outside], strict=True))
    cut = np.isin(u, inside) != np.isin(v, inside)
    assert last < first and np.isin(homes, inside).all()
    assert not np.isin(uv[cut], homes).any()
    return {
        "center": center,
        "volume": volume,
        "largest_inside_upper_time_and_address": list(map(int, last)),
        "smallest_outside_lower_time_and_address": list(map(int, first)),
        "constant_one_period_boundary_links": ports[uv[cut]].tolist(),
    }


def prices(labels, ref, bands):
    _, g = np.unique(labels, return_inverse=True)
    g = g.ravel()

    def freq(p):
        return np.bincount(g, weights=np.outer(p, p).ravel())

    lo, hi = freq(bands[:, 0]), freq(bands[:, 1])
    assert min(lo) > 0
    K = -np.log(freq(ref))
    return g, K, -np.log(hi) - 1e-12, -np.log(lo) + 1e-12


def matched_home_reading(ports, uv, cfg, periods, refs):
    """Time four internal links on matched arms, keeping only station totals."""
    w = pairs.chart.source.thermal.apparatus()[0]
    homes = np.searchsorted(ports, [8, 20, 32, 44])
    probe_links = [
        np.flatnonzero(((uv == h).any(1)) & (np.isin(uv, homes).sum(1) == 1))[0]
        for h in homes
    ]
    src = json.loads((HERE / "sector_passage.json").read_text())[
        "spatial_modular_source"
    ]["reference"]
    (_, zero), (only_x, only_y) = src["control_word_ids"]
    # Member A switches only the first exterior clock, member B only the
    # second. The full local timing record may also resolve the other homes.
    arms = [[(zero, zero), (only_x, only_y)], [(only_x, zero), (zero, only_y)]]
    ticks = []
    for arm in arms:
        rows = [
            sum(periods[tuple(w[i, :4] > 0)][probe_links] for i in members)
            for members in arm
        ]
        ticks.append(np.array(rows))
    # Decode using only protected link times, with one empty period per member.
    ids = ((cfg[:, None, :] + cfg[None, :, :]) @ np.array([1, 3, 9, 27])).ravel()
    outcomes = np.unique(ids)
    signed = np.zeros(len(outcomes))
    for sign, arm in zip([1, -1], ticks, strict=True):
        keys = (arm - 2) @ np.array([1, 3, 9, 27])
        loc = np.searchsorted(outcomes, keys)
        np.testing.assert_array_equal(outcomes[loc], keys)
        np.add.at(signed, loc, sign / 2)
    readings = {}
    for name, ref in refs.items():
        _, k, lower, upper = prices(ids, *ref)
        readings[name] = [
            float(signed @ k),
            [
                float(signed.clip(min=0) @ lower + signed.clip(max=0) @ upper),
                float(signed.clip(min=0) @ upper + signed.clip(max=0) @ lower),
            ],
        ]
    lo, hi = readings["retained"][1]
    l2, h2 = readings["exchanged"][1]
    return {
        "word_pairs_by_arm": arms,
        "links": ports[uv[probe_links]].tolist(),
        "summed_periods_by_arm": np.array(ticks).tolist(),
        "entropy_price_by_reference": readings,
        "relative_modular_contrast_nats": [lo - h2, hi - l2],
    }


def tangent(g, K, lo, hi, delta, error):
    gain = np.bincount(g, weights=delta.ravel())
    val = gain @ K
    off = error * (hi.max() - lo.min()) + 2e-13
    lower = np.maximum(gain, 0) @ lo + np.minimum(gain, 0) @ hi - off
    upper = np.maximum(gain, 0) @ hi + np.minimum(gain, 0) @ lo + off
    return float(val), [float(lower), float(upper)]


def survey():
    (ports, clear, adj, uv, cfg, _ids, periods), refs, mu, delta, errors = (
        pairs.apparatus()
    )
    R, group, time = timed_charts(adj, uv, cfg, periods)
    homes = np.searchsorted(ports, [8, 20, 32, 44])
    ext = np.flatnonzero(clear >= 2)
    ordinary = shortest_path(csr_matrix(adj), unweighted=True)
    ranking = np.argsort(np.argsort(ordinary, axis=1, kind="stable"), axis=1)
    u, v = uv.T
    m = len(cfg)
    n = len(ports)
    A, B = np.repeat(cfg, m, axis=0), np.tile(cfg, (m, 1))
    rows = []
    for center in ext:
        r = 1 + max(ranking[center, homes[[0, 3]]])
        all_homes = 1 + R[:, center, homes].max()
        for volume in sorted({r - 1, r, r + 1, all_homes} & set(range(1, n))):
            member = R[:, center] < volume
            assert np.all(member.sum(1) == volume)
            cross = member[:, u] != member[:, v]
            area = cross.sum(1)[group]
            # Six times the service rate. The mean of two link times is t/2.
            capacity = (cross * (12 // time)).sum(1)[group] / 6
            geom = np.unique(member, axis=0, return_inverse=True)[1][group]
            visible = member[group][:, homes]
            fine = (
                geom * 256
                + (visible * A) @ np.array([1, 2, 4, 8])
                + 16 * ((visible * B) @ np.array([1, 2, 4, 8]))
            )
            coarse = geom * 81 + (visible * (A + B)) @ np.array([1, 3, 9, 27])
            row = {
                "center": int(ports[center]),
                "volume": int(volume),
                "fixed_members": bool(np.ptp(geom) == 0),
                "always_all_homes": bool(visible.all()),
                "distinct_boundaries": int(np.ptp(geom) + 1),
                "unchanged_cut_transcript": bool(
                    np.ptp(geom) == 0 and np.all(np.ptp(time[:, cross[0]], axis=0) == 0)
                ),
            }
            for tag, d in delta.items():
                row[tag] = {
                    "area": float(area @ d.ravel()),
                    "capacity": float(capacity @ d.ravel()),
                    "area_error": float(np.ptp(area) * errors[tag] + 2e-13),
                    "capacity_error": float(np.ptp(capacity) * errors[tag] + 2e-13),
                }
                for metric, readings in (("area", area), ("capacity", capacity)):
                    if np.ptp(readings) == 0:
                        row[tag][metric] = row[tag][metric + "_error"] = 0.0
            for label, reading in (("ordered", fine), ("pooled", coarse)):
                for key, (p, bands) in refs.items():
                    price = prices(reading, p, bands)
                    for tag, d in delta.items():
                        value, band = tangent(*price, d, errors[tag])
                        if (
                            label == "ordered"
                            and row["fixed_members"]
                            and tag == "cross"
                        ):
                            assert abs(value) < errors[tag] * 30 + 1e-12
                            # Exact factorization holds for every faithful
                            # one-member price, including the calibration box.
                            value, band = 0.0, [0.0, 0.0]
                        row[tag][f"{label}_{key}"] = [value, band]
            rows.append(row)
    sealed = [
        r for r in rows if r["always_all_homes"] and r["unchanged_cut_transcript"]
    ]
    witness = min(sealed, key=lambda r: (r["volume"], r["center"]))
    certificate = fixed_for_every_register(
        ports, adj, uv, witness["center"], witness["volume"]
    )
    return (
        rows,
        cfg,
        mu,
        delta,
        errors,
        certificate,
        matched_home_reading(ports, uv, cfg, periods, refs),
    )


def summarize(rows):
    out = {}
    for label in ("ordered", "pooled"):
        key = label + "_ordinary"
        for metric in ("area", "capacity"):
            flags = {
                k: []
                for k in (
                    "fixed_metric_nonzero_entropy",
                    "ordinary_positive_coefficient",
                    "calibrated_stationarity_rejected",
                )
            }
            for r in rows:
                c, o = r["cross"], r["ordinary"]
                sig = c[key][1][0] > 0 or c[key][1][1] < 0
                if c[metric + "_error"] == 0 and sig:
                    flags["fixed_metric_nonzero_entropy"].append(r)
                if abs(o[metric]) <= o[metric + "_error"]:
                    continue
                signs = np.array([-1, 1])
                ordinary_metric = o[metric] + signs * o[metric + "_error"]
                eta = -np.divide.outer(o[key][1], ordinary_metric)
                if eta.min() <= 0:
                    continue
                flags["ordinary_positive_coefficient"].append(r)
                source_metric = c[metric] + signs * c[metric + "_error"]
                products = np.multiply.outer([eta.min(), eta.max()], source_metric)
                residual = c[key][1] + np.array([products.min(), products.max()])
                if residual[0] > 1e-12 or residual[1] < -1e-12:
                    flags["calibrated_stationarity_rejected"].append(r)
            out[label + "/" + metric] = {k: len(x) for k, x in flags.items()}
    return out


def run():
    rows, cfg, mu, delta, errors, certificate, control = survey()
    key = "pooled_ordinary"

    def strength(r):
        z = r["cross"][key][1]
        return max(0, z[0], -z[1])

    fixed = [r for r in rows if r["fixed_members"] and r["always_all_homes"]]
    sealed = [r for r in fixed if r["unchanged_cut_transcript"]]
    examples = [
        max(fixed, key=strength),
        max((r for r in rows if not r["fixed_members"]), key=strength),
    ]
    if sealed:
        examples.append(min(sealed, key=lambda r: (r["volume"], r["center"])))
    # Identical membership and covered homes give identical measured outcomes.
    # Keep one complete table of cases for every center in the replay code
    # and archive the ranges, counts, ordinary calibrations and named witnesses.
    shared = sorted({(r["center"], r["volume"]) for r in examples})
    for r in rows:
        c, o = r["cross"], r["ordinary"]
        if o["area"] * o[key][0] < -1e-12 and abs(c["area"]) > 1e-10:
            examples.append(r)
            break
    return {
        "input_sha256": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "SECTOR_REGISTER_BALLS.md",
                "sector_fronts_collective.json",
                "sector_calorimetry.json",
                "sector_passage.json",
            )
        },
        "all_register_size_certificate": certificate,
        "matched_local_clock_control": control,
        "cases": len(rows),
        "centers": len({r["center"] for r in rows}),
        "stable_enclosing_all": len(fixed),
        "stable_enclosing_all_with_unchanged_cut_transcript": len(sealed),
        "unique_primary_examples": shared,
        "summary": summarize(rows),
        "examples": examples,
        "mean_occupancy_configs": mu.tolist(),
        "configs": cfg.tolist(),
        "distinct_minus_independent": delta["cross"].tolist(),
        "source_errors_l1": errors,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    s = run()
    dst = HERE / "sector_register_balls.json"
    if opt.write:
        dst.write_text(json.dumps(s, indent=2, allow_nan=False) + "\n")
    else:
        pairs.chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dst.read_text()), s
        )
    print(s["cases"], s["stable_enclosing_all"], s["summary"])
