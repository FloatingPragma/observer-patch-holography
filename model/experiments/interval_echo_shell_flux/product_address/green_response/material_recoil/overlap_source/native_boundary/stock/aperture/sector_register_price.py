"""Determine the modular price of summed exterior clocks from the actual reference law."""

import argparse
import itertools
import json
import math
from functools import cache
from hashlib import sha256

import numpy as np
from scipy.special import logsumexp

from . import sector_fronts_lab
from . import sector_register_pairs as pairs

HERE = pairs.collective.HERE


@cache
def logchoose(n):
    terms = [1]
    for i in range(n):
        terms.append(terms[-1] * (n - i) // (i + 1))
    return np.log(np.array([float(z) for z in terms]))


def kernel(n, s, t):
    k = np.arange(max(0, s + t - n), min(s, t) + 1)
    return k, logchoose(s)[k] + logchoose(n - s)[t - k]


def score(k, weights, ell):
    v = weights + ell * k
    norm = logsumexp(v)
    return -float(norm), float(k @ np.exp(v - norm))


def mixed(n, s, t, ell):
    rect = [kernel(n, s + a, t + b) for a, b in ((0, 0), (0, 1), (1, 0), (1, 1))]

    def eval_at(e):
        return np.array([score(k, w, e) for k, w in rect]).T

    sign = np.array([1, -1, -1, 1])
    # For each cell d(log Z)/d(log odds)=E[k], increasing with the log
    # odds since its derivative is Var(k). Enclose the derivative of the
    # signed four-cell contrast throughout every calibration subinterval.
    grid = np.linspace(ell[0], ell[1], 65)
    previous = eval_at(grid[0])
    bounds = []
    for x, y in itertools.pairwise(grid):
        current = eval_at(y)
        lo, hi = previous[1] - 1e-8, current[1] + 1e-8
        rate_lo = (-sign).clip(min=0) @ lo + (-sign).clip(max=0) @ hi
        rate_hi = (-sign).clip(min=0) @ hi + (-sign).clip(max=0) @ lo
        base = sign @ previous[0]
        bounds.append(
            [
                base + min(0, rate_lo) * (y - x) - 2e-9,
                base + max(0, rate_hi) * (y - x) + 2e-9,
            ]
        )
        previous = current
    nominal = float(sign @ eval_at(np.mean(ell))[0])
    return nominal, [min(z[0] for z in bounds), max(z[1] for z in bounds)]


def observe(data, n, ell):
    raw = np.array(data, int)
    xy = (raw[:, :2] - 12 * n) // 2
    np.testing.assert_array_equal(raw[:, :2], 12 * n + 2 * xy)
    x, ix = np.unique(xy[:, 0], return_inverse=True)
    y, iy = np.unique(xy[:, 1], return_inverse=True)
    counts = np.zeros((len(x), len(y)), int)
    counts[ix, iy] = raw[:, 2]
    K = np.array([[score(*kernel(n, int(a), int(b)), ell)[0] for b in y] for a in x])
    # Remove two arbitrary marginal prices before subtracting large means.
    K = K - K[:, [len(y) // 2]] - K[[len(x) // 2], :] + K[len(x) // 2, len(y) // 2]
    return connected(K, counts), K, counts


def connected(K, counts):
    f, g = counts.sum(1), counts.sum(0)
    m = int(counts.sum())
    d = float(np.sum(K * counts))
    u = f @ K @ g
    val = (m * d - u) / (m * (m - 1))
    # Delete one register, including it once less in each marginal. This
    # gives the jackknife sampling standard error for the U-statistic.
    jk = ((m - 1) * (d - K) - (u - (K @ g)[:, None] - (f @ K)[None, :] + K)) / (
        (m - 1) * (m - 2)
    )
    mean = np.sum(counts * jk) / m
    se = np.sqrt((m - 1) / m * np.sum(counts * (jk - mean) ** 2))
    return {
        "n_registers": m,
        "connected_modular_nats": float(val),
        "jackknife_standard_error": float(se),
    }


def native_control(flags, blank, n, calibration, center):
    # One member switches only X and another switches only Y. Every other
    # physical word is held. Both arms have the same law at every member.
    a, b = center
    ids = calibration["control_word_ids"]
    both, zero = ids[0]
    onlyx, onlyy = ids[1]
    assert flags[[zero, onlyx, onlyy, both]].tolist() == [
        [0, 0],
        [1, 0],
        [0, 1],
        [1, 1],
    ]
    shared = min(a, b)
    stock = (
        [both] * shared
        + [onlyx] * (a - shared)
        + [onlyy] * (b - shared)
        + [zero] * (n - 2 - max(a, b))
    )
    assert len(stock) == n - 2
    out = []
    for v in ((zero, zero), (onlyx, onlyy), (onlyx, zero), (zero, onlyy)):
        totals = flags[stock + list(v)].sum(0)
        out.append((blank * n + 2 * totals).tolist())
    return out


def study():
    input = HERE / "sector_fronts_collective.json"
    passage = HERE / "sector_passage.json"
    raw = json.loads(input.read_bytes())
    src = json.loads(passage.read_bytes())["spatial_modular_source"]
    n = raw["size"]
    cal = src["thermal_source"]["references"]

    flags, routes = sector_fronts_lab.two_delay_channels(pairs.chart)
    cells = 2 * flags[:, 0] + flags[:, 1]
    ref = {"ordinary": np.bincount(cells) / len(cells)}
    bath = pairs.chart.source.thermal.bath()
    ref.update(
        {
            tag: np.bincount(cells, weights=bath[tag].diagonal())
            for tag in ("retained", "exchanged")
        }
    )

    def logodds(v):
        return math.log(v[0] * v[3] / (v[1] * v[2]))

    bands = {
        "ordinary": [logodds(ref["ordinary"])] * 2,
        **{k: cal[k]["log_odds"] for k in cal},
    }
    centers = [
        (int(n * flags[:, 0].mean()), int(n * flags[:, 1].mean())),
        (n // 2, n // 2),
    ]
    rects = []
    for s, t in centers:
        rows = {k: mixed(n, s, t, e)[1] for k, e in bands.items()}
        relative = (np.array(rows["retained"]) - np.array(rows["exchanged"])[::-1]) / 2
        rects.append(
            {
                "base_total": [s, t],
                "ordinary_clock_control": native_control(
                    flags, routes["blank_duration"], n, src["reference"], (s, t)
                ),
                "mixed_price_per_covariance": rows,
                "matched_relative_modular_change": relative.tolist(),
            }
        )
    results = {}
    for arm in ("peer", "ordinary"):
        joint = raw["native_acquisition"][arm + "_summed_times_and_multiplicity"]
        outs = {}
        tabs = {}
        for tag, p in ref.items():
            row, K, counts = observe(joint, n, logodds(p))
            outs[tag] = row
            tabs[tag] = K
        # The same terminal records enter both reference scores; their
        # uncertainty must be computed on the difference for each deletion.
        diff = tabs["retained"] - tabs["exchanged"]
        outs["relative"] = connected(diff, counts)
        results[arm] = outs
    contrast = {
        k: {
            "nats_per_register": results["peer"][k]["connected_modular_nats"]
            - results["ordinary"][k]["connected_modular_nats"],
            "sampling_standard_error": math.hypot(
                results["peer"][k]["jackknife_standard_error"],
                results["ordinary"][k]["jackknife_standard_error"],
            ),
        }
        for k in (*ref, "relative")
    }
    return {
        "sha256": {
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (input, passage, HERE / "SECTOR_REGISTER_BALLS.md")
        },
        "exterior_receiver": routes,
        "matched_corners": rects,
        "whole_register_measured": results,
        "peer_minus_ordinary": contrast,
        "single_member_mixed_relative_price": src["thermal_source"]["mixed_price"],
        "covariance_contrast_per_member": raw["native_acquisition"][
            "both_arms_measured_covariance_excess"
        ],
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    out = study()
    file = HERE / "sector_register_price.json"
    if p.parse_args().write:
        file.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        pairs.chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(file.read_text()), out
        )
    print(out["matched_corners"], out["peer_minus_ordinary"])
