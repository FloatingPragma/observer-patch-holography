"""Visible forward and reverse paths of two finite material repairs."""

import argparse
from collections import Counter
from fractions import Fraction as Fr
import json
from math import comb, log, pi, sqrt

import numpy as np
from scipy.stats import beta

from .gibbs import (
    INDEX, LOCAL, PACKED, REFERENCE, SOURCE, native_exchange, run_shots, transposition,
)
from .measure import HERE, physical_path

SHOTS = 2**16
# In total the Hoeffding and simultaneous marginal-bin intervals fail with
# probability at most 0.01 under independent launches of the declared model.
EPS = sqrt(log(2 * 1024 / 0.005) / (2 * SHOTS))


def make_repair(port):
    route = [] if port == 3 else [(e, 1, 1) for e in transposition(port, 3)]
    return route + LOCAL + route[::-1]


def visible_repair(x, port):
    held = x & ~(1 << port)
    p = Fr(4 - held.bit_count(), 5)
    return ((held, 1 - p), (held | (1 << port), p))


def forecast():
    r = [Fr(comb(4, 4 - x.bit_count()), 70) for x in range(16)]
    p = [a * (Fr(7, 16) if (x & 3).bit_count() == 1 else Fr(7, 4))
         for x, a in enumerate(r)]
    q = [Fr(0)] * 16
    forward = {}
    for x, a in enumerate(p):
        for m, b in visible_repair(x, 0):
            for y, c in visible_repair(m, 2):
                forward[x, m, y] = a * b * c
                q[y] += a * b * c
    reverse = {}
    for y, a in enumerate(q):
        for m, b in visible_repair(y, 2):
            for x, c in visible_repair(m, 0):
                reverse[x, m, y] = a * b * c
    odds = {k: p[x] * r[y] / (q[y] * r[x]) for k in forward for x, m, y in [k]}
    return r, p, q, forward, reverse, odds


def prepare(rng, state):
    if state == "R":
        return run_shots(rng, np.full(SHOTS, INDEX[60]), REFERENCE)
    selector = rng.random(SHOTS) < native_exchange(pi / 3)
    source = run_shots(rng, np.where(selector, INDEX[29], INDEX[60]), SOURCE, complement=True)
    if state == "Q":
        for port in (0, 2):
            source = run_shots(rng, source, make_repair(port))
    return source


def acquire(seed, start, schedule, parked=False):
    rng = np.random.default_rng(seed)
    indices = prepare(rng, start)
    codes = PACKED[indices] & 15
    for port in schedule:
        if parked:
            if native_exchange(pi, disabled=True) != 0:
                raise ValueError("Disabled exchange moved an occupation")
        else:
            indices = run_shots(rng, indices, make_repair(port))
        codes = codes * 16 + (PACKED[indices] & 15)
    hist = np.bincount(codes, minlength=16 ** (len(schedule) + 1))
    return [[int(k), int(hist[k])] for k in np.flatnonzero(hist)]


def encode(x, m, y):
    return 256 * x + 16 * m + y


def interval(n):
    alpha = 0.005 / 176
    lower = 0.0 if n == 0 else beta.ppf(alpha / 2, n, SHOTS - n + 1)
    upper = 1.0 if n == SHOTS else beta.ppf(1 - alpha / 2, n + 1, SHOTS - n)
    return float(lower), float(upper)


def summarize(banks):
    r, p, q, forward, reverse, odds = forecast()
    hist = {k: dict(v) for k, v in banks.items()}
    cal = {lab: [interval(hist[lab].get(x, 0)) for x in range(16)] for lab in ("P", "Q", "R")}
    f = {k: hist["forward"].get(encode(*k), 0) for k in forward}
    b = {k: hist["reverse"].get(encode(*k[::-1]), 0) for k in forward}
    groups = {}
    negative_by_calibration = 0
    overlap = 0
    for k, v in odds.items():
        x, _, y = k
        entry = groups.setdefault(v, [Fr(0), Fr(0), 0, 0])
        entry[0] += forward[k]
        entry[1] += reverse[k]
        entry[2] += f[k]
        entry[3] += b[k]
        px, qy, rx, ry = cal["P"][x], cal["Q"][y], cal["R"][x], cal["R"][y]
        s = [px[0] * ry[0] / (qy[1] * rx[1]), px[1] * ry[1] / (qy[0] * rx[0])]
        if x == y:  # The reference cancels identically.
            s = [px[0] / qy[1], px[1] / qy[0]]
        l, u = interval(f[k])
        a, z = interval(b[k])
        overlap += int(max(l / z, s[0]) <= min(u / a, s[1]))
        if s[1] < 1:
            negative_by_calibration += f[k]
    bad_reverse = sum(n for c,n in banks["wrong_reverse"] if (c % 16, c // 16 % 16, c // 256) not in forward)
    neg = sum(f[k] for k,v in odds.items() if v < 1)
    neg_exact = sum(a for k,a in forward.items() if odds[k] < 1)
    integral = sum(f[k] / float(v) for k,v in odds.items()) / SHOTS
    avg = sum(f[k] * log(v) for k,v in odds.items()) / SHOTS
    kl = sum(float(a)*log(a/c) for a,c in zip(p,r)) - sum(float(a)*log(a/c) for a,c in zip(q,r))
    # Sum over the event selected by the sign of the predicted flat-reference error.
    plus = [k for k in forward if odds[k] > p[k[0]] / q[k[2]]]
    foil = sum(f[k] - b[k]*float(p[k[0]]/q[k[2]]) for k in plus) / SHOTS
    foil_radius = EPS * (1 + max(float(p[k[0]]/q[k[2]]) for k in plus))
    return dict(
        by_odds={str(v): {"forward": str(a), "reverse": str(h), "forward_count": i, "reverse_count": j}
                 for v,(a,h,i,j) in sorted(groups.items())},
        negative_forecast=str(neg_exact), negative_count=neg,
        negative_from_independent_calibration=negative_by_calibration,
        independent_odds_intervals_overlap=overlap,
        reverse_wrong_order_outside_support=bad_reverse,
        exp_minus_sigma=round(integral,9), mean_sigma=round(avg,9),
        kl_descent=round(kl,9), flat_reference_residual=round(foil,9),
        flat_reference_radius=round(foil_radius,9),
    )


def run():
    banks = {state: acquire(110721 + i, state, ()) for i,state in enumerate(("R","P","Q"))}
    for i, (label, source, schedule) in enumerate((
        ("forward","P",(0,2)), ("reverse","Q",(2,0)), ("wrong_reverse","Q",(0,2)),
        ("equilibrium","R",(0,2)), ("parked","P",(0,2)),
    )):
        banks[label] = acquire(810734 + i, source, schedule, parked=label == "parked")
    return dict(shots_per_bank=SHOTS, path=physical_path(0), word_bins=16, path_bins=4096,
                histograms=banks, statistics=summarize(banks))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = run()
    data = json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    dest = HERE / "fluctuation.json"
    if args.write:
        dest.write_text(data)
    elif dest.read_text() != data:
        raise ValueError("Fluctuation replay mismatch")
    print(receipt["statistics"])
