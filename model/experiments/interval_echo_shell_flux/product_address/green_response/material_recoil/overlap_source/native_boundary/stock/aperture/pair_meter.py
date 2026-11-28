"""Infer an unknown connected charge from local ordinary-pair interference alone."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import expm

from . import pair_contact as source

HERE = source.HERE
SHOTS = 10**6
# Covers all four reference bins, calibration quadratures, source comparisons
# and independent head checks at simultaneous 99% confidence.
LOG_RISK = np.log(2 * 128 / 0.01)


def estimate(ref, clocks, native, held):
    assert (
        sum(ref) == SHOTS
        and min(ref) >= 0
        and 0 <= min(native, held) <= max(native, held) <= SHOTS
    )
    n = sum(ref[:4])
    f = np.asarray(ref[:4]) / n
    dr = np.sqrt(LOG_RISK / (2 * n))
    lo, hi = f - dr, f + dr
    assert min(lo) > 0
    odds = np.log([lo[0] * lo[3] / (hi[1] * hi[2]), hi[0] * hi[3] / (lo[1] * lo[2])])
    err = np.sqrt(2 * LOG_RISK / SHOTS)
    phases = []
    for x, y in clocks:
        assert min(x, y) >= 0 and max(x, y) <= SHOTS
        phases.append(complex(2 * x / SHOTS - 1, 2 * y / SHOTS - 1))
    distance = abs(1 - np.mean(phases))
    gain = (distance + np.array([-1, 1]) * np.sqrt(2) * err) ** 2 / 2
    assert distance > np.sqrt(2) * err and max(odds) < 0
    delta = (held - native) / SHOTS + np.array([-1, 1]) * err
    samples = [-o * d / g for o in odds for d in delta for g in gain]
    return [float(min(samples)), float(max(samples))]


def run():
    _, table, rejected, total = source.source_counts()
    p = table.ravel() / table.sum()
    q = np.outer(table.sum(1), table.sum(0)).ravel() / table.sum() ** 2
    rng = np.random.default_rng(537103)
    reference = rng.multinomial(SHOTS, np.r_[table.ravel(), rejected] / total).tolist()
    B, D = source.blocks()["full"][2]
    phase = expm(1j * source.TAU * B / 4) @ expm(-1j * source.TAU * (B + D) / 4)
    clocks = []
    for coin in (1, -1):
        for sign in (1, -1):
            u = np.array([1, sign, coin, sign * coin]) / 2
            visibility = u @ phase @ u
            clocks.append(
                [
                    int(rng.binomial(SHOTS, np.clip((1 + y) / 2, 0, 1)))
                    for y in (visibility.real, visibility.imag)
                ]
            )
    g = phase[0, 0]
    bright = np.array([1, (1 + g.real) / 2, (1 + g.real) / 2, (1 + abs(g) ** 2) / 2])
    readings = []
    for frac in (0, 0.5, 1):
        held = (1 - frac) * p + frac * q
        old, new = [int(rng.binomial(SHOTS, d @ bright)) for d in (p, held)]
        interval = estimate(reference, clocks, old, new)
        # Scoring head readings are available only after estimating the charge.
        head = rng.multinomial(SHOTS, held).tolist()
        expected = float((held - p) @ np.log(q / p))
        assert interval[0] <= expected <= interval[1]
        # Every four-bin score is bounded by the range of log(Q/P).
        weights = np.log(q / p)
        score = float((np.array(head) / SHOTS - p) @ weights)
        width = np.ptp(weights) * np.sqrt(LOG_RISK / (2 * SHOTS))
        readings.append(
            dict(
                frac=frac,
                bright=[old, new],
                from_ordinary_fringe=interval,
                scoring_head_counts=head,
                known_reference_head_interval=[score - width, score + width],
                exact_nats=expected,
            )
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "PAIR_METER.md", HERE / "pair_contact.json")
        },
        trials=SHOTS,
        reference=reference,
        clocks=clocks,
        readings=readings,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "pair_meter.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        source.timing.lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print([(r["frac"], r["from_ordinary_fringe"]) for r in out["readings"]])
