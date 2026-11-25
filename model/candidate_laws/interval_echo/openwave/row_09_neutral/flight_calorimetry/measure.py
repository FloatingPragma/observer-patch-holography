"""Read the moving pair's absolute formation energy and its subsequent charge contrast."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta

from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release import (
    measure as release,
)
from .certificate import exact_support
from .apparatus import (
    EDGES,
    RING,
    SENSORS,
    clock,
    frame,
    launch,
    loads,
    local_meter,
    sector,
    unfiltered,
)

HERE = Path(__file__).resolve().parent
STARTS = 2**46
LOCAL_STARTS = 2**38
FAMILY = 4096


def interval(k, n):
    tail = 0.01 / (2 * FAMILY)
    return [
        float(beta.ppf(tail, k, n - k + 1)) if k else 0.0,
        float(beta.isf(tail, k + 1, n - k)) if k < n else 1.0,
    ]


def draw(p, source, setting, n, rng, bench=False):
    envelope = 1.0 if bench else 1 / 3
    all_p = np.r_[
        source * p.ravel(), source * (1 - setting), envelope - source, 1 - envelope
    ]
    assert min(all_p) >= -2e-14 and abs(sum(all_p) - 1) < 2e-12
    all_p = np.maximum(all_p, 0)
    return rng.multinomial(n, all_p / all_p.sum())


def timing_bank(
    archive, name, kind, contrast, zero, sign, k, size, prior, rng, bench=False
):
    accept, readings = clock(kind, contrast, zero, sign, size)
    bins = np.arange(size) * 128 / size
    absolute = np.minimum(bins, 128 - bins)
    out = []
    for detector, p in readings.items():
        key = f"{name}_{sign}_{k}_{int(zero)}_{size}_{detector}"
        bank = draw(p, prior, accept, STARTS, rng, bench)
        archive[key] = bank
        observed = bank[:-3].reshape(size, 5)
        if kind == "pair":
            joint = observed[:, [0, 4]].sum(axis=1)
        elif kind == "host":
            joint = observed[:, 2]
        else:
            joint = observed[:, [1, 3]].sum(axis=1)
        total = int(joint.sum())
        soft = int(joint[absolute <= 1.6].sum())
        high = int(joint[(bins >= 44) & (bins <= 58)].sum())
        top = int(np.argmax(joint))
        eligible = np.minimum(abs(bins - bins[top]), 128 - abs(bins - bins[top])) >= 0.5
        second = max(np.flatnonzero(eligible), key=joint.__getitem__)
        out.append(
            dict(
                key=key,
                sensor=detector,
                starts=STARTS,
                accepted=int(observed.sum()),
                joint=total,
                contrast_two=int(observed[:, [0, 4]].sum()),
                soft=soft,
                high=high,
                soft_share_99=interval(soft, total),
                high_share_99=interval(high, total),
                peak_bins=[
                    [float(bins[i]), int(joint[i])]
                    for i in (top, second)
                    if joint[i] and (i == top or eligible[i])
                ],
                refusals=bank[-3:].tolist(),
            )
        )
    return out


def near_meter(sign, prior, k, rng):
    words, p = unfiltered("pair", sign)
    out = []
    for sensor, (a, b) in SENSORS.items():
        ordered = [(n, loads(w, (a, b))) for n, w in enumerate(words)]
        hit = np.array([abs(x - y) == 2 for _, (x, y) in ordered])
        local = [e for e in EDGES if a in e or b in e]
        assert len(local) == 3
        cohort = []
        for v, z in local:
            d = np.array([np.subtract(*loads(w, (v, z))) for w in words])
            prob = np.array([sum(p[(d == i) & hit]) for i in range(-2, 3)])
            full = draw(np.r_[prob, p[~hit].sum()], prior, 1, LOCAL_STARTS, rng)
            assert full[-3] == 0
            cohort.append(dict(edge=[v, z], counts=full.tolist()))
        nominal = [4 if set(rec["edge"]) == {a, b} else 1 for rec in cohort]
        lower, upper = 0.0, 0.0
        for rec, value in zip(cohort, nominal, strict=True):
            c = rec["counts"][:5]
            hit = sum(n for d, n in zip(range(-2, 3), c, strict=True) if d * d == value)
            q = interval(hit, sum(c))[0]
            lower += 8 * value * q
            upper += 8 * (value * q + 4 * (1 - q))
        out.append(
            dict(
                sensor=sensor,
                sign=sign,
                k=k,
                cohorts=cohort,
                post_tag_contrast_99=[lower, upper],
            )
        )
    return out


def run():
    local_meter()
    weights = release.source_records().sum(axis=0)
    rng = np.random.default_rng(909244)
    counts = {}
    timing = []
    direct = []
    for sign in (-1, 1):
        for k, prior in enumerate(weights):
            direct.extend(near_meter(sign, prior, k, rng))
            for size in (1024, 4096):
                for zero in (False, True):
                    timing.extend(
                        timing_bank(
                            counts, "native", "pair", 8, zero, sign, k, size, prior, rng
                        )
                    )
                timing.extend(
                    timing_bank(
                        counts, "exchange", "pair", 0, True, sign, k, size, prior, rng
                    )
                )
            kind = "vacancy" if sign < 0 else "excess"
            timing.extend(
                timing_bank(
                    counts, "depleted", kind, 8, True, sign, k, 4096, prior, rng
                )
            )
    for kind, sign in (("vacancy", -1), ("excess", 1), ("host", 0)):
        # Fresh sharp compensated stock, with the host itself as the energy zero.
        timing.extend(
            timing_bank(
                counts, "standard", kind, 8, True, sign, 0, 4096, 1.0, rng, True
            )
        )
    full = {}
    for kind, contrast in (
        ("pair", 8),
        ("pair", 0),
        ("vacancy", 8),
        ("excess", 8),
        ("host", 8),
    ):
        words, h, cost, _ = sector(kind, contrast)
        ev, modes = frame(kind, contrast, True)
        full[f"{kind}_{contrast}"] = dict(
            states=len(words),
            contrast_values=sorted(set(cost.tolist())),
            zero_phase_values=ev.tolist(),
            source_spectral_weights=(
                abs(modes[words.index(launch(kind, -1))]) ** 2
            ).tolist(),
            smallest_all_phase_energy=float(np.linalg.eigvalsh(h)[0]),
        )
    return (
        dict(
            finite_certificate=exact_support(),
            ring=RING,
            source_acceptance=weights.tolist(),
            full=full,
            direct=direct,
            timing=timing,
            seeds=[909244],
            clock_columns="for each increasing clock bin, signed d = -2,-1,0,1,2; then rejected spatial phase, original workstation miss, original source refusal",
            direct_columns="d=-2,-1,0,1,2 after magnitude-two coincidence, no coincidence, spatial rejection (unused), workstation miss, source refusal",
        ),
        counts,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    summary, counts = run()
    if args.write:
        np.savez_compressed(HERE / "counts.npz", **counts)
        (HERE / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"{len(counts)} timing banks, {len(summary['direct'])} three-contact calorimetry banks"
    )
