"""Cubical rods, bulk face registers and an ordinary material excitation."""

import json
from collections import Counter
from fractions import Fraction
from itertools import product
from math import log, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SIDE = 6
ORIGIN = (0, 0, 0)
TAGS = tuple(product((0, 1), repeat=3))
AXES = tuple(np.eye(3, dtype=int))
SWAPS = {
    (2, 0, 0): (0, 1, 0),
    (0, 1, 0): (2, 0, 0),
    (3, 0, 0): (1, 1, 0),
    (1, 1, 0): (3, 0, 0),
}
CUTS = {
    "ordinary_line": tuple((x, 0, 0) for x in range(4)),
    "held_sheet": tuple((x, y, 0) for x in (1, 2, 3) for y in (0, 1)),
    "overlap_pair": ((1, 0, 0), (2, 0, 0)),
}
STEP = Fraction(1, 8)
TRIALS = 2**32


def sector(tag):
    h, y, z = tag
    return h ^ (y & z)


def positions(cut, g):
    return tuple(SWAPS.get(x, x) if g else x for x in cut)


def interior_registers(cells):
    occupied = set(cells)
    return sorted(
        (a, axis)
        for a in cells
        for axis, step in enumerate(AXES)
        if tuple((a + step) % SIDE) in occupied
    )


def rod_survey(cells):
    # Each positively oriented material cube carries three unit rod vectors.
    # Compare complete corner point sets to remove faces with a second owner.
    faces = Counter()
    areas = {}
    rods = np.eye(3, dtype=int)
    volume = 0
    for corner in cells:
        volume += int(round(abs(np.linalg.det(rods))))
        for normal in range(3):
            a, b = [rods[j] for j in range(3) if j != normal]
            for sign in (0, 1):
                start = np.array(corner) + sign * rods[normal]
                face = tuple(
                    sorted(
                        tuple((start + i * a + j * b) % SIDE)
                        for i, j in product((0, 1), repeat=2)
                    )
                )
                faces[face] += 1
                areas[face] = int(round(np.linalg.norm(np.cross(a, b))))
    return sum(areas[k] for k, n in faces.items() if n == 1), volume


def tag_weights(arm, t=Fraction(0), biased=False):
    if biased:
        return [Fraction(2 ** sector(k), 12) for k in TAGS]
    if arm == "ordinary":
        return [(Fraction(1, 2) + (2 * h - 1) * t) / 4 for h, _, _ in TAGS]
    if arm == "held":
        return [(1 + t * (-1) ** sum(k)) / 8 for k in TAGS]
    raise ValueError(arm)


def controller_step(tag, axis, biased=False):
    """One reversible fair proposal on the three physical controller bits."""
    moved = tuple(value ^ (i == axis) for i, value in enumerate(tag))
    ratio = Fraction(2 ** sector(moved), 2 ** sector(tag)) if biased else Fraction(1)
    forward = min(Fraction(1), ratio) / 2
    return {tag: 1 - forward, moved: forward}


def histogram(cut, weights):
    # A geometry outcome records the locations of the enrolled material labels
    # in label order. All face bits interior to that geometry are recorded.
    probs = Counter()
    for tag, mass in zip(TAGS, weights, strict=True):
        cells = positions(cut, sector(tag))
        registers = interior_registers(cells)
        ctl = tag if ORIGIN in cut else ()
        for word in range(2 ** len(registers)):
            probs[cells, ctl, word] += mass / 2 ** len(registers)
    if sum(probs.values()) != 1:
        raise ValueError("incomplete regional probability law")
    return probs


def projection_slopes(cut, arm, biased=False):
    ref = histogram(cut, tag_weights(arm, biased=biased))
    plus = tag_weights(arm, STEP)
    minus = tag_weights(arm, -STEP)
    hi, lo = histogram(cut, plus), histogram(cut, minus)
    entropy = -sum(
        float((hi[k] - lo[k]) / (2 * STEP)) * log(float(p)) for k, p in ref.items()
    )
    dg = sum(sector(k) * (a - b) for k, a, b in zip(TAGS, plus, minus, strict=True)) / (
        2 * STEP
    )
    geo = [rod_survey(positions(cut, g)) for g in (0, 1)]
    area = dg * (geo[1][0] - geo[0][0])
    volume = dg * (geo[1][1] - geo[0][1])
    return entropy, area, volume, dg


def terminal_counts(cut, arm, rng):
    sectors = [rod_survey(positions(cut, g)) for g in (0, 1)]
    raw = {}
    for sign in (-1, 1):
        weights = tag_weights(arm, sign * STEP)
        distribution = histogram(cut, weights)
        counts = rng.multinomial(TRIALS, list(map(float, distribution.values())))
        area_hits = rng.binomial(
            TRIALS,
            float(sum(p * sector(tag) for p, tag in zip(weights, TAGS, strict=True))),
        )
        energy_hits = rng.binomial(
            TRIALS, float(sum(p * tag[0] for p, tag in zip(weights, TAGS, strict=True)))
        )
        raw[str(sign)] = {
            "bulk_counts_in_tag_then_word_insertion_order": counts.tolist(),
            "rod_sector_hits": int(area_hits),
            "ordinary_energy_hits": int(energy_hits),
        }
    return {"sector_area_volume": sectors, "records": raw}


def concentrated_entropy(counts):
    # Bounded differences for the plug-in estimator, plus its standard
    # elementary bias bound 0<=H(p)-E H(phat)<=log(1+(k-1)/N).
    n = sum(counts)
    freq = np.asarray([c for c in counts if c]) / n
    estimate = -float(freq @ np.log(freq))
    radius = (log(n) + 1) * sqrt(5 / n)
    bias = np.log1p((len(counts) - 1) / n)
    return [estimate - radius, estimate + radius + float(bias)]


def interval_subtract(a, b, scale):
    return [(a[0] - b[1]) / scale, (a[1] - b[0]) / scale]


def read(packet):
    minus, plus = (packet["records"][str(i)] for i in (-1, 1))
    step = float(2 * STEP)
    s = interval_subtract(
        concentrated_entropy(plus["bulk_counts_in_tag_then_word_insertion_order"]),
        concentrated_entropy(minus["bulk_counts_in_tag_then_word_insertion_order"]),
        step,
    )
    e = sqrt(5 / TRIALS)

    def reading(key, contrast):
        d = contrast * (plus[key] - minus[key]) / (TRIALS * step)
        err = 2 * abs(contrast) * e / step
        return [d - err, d + err]

    a0, a1 = [r[0] for r in packet["sector_area_volume"]]
    return {
        "entropy": s,
        "area": reading("rod_sector_hits", a1 - a0),
        "ordinary_energy": reading("ordinary_energy_hits", 1),
    }


def run():
    rng = np.random.Generator(np.random.PCG64(63019))
    out = {}
    for name, cut in CUTS.items():
        arm = "ordinary" if name == "ordinary_line" else "held"
        report = terminal_counts(cut, arm, rng)
        s, a, v, _ = projection_slopes(cut, arm)
        cs, _, _, _ = projection_slopes(cut, arm, biased=True)
        report.update(
            {
                "controller_visible": ORIGIN in cut,
                "measured_secant_99": read(report),
                "entropy_derivative": s,
                "area_derivative": str(a),
                "volume_derivative": str(v),
                "biased_reference_residual": cs + log(2) * float(a) / 2,
            }
        )
        out[name] = report
    so, ao = (
        out["ordinary_line"]["measured_secant_99"][k] for k in ("entropy", "area")
    )
    eta = [so[0] / (-ao[0]), so[1] / (-ao[1])]
    hs, ha = (out["held_sheet"]["measured_secant_99"][k] for k in ("entropy", "area"))
    return {
        "seed": 63019,
        "trials_per_bank": TRIALS,
        "ordinary_eta_99": eta,
        "held_residual_99": [hs[0] + eta[0] * ha[0], hs[1] + eta[1] * ha[1]],
        "regions": out,
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
