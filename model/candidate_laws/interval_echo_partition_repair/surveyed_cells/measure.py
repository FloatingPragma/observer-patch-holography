"""Cubical rods, bulk face registers and an ordinary material excitation."""

import argparse
import hashlib
import json
from collections import Counter
from decimal import Decimal, localcontext
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
SEED = 63019
BULK = "bulk_counts_in_tag_then_word_insertion_order"
# SHA-256 of only the six integer banks in sorted, compact ASCII JSON.
FROZEN_BANKS_SHA256 = "e8a9ddb4f80a9c8b412f57afb276b91824e6b88672a93f1b1aa34a5bd7e87eff"


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
    with localcontext() as ctx:
        ctx.prec = 50
        def decimal(q):
            return Decimal(q.numerator)/Decimal(q.denominator)
        entropy = -float(sum(decimal((hi[k]-lo[k])/(2*STEP))*decimal(p).ln()
                             for k,p in ref.items()))
    dg = sum(sector(k) * (a - b) for k, a, b in zip(TAGS, plus, minus, strict=True)) / (
        2 * STEP
    )
    geo = [rod_survey(positions(cut, g)) for g in (0, 1)]
    area = dg * (geo[1][0] - geo[0][0])
    volume = dg * (geo[1][1] - geo[0][1])
    return entropy, area, volume, dg


def terminal_probabilities(cut, arm, sign):
    weights = tag_weights(arm, sign * STEP)
    words = histogram(cut, weights)
    sector_probability = sum(p * sector(tag) for p, tag in zip(weights, TAGS, strict=True))
    energy_probability = sum(p * tag[0] for p, tag in zip(weights, TAGS, strict=True))
    return tuple(words.values()), sector_probability, energy_probability


def terminal_counts(cut, arm, rng):
    raw = {}
    for sign in (-1, 1):
        words, sector_probability, energy_probability = terminal_probabilities(cut, arm, sign)
        counts = rng.multinomial(TRIALS, list(map(float, words)))
        area_hits = rng.binomial(TRIALS, float(sector_probability))
        energy_hits = rng.binomial(TRIALS, float(energy_probability))
        raw[str(sign)] = {
            BULK: counts.tolist(),
            "rod_sector_hits": int(area_hits),
            "ordinary_energy_hits": int(energy_hits),
        }
    return {"sector_area_volume": [rod_survey(positions(cut, g)) for g in (0, 1)], "records": raw}


def concentrated_entropy(counts):
    # Bounded differences for the plug-in estimator, plus its standard
    # elementary bias bound 0<=H(p)-E H(phat)<=log(1+(k-1)/N).
    # Subtracting two bulk entropies and comparing the near-zero held
    # residual magnifies platform differences in vectorized log/dot.
    # Reduce the frozen integer counts with correctly rounded decimal logs.
    with localcontext() as ctx:
        ctx.prec = 50
        n = Decimal(sum(counts))
        probabilities = (Decimal(c)/n for c in counts if c)
        estimate = -sum(p*p.ln() for p in probabilities)
        radius = (n.ln()+1)*(5/n).sqrt()
        bias = (1+Decimal(len(counts)-1)/n).ln()
        return [float(estimate-radius), float(estimate+radius+bias)]


def verify_model_bank(cut, arm, sign, bank):
    """Test one observed integer bank against the law's exact probabilities."""
    if set(bank) != {BULK, "rod_sector_hits", "ordinary_energy_hits"}:
        raise ValueError("unexpected cubical detector keys")
    words, sector_probability, energy_probability = terminal_probabilities(cut, arm, sign)
    counts = bank[BULK]
    if not isinstance(counts, list) or len(counts) != len(words):
        raise ValueError("bulk counts disagree with the model alphabet")

    def check_hits(hits, probability, squared_radius_factor, channel):
        if type(hits) is not int or not 0 <= hits <= TRIALS:
            raise ValueError(f"{channel} hits must be integer counts in the bank")
        if (probability in (0, 1) and hits != TRIALS * probability) or (
            (hits - TRIALS * probability) ** 2 > squared_radius_factor * TRIALS
        ):
            raise ValueError(f"{channel} counts lie outside the model probability bound")

    for hits, probability in zip(counts, words, strict=True):
        check_hits(hits, probability, 8, "bulk word")
    if sum(counts) != TRIALS:
        raise ValueError("bulk counts do not sum to the trial total")
    check_hits(bank["rod_sector_hits"], sector_probability, 5, "rod")
    check_hits(bank["ordinary_energy_hits"], energy_probability, 5, "ordinary energy")
    with localcontext() as ctx:
        ctx.prec = 50
        probabilities = (Decimal(p.numerator) / Decimal(p.denominator) for p in words if p)
        population_entropy = float(-sum(p * p.ln() for p in probabilities))
    lo, hi = concentrated_entropy(counts)
    if not lo <= population_entropy <= hi:
        raise ValueError("bulk confidence interval excludes the model entropy")


def interval_subtract(a, b, scale):
    return [(a[0] - b[1]) / scale, (a[1] - b[0]) / scale]


def read(packet):
    minus, plus = (packet["records"][str(i)] for i in (-1, 1))
    step = float(2 * STEP)
    s = interval_subtract(
        concentrated_entropy(plus[BULK]),
        concentrated_entropy(minus[BULK]),
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


def analyze(banks, seed):
    out = {}
    for name, cut in CUTS.items():
        arm = "ordinary" if name == "ordinary_line" else "held"
        report = {"sector_area_volume": [rod_survey(positions(cut, g)) for g in (0, 1)], "records": banks[name]}
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
        "seed": seed,
        "trials_per_bank": TRIALS,
        "ordinary_eta_99": eta,
        "held_residual_99": [hs[0] + eta[0] * ha[0], hs[1] + eta[1] * ha[1]],
        "regions": out,
    }


def replay(saved):
    """Authenticate the original observations, then recompute the entire receipt."""
    settings = (saved["seed"], saved["trials_per_bank"])
    if any(type(value) is not int for value in settings) or settings != (SEED, TRIALS):
        raise ValueError("cubical frozen acquisition has unexpected seed or trial count")
    if set(saved["regions"]) != set(CUTS):
        raise ValueError("cubical frozen acquisition has unexpected settings or regions")
    banks = {name: row["records"] for name, row in saved["regions"].items()}
    encoded = json.dumps(banks, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("ascii")
    if hashlib.sha256(encoded).hexdigest() != FROZEN_BANKS_SHA256:
        raise ValueError("cubical integer bank SHA-256 differs from the frozen acquisition")
    for name, cut in CUTS.items():
        if set(banks[name]) != {"-1", "1"}:
            raise ValueError("cubical bank requires both signed settings")
        arm = "ordinary" if name == "ordinary_line" else "held"
        for sign in (-1, 1):
            verify_model_bank(cut, arm, sign, banks[name][str(sign)])
    return analyze(banks, SEED)


def run():
    return replay(json.loads((HERE / "results.json").read_text()))


def generate(seed=SEED):
    """Make an independent acquisition using the current NumPy sampling routines."""
    rng = np.random.Generator(np.random.PCG64(seed))
    banks = {}
    for name, cut in CUTS.items():
        arm = "ordinary" if name == "ordinary_line" else "held"
        banks[name] = terminal_counts(cut, arm, rng)["records"]
    return analyze(banks, seed)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", nargs="?", const=SEED, type=int, metavar="SEED",
                        help="draw fresh counts using PCG64 (default seed 63019)")
    args = parser.parse_args()
    print(json.dumps(run() if args.generate is None else generate(args.generate), indent=2, allow_nan=False))
