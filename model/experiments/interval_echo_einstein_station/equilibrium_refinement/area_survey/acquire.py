"""Independent source-word and complete three-outcome return acquisitions on 80 cells."""

import json
from math import log

import numpy as np
from scipy.optimize import brentq
from mpmath import iv

from . import measure

SHOTS = 2**55
BATCH = 2**49
STEP = 1 / 1024
TAIL = log(2 * 2000 / 0.01)
# Population error for the 80-cell conditional return, independently certified.
NATIVE_ERROR = 0.000000053


def count(rng, probs):
    probs = np.array(probs, dtype=float, copy=True)
    probs[-1] = 1 - sum(probs[:-1])
    return rng.multinomial(BATCH, probs, size=SHOTS // BATCH).sum(0)


def band(x):
    y = float(x)
    d = np.sqrt(2 * y * TAIL) + 2 * TAIL
    return iv.mpf([max(0.0, (y - d) / SHOTS - 1e-14), (y + d) / SHOTS + 1e-14])


def pair(z):
    return [
        float(np.nextafter(float(z.a), -np.inf)),
        float(np.nextafter(float(z.b), np.inf)),
    ]


def interval_geometry(edges):
    squared = np.empty((4, 4), dtype=object)
    squared[:] = iv.mpf(0)
    for (i, j), d in zip(measure.EDGES, edges, strict=True):
        squared[i, j] = squared[j, i] = d * d
    g = (squared[0, 1:, None] + squared[0, None, 1:] - squared[1:, 1:]) / 2
    a, b, c = g[0]
    _, d, e = g[1]
    _, _, f = g[2]
    volume = iv.sqrt(a * d * f + 2 * b * c * e - a * e**2 - d * c**2 - f * b**2) / 6
    area = iv.mpf(0)
    for i, j, k in measure.combinations(range(4), 3):
        a, b, c = squared[i, j], squared[i, k], squared[j, k]
        area += iv.sqrt(4 * a * b - (a + b - c) ** 2) / 4
    return area, volume


def acquisition(s, pop, labels, rng):
    counts = count(rng, np.bincount(labels, weights=pop))
    entropy = iv.mpf(0)
    for n in counts:
        b = band(n)
        entropy -= b * iv.log(b)
    enrolled = pop @ s.starts
    success = enrolled - pop @ s.scores / (
        measure.MEAN_TIMEOUT * measure.source.EXPOSURE / 2
    )
    events = []
    ranges = []
    for h, a in zip(enrolled, success, strict=True):
        row = count(rng, [a, h - a, 1 - h])
        events.append(row.tolist())
        a = band(int(row[0]))
        h = band(int(row[0] + row[1]))
        ranges.append(
            (1 - a / h + iv.mpf([-1, 1]) * iv.mpf(str(NATIVE_ERROR)))
            * (measure.MEAN_TIMEOUT * measure.source.EXPOSURE / 2)
        )
    edges = [
        (ranges[measure.RAYS.index((i, j))] + ranges[measure.RAYS.index((j, i))]) / 2
        for i, j in measure.EDGES
    ]
    area, volume = interval_geometry(edges)
    return dict(interior_counts=counts.tolist(), return_counts=events), [
        entropy,
        area,
        volume,
    ]


def held_lottery(s, labels, rng):
    """A dyadic approximation to the analytically matched finite source arm."""
    base = s.uniform
    returns = s.starts - s.scores / (measure.MEAN_TIMEOUT * measure.source.EXPOSURE / 2)
    design = np.vstack(
        [
            np.ones(len(base)),
            s.h.diagonal() / 256,
            s.starts.T,
            returns.T,
            s.cut_count == 2,
        ]
    )
    target = np.zeros(len(design))
    target[-1] = 1
    direction = np.linalg.lstsq(design, target, rcond=1e-11)[0]
    tickets = np.rint((base - STEP * direction) * 2**50).astype(np.int64)
    tickets[tickets.argmax()] += 2**50 - tickets.sum()
    assert tickets.min() > 0
    loss = tickets / 2**50
    change = loss - base
    ref, xr = acquisition(s, base, labels, rng)
    arm, yr = acquisition(s, loss, labels, rng)
    native_offset = (
        np.max(abs(change @ returns)) + NATIVE_ERROR * np.linalg.norm(change, 1) + 1e-13
    )
    return dict(
        source_word_tickets=tickets.tolist(),
        total_tickets=2**50,
        baseline=ref,
        held=arm,
        entropy_difference=pair(yr[0] - xr[0]),
        area_difference=pair(yr[1] - xr[1]),
        volume_difference=pair(yr[2] - xr[2]),
        max_native_success_difference_upper=float(native_offset),
        largest_herald_probability_difference=float(np.max(abs(change @ s.starts))),
        measured_energy_difference=float(change @ s.h.diagonal()),
    )


def run():
    iv.dps = 35
    s = measure.survey(1)
    p = s.uniform
    grad, _ = measure.differentials(p, s.starts, s.scores)
    directions, ctrl, words = measure.survey_directions(s, grad)
    groups = [tuple(x for x in w if x not in s.rim) for w in s.words]
    keys = {word: i for i, word in enumerate(sorted(set(groups)))}
    labels = np.array([keys[w] for w in groups])

    def volume(v):
        return measure.geometry(measure.distances(v, s.starts, s.scores))[0][1]

    target = volume(p)
    rng = np.random.default_rng(763001)
    seen = {}
    obs = {}
    for key in ("zero_to_one", "zero_to_two", "matched_energy_and_cut_load"):
        raw = directions[key]
        ends = []
        hist = {}
        for sign in (-1, 1):
            t = sign * STEP

            def err(c):
                return volume(p + t * raw + c * ctrl) - target

            c = brentq(err, -0.001, 0.001, xtol=1e-15)
            pop = p + t * raw + c * ctrl
            assert pop.min() > 0 and abs(pop.sum() - 1) < 1e-12
            bank, readings = acquisition(s, pop, labels, rng)
            hist[str(sign)] = dict(
                compensation=c, computed_volume_offset=float(err(c)), **bank
            )
            ends.append(readings)
        delta = [(b - a) / (2 * STEP) for a, b in zip(*ends, strict=True)]
        if not delta[2].a <= 0 <= delta[2].b:
            raise ValueError("held clock-volume intervals miss the matched target")
        seen[key] = hist
        obs[key] = delta
    eta = -obs["zero_to_one"][0] / obs["zero_to_one"][1]
    held = {
        k: dict(
            entropy=pair(s), area=pair(a), volume=pair(v), residual=pair(s + eta * a)
        )
        for k, (s, a, v) in obs.items()
    }
    return dict(
        starts_per_setting=SHOTS,
        batch_size=BATCH,
        source_probability_step=STEP,
        confidence=0.99,
        enrolled_words=len(p),
        interior_records=len(keys),
        ordinary_energy_calibrating_transfer="zero_to_one",
        volume_control_words=words,
        calibrated_eta=pair(eta),
        held=held,
        counts=seen,
        held_lottery=held_lottery(s, labels, rng),
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
