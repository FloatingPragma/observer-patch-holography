"""Predict terminal charge-probe counts using only a held neutral field reading."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import ceil, exp, log, sqrt

import numpy as np

from . import measure as m

HERE = m.HERE
NFIELD = 2**38
NCHARGE = 2**32
ORIGINS = (10, 11)
RECEIVER = [2, 3, 4]
WEIGHTS = (F(1, 4), F(3, 4))
NBINS = 8 + 2 * 48 + 2 * 48 + 48


def count_error(n):
    k = ceil(sqrt(n * log(2 * NBINS / 0.01) / 2))
    assert 2 * NBINS * exp(-2 * k * k / n) < 0.01
    return F(k, n)


def generate(parent):
    rng = np.random.default_rng(320817)
    faces, probs, _ = m.laws()
    f = [list(face) for face in faces].index(RECEIVER)
    pure_field = probs[-1, f, :4].reshape(2, 2)
    h = parent["selected"]["head_index"]
    pure_ray = np.asarray(parent["phase_flights"])[:, h, :2][:, :, ORIGINS, :]
    assert pure_ray.shape == (2, 2, 2, 12)

    def mix(pure, weight):
        return (1 - float(weight)) * pure[0] + float(weight) * pure[1]

    field = [
        rng.binomial(NFIELD, mix(pure_field, w)).tolist() for w in (0, 1, *WEIGHTS)
    ]
    rays = []
    for w in WEIGHTS:
        ps = mix(pure_ray.swapaxes(0, 1), w)
        rays.append(
            [
                [rng.multinomial(NCHARGE, p / p.sum()).tolist() for p in sign]
                for sign in ps
            ]
        )
    return dict(
        neutral_clicks=field,
        charged_endpoints=rays,
        neutral_shots=NFIELD,
        charged_shots=NCHARGE,
        seed=320817,
        protocol_sha256=sha256(
            (HERE / "TRANSFER_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        original_calibration_sha256=sha256(
            (HERE.parent / "results.json").read_bytes()
        ).hexdigest(),
    )


def detector_weight(neutral, index):
    offsets = [F(v[1] - v[0], NFIELD) for v in neutral]
    d = offsets[1] - offsets[0]
    # The two endpoint differences and the held difference are each within 2e.
    assert d > 4 * count_error(NFIELD)
    n = offsets[index] - offsets[0]
    err = 4 * count_error(NFIELD)
    lo, hi = max(F(0), (n - err) / d), min(F(1), (n + err) / d)
    assert lo <= hi
    return (lo, hi)


def assess(cal, held, bounds):
    assert cal.shape == (2, 2, 2, 12) and held.shape == (2, 2, 12)
    assert (cal.sum(axis=-1) == NCHARGE).all() and (held.sum(axis=-1) == NCHARGE).all()
    error = NCHARGE * count_error(NCHARGE)
    misses, max_gap = 0, F(0)
    for sign, origin, dest in np.ndindex(held.shape):
        y0, y1 = (int(a) for a in cal[sign, :, origin, dest])
        pred = sorted(y0 + p * (y1 - y0) for p in bounds)
        y = int(held[sign, origin, dest])
        gap = max(F(0), pred[0] - 2 * error - y, y - pred[1] - 2 * error)
        misses += gap > 0
        max_gap = max(max_gap, gap)
    return dict(
        bin_failures=misses, worst_excluded_probability=float(max_gap / NCHARGE)
    )


def infer(parent, bank):
    assert bank["neutral_shots"] == NFIELD and bank["charged_shots"] == NCHARGE
    assert (
        bank["protocol_sha256"]
        == sha256((HERE / "TRANSFER_PROTOCOL.md").read_bytes()).hexdigest()
    )
    assert (
        bank["original_calibration_sha256"]
        == sha256((HERE.parent / "results.json").read_bytes()).hexdigest()
    )
    charged = np.array(parent["selected"]["held_source_counts"], dtype=np.int64)[
        :, :, ORIGINS, :
    ]
    assert charged.shape == (2, 3, 2, 12) and parent["shots"] == NCHARGE
    cal = charged[:, :2]
    reports = []
    for j, cohort in enumerate(bank["charged_endpoints"]):
        bounds = detector_weight(bank["neutral_clicks"], j + 2)
        held = np.array(cohort, dtype=np.int64)
        reports.append(
            dict(
                field_inferred_weight=[float(p) for p in bounds],
                ordinary_assignment=assess(cal, held, bounds),
                swapped_assignment=assess(cal[:, ::-1], held, bounds),
            )
        )
    # Vacuum has exactly zero signed neutral response. Find every mixture
    # consistent with that same difference, then challenge it on the vacuum
    # charge data. This separate test uses 148 marginals, below NBINS.
    a, b = [F(plus - minus, NFIELD) for minus, plus in bank["neutral_clicks"][:2]]
    error = 2 * count_error(NFIELD)
    assert a + error < 0 < b - error
    same_neutral = ((-a - error) / (b - a), (-a + error) / (b - a))
    quiet = np.array(parent["selected"]["held_vacuum_deadline_counts"], dtype=np.int64)
    return dict(
        total_checked_marginals=NBINS,
        pure_field_counts=bank["neutral_clicks"][:2],
        field_probability_radius=float(count_error(NFIELD)),
        charged_probability_radius=float(count_error(NCHARGE)),
        held=reports,
        archived_equal_mixture=assess(cal, charged[:, 2], (F(1, 2), F(1, 2))),
        vacuum_with_same_neutral_difference=dict(
            supplementary_marginals=148,
            matching_weight=[float(p) for p in same_neutral],
            charged_challenge=assess(cal, quiet[:, ORIGINS, :], same_neutral),
        ),
    )


def run(write=False):
    parent = json.loads((HERE.parent / "results.json").read_text())
    path = HERE / "transfer_counts.json"
    if write:
        path.write_text(json.dumps(generate(parent), indent=1) + "\n")
    return infer(parent, json.loads(path.read_text()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    print(json.dumps(run(parser.parse_args().write), indent=2))
