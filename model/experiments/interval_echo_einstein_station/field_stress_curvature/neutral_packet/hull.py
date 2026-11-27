"""Identify all three prepared fractions from two neutral field probabilities."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import exp, isqrt
import numpy as np

from . import measure as m, transfer as t

HERE = m.HERE
HALF_PLANES = ((-1, 0, 0), (0, -1, 0), (1, 1, 1))
MARGINALS = 640 + 4 + 6 + 144 + 144


def band(k, n):
    """Invert a Bernstein event with true binomial mean bounded by k+error."""
    k, n = int(k), int(n)
    assert 0 <= k <= n and 2 * MARGINALS * exp(-15) < 0.01
    d = 21 + isqrt(400 + 30 * k)
    return F(max(0, k - d), n), F(min(n, k + d), n)


def volumes(refs):
    # The determinant is affine in each of these six probabilities separately.
    vals = []
    for a, b, c, d, e, f in product(*(edge for state in refs for edge in state)):
        vals.append((a - e) * (d - f) - (b - f) * (c - e))
    return min(vals), max(vals)


def vertices(refs, held):
    restrictions = list(HALF_PLANES)
    for j, (lo, hi) in enumerate(held):
        (a, b, v) = [r[j][0] for r in refs]
        restrictions.append((a - v, b - v, hi - v))
        (a, b, v) = [r[j][1] for r in refs]
        restrictions.append((v - a, v - b, v - lo))
    corners = set()
    for (a, b, c), (d, e, f) in combinations(restrictions, 2):
        det = a * e - b * d
        if not det:
            continue
        x, y = F(c * e - b * f) / det, F(a * f - c * d) / det
        if all(a * x + b * y <= c for a, b, c in restrictions):
            corners.add((x, y, 1 - x - y))
    assert corners
    return tuple(sorted(corners))


def read_banks():
    parent = json.loads((HERE.parent / "results.json").read_text())
    old = json.loads((HERE / "transfer_counts.json").read_text())
    first = json.loads((HERE / "counts.json").read_text())
    face = first["faces"].index(t.RECEIVER)
    refs = [[band(k, t.NFIELD) for k in arm] for arm in old["neutral_clicks"][:2]]
    refs.append([band(k, m.SHOTS) for k in first["clicks"][-1][face][4:6]])
    cal = np.array(parent["selected"]["held_source_counts"], dtype=np.int64)[
        :, :2, t.ORIGINS, :
    ].transpose(1, 0, 2, 3)
    vac = np.array(parent["selected"]["held_vacuum_deadline_counts"], dtype=np.int64)[
        :, t.ORIGINS, :
    ]
    return refs, np.concatenate((cal, vac[None, :, :, :])), old, parent


def acquire(parent):
    weights = (0.25, 0.125, 0.625)
    seed = 710493
    rng = np.random.default_rng(seed)
    faces, p, _ = m.laws()
    face = [list(f) for f in faces].index(t.RECEIVER)
    neutral = p[-1, face, :6].reshape(3, 2)
    phase = np.array(parent["phase_flights"])[:, parent["selected"]["head_index"], :2][
        :, :, t.ORIGINS, :
    ].transpose(1, 0, 2, 3)
    quiet = np.array(parent["vacuum_deadline_kernels"])[:, t.ORIGINS, :]
    ray = np.concatenate((phase, quiet[None, :, :, :]))
    field = np.einsum("s,si->i", weights, neutral)
    charged = np.einsum("s,sijk->ijk", weights, ray)
    return dict(
        seed=seed,
        protocol_sha256=sha256((HERE / "HULL_PROTOCOL.md").read_bytes()).hexdigest(),
        neutral_shots=t.NFIELD,
        ray_shots=t.NCHARGE,
        neutral=rng.binomial(t.NFIELD, field).tolist(),
        charged=[
            [rng.multinomial(t.NCHARGE, row / row.sum()).tolist() for row in sign]
            for sign in charged
        ],
    )


def assess(corners, cal, held):
    assert cal.shape == (3, 2, 2, 12) and held.shape == (2, 2, 12)
    assert (held.sum(axis=-1) == t.NCHARGE).all() and (
        cal.sum(axis=-1) == t.NCHARGE
    ).all()
    errors = []
    for idx in np.ndindex(held.shape):
        intervals = [band(arm[idx], t.NCHARGE) for arm in cal]
        lo = min(sum(w[i] * intervals[i][0] for i in range(3)) for w in corners)
        hi = max(sum(w[i] * intervals[i][1] for i in range(3)) for w in corners)
        lh, uh = band(held[idx], t.NCHARGE)
        errors.append(max(0, lo - uh, lh - hi))
    return dict(
        excluded=sum(e > 0 for e in errors), largest_exclusion=float(max(errors))
    )


def run(write=False):
    refs, cal, old, parent = read_banks()
    path = HERE / "hull_counts.json"
    if write:
        path.write_text(json.dumps(acquire(parent), indent=1) + "\n")
    new = json.loads(path.read_text())
    assert (
        new["protocol_sha256"]
        == sha256((HERE / "HULL_PROTOCOL.md").read_bytes()).hexdigest()
    )
    assert new["neutral_shots"] == t.NFIELD and new["ray_shots"] == t.NCHARGE
    det = volumes(refs)
    assert det[0] > 0
    # Magnitudes of the two columns of the two-rate calibration Jacobian.
    entries = [
        [
            max(abs(r[j][1] - refs[2][j][0]), abs(r[j][0] - refs[2][j][1]))
            for j in range(2)
        ]
        for r in refs[:2]
    ]
    inverse = max(map(sum, entries)) / det[0]
    condition = inverse * max(map(sum, zip(*entries)))
    reference_band = [[list(map(float, b)) for b in arm] for arm in refs]
    reports = []
    for field, charged in zip(
        old["neutral_clicks"][2:] + [new["neutral"]],
        old["charged_endpoints"] + [new["charged"]],
        strict=True,
    ):
        corners = vertices(refs, [band(k, t.NFIELD) for k in field])
        held = np.array(charged, dtype=np.int64)
        reports.append(
            dict(
                weights=[[float(min(col)), float(max(col))] for col in zip(*corners)],
                charged=assess(corners, cal, held),
                cyclic=assess(corners, cal[[1, 2, 0]], held),
            )
        )
    return dict(
        weight_order=["negative_phase", "positive_phase", "vacuum"],
        determinant=list(map(float, det)),
        infinity_inverse_norm_upper=float(inverse),
        infinity_condition_upper=float(condition),
        calibrated_neutral_rates=reference_band,
        coverage_marginals=MARGINALS,
        held=reports,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    print(json.dumps(run(parser.parse_args().write), indent=2))
