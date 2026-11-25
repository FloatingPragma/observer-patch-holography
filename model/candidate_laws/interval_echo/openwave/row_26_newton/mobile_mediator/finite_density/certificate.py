"""Certified free energies on the finite vacancy medium, including near cancellation."""

import argparse
from itertools import product
import json
from math import ceil, floor

from mpmath import iv
import numpy as np
from scipy.linalg import eigh

from . import measure

SCALE = 10**10
EIGEN_ERROR = "0.00000001"


def free_energy(kappa, distance, a, b, cut=None, number=None):
    total = iv.mpf(0)
    error = iv.mpf([-iv.mpf(EIGEN_ERROR).b, iv.mpf(EIGEN_ERROR).b])
    for n, (h, bits) in enumerate(measure.sectors(kappa, cut)):
        if number is not None and n != number:
            continue
        h = h + np.diag(a * bits[:, 0] + b * bits[:, distance])
        e, v = eigh(h, driver="evd")
        # Each input and source field is an exact integer or half-integer.
        assert np.array_equal(2 * h, (2 * h).astype(int))
        assert len(h) <= 126 and max(abs(h).sum(axis=1)) < 511
        assert max(abs(e)) < 511 and np.linalg.norm(v) < 11.9
        assert np.linalg.norm(h @ v - v * e) < 5e-11
        assert np.linalg.norm(v.T @ v - np.eye(len(e))) < 5e-13
        total += sum(iv.exp(-(iv.mpf(float(x)) + error) / 16) for x in e)
    return -16 * iv.ln(total)


def interval(x):
    # Endpoints in this experiment have magnitude <2 in source energy units.
    return [floor(float(x.a) * SCALE) - 1, ceil(float(x.b) * SCALE) + 1]


def panels(kappa, d, cut=None, number=None):
    v = measure.mediator.contact(kappa)[0] / 2
    arms = {
        (a, b): free_energy(kappa, d, a, b, cut, number)
        for a, b in product((0, v[0], v[2]), repeat=2)
    }
    values = [
        arms[a, b] + arms[0, 0] - arms[a, 0] - arms[0, b]
        for a, b in product((v[0], v[2]), repeat=2)
    ]
    return dict(signed=[interval(x) for x in values], mean=interval(sum(values) / 4))


def certify():
    rows = [panels(8, d) for d in measure.DISTANCES]
    # For a potential u-C/(d+c), C>0 and d+c>0, this curvature is negative.
    # Strictly positive directed lower endpoint rejects every origin shift.
    curvature = [
        rows[0]["mean"][i] + rows[2]["mean"][i] - 2 * rows[1]["mean"][1 - i]
        for i in (0, 1)
    ]
    assert curvature[0] > 0 and all(p["mean"][1] < 0 for p in rows)
    return dict(
        scale=SCALE,
        active=rows,
        second_difference=curvature,
        changed_law_first=panels(4, 1),
        cut_grand=panels(8, 3, 1),
        cut_fixed_number=panels(8, 3, 1, 4),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    report = certify()
    path = measure.HERE / "bounds.json"
    if args.write:
        path.write_text(json.dumps(report, indent=2) + "\n")
    else:
        old = json.loads(path.read_text())

        def compare(now, saved):
            if isinstance(now, dict):
                assert now.keys() == saved.keys()
                for k in now:
                    compare(now[k], saved[k])
            elif isinstance(now, list):
                assert len(now) == len(saved)
                if len(now) == 2 and isinstance(now[0], int):
                    assert max(now[0], saved[0]) <= min(now[1], saved[1])
                    assert max(now[1] - now[0], saved[1] - saved[0]) <= 4000
                else:
                    for x, y in zip(now, saved, strict=True):
                        compare(x, y)
            else:
                assert now == saved

        compare(report, old)
        print("Replayed the native all-sector connected free-energy enclosures.")
