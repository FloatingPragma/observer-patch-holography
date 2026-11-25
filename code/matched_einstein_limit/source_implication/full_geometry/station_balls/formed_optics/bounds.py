"""Outward interval propagation of complete signed camera responses and their invariants."""

import json

from mpmath import iv
import numpy as np
from scipy.stats import beta

from . import measure


def ranges(counts):
    iv.dps = 34
    rows = np.asarray(counts)
    obs = rows[:, :, :5].transpose(0, 2, 1)
    lo = np.nan_to_num(beta.ppf(measure.ALPHA, obs, measure.SHOTS - obs + 1))
    hi = np.nan_to_num(
        beta.ppf(1 - measure.ALPHA, obs + 1, measure.SHOTS - obs), nan=1.0
    )
    # Leave room for floating point evaluation of tail endpoints.
    return np.array(
        [
            iv.mpf([max(0.0, float(low) - 1e-14), float(h) + 1e-14])
            for low, h in zip(lo.flat, hi.flat, strict=True)
        ],
        dtype=object,
    ).reshape(obs.shape)


def quotient(x, y):
    return None if y.a <= 0 <= y.b else x / y


def pair(v):
    if v is None:
        return None
    return [np.nextafter(float(v.a), -np.inf), np.nextafter(float(v.b), np.inf)]


def area(a, b):
    return a[0, 0] * b[1, 1] + a[1, 1] * b[0, 0] - a[1, 0] * b[0, 1] - a[0, 1] * b[1, 0]


def determinant(x):
    return x[0, 0] * x[1, 1] - x[1, 0] * x[0, 1]


def analyze(tables, *, formed=False):
    iv.dps = 34
    independent, h, full, ordinary = tables[:4]
    if formed:
        p = full
        change = 3 * full - 4 * h + independent
        o = tables[4]
    else:
        p = independent
        change = 4 * h - full - 3 * independent
        o = ordinary
    w = p.sum()
    if w.a <= 0:
        raise ValueError("no resolved image baseline")
    b = np.asarray(
        [
            [iv.mpf(int(v)) / measure.camera.UNIT for v in row]
            for row in measure.camera.PIXELS
        ],
        dtype=object,
    )
    # Subtract fixed center/scale pivots before interval evaluation.
    pm = np.vectorize(lambda a: float(a.mid))(p)
    wm = pm.sum()
    x = b - np.asarray(pm.sum(0) @ measure.BEARINGS / wm, dtype=object)
    y = b - np.asarray(pm.sum(1) @ measure.BEARINGS / wm, dtype=object)
    xs = x - p.sum(0) @ x / w
    ys = y - p.sum(1) @ y / w
    cov = ys.T @ p @ xs / w

    def derivative(dp):
        return (ys.T @ dp @ xs - cov * dp.sum()) * (2 / w)

    a = derivative(change)
    ordinary = derivative(2 * (o - p))
    ad = determinant(a)
    od = determinant(ordinary)
    af = area(2 * cov, a)
    ao = area(2 * cov, ordinary)
    mixed = area(a, ordinary)
    return dict(
        formed_area_slope=pair(af),
        ordinary_area_slope=pair(ao),
        energy_from_area=pair(quotient(af, ao)),
        relative_trace=pair(quotient(mixed, od)),
        relative_determinant=pair(quotient(ad, od)),
        ordinary_det=pair(od),
        formed_det=pair(ad),
        pencil_discriminant=pair(mixed**2 - 4 * ad * od),
    )


def report():
    receipt = json.loads((measure.HERE / "results.json").read_text())
    return {
        charge: {
            g: {s: analyze(ranges(bank["counts"])) for s, bank in dic.items()}
            for g, dic in reports.items()
        }
        for charge, reports in receipt["charges"].items()
    }


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
