"""Bound the measured overlap relay against the SPARC outer mass requirement."""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from os.path import relpath

from mpmath import iv
import numpy as np
import sympy as sp

from . import reference_arrival as source
from .contact_lapse_sky import sky

HERE = source.local.HERE
# Conservative upper conversion in kpc (km/s)^2 / solar mass.
G_HIGH = F("0.00000431")


def source_caps():
    z = source.local.inputs("reference_arrival")
    ref = source.local.survey.standards(source.local.survey.bulk.federation(0))[0][0]
    groups = [np.flatnonzero(ref[0] == n) for n in (32, 48, 16)]
    p = sp.Matrix([sp.Rational(int(ref[ids, :].sum()), 2592) for ids in groups])
    matrix = sp.zeros(5, 9)
    for i, j in np.ndindex(3, 3):
        matrix[i, 3 * i + j] = 1
        if j < 2:
            matrix[3 + j, 3 * i + j] = 1
    rhs = p.col_join(p[:2, :])
    peaks = set()
    for cols in combinations(range(9), 5):
        m = matrix[:, cols]
        if m.det() == 0:
            continue
        v = m.inv() * rhs
        if min(v) < 0:
            continue
        out = [F(0)] * 9
        for j, a in zip(cols, v, strict=True):
            out[j] = F(a)
        peaks.add(tuple(out))
    vertices = sorted(peaks)
    h = iv.mpf(z["harmonic_normalization"])
    eps = iv.mpf([0, z["max_charge_error_nats"]])
    b = [
        (-iv.log(iv.mpf(int(ref[i[0], j[0]])) / 2592) - eps) / h
        for i in groups
        for j in groups
    ]
    values = [
        sum(sky.sample.box(v) * x for v, x in zip(row, b, strict=True))
        for row in vertices
    ]
    lo = min(float(x.a) for x in values)
    hi = max(float(x.b) for x in values)
    full = max(float(x.b) for x in b) - min(float(x.a) for x in b)
    # Bounds rounded outwards compared with all candidate vertices.
    return {
        "class_masses": [str(x) for x in p],
        "vertices": len(vertices),
        "min_vertex": [
            str(x)
            for x in vertices[min(range(len(values)), key=lambda i: float(values[i].a))]
        ],
        "max_vertex": [
            str(x)
            for x in vertices[max(range(len(values)), key=lambda i: float(values[i].b))]
        ],
        "fixed_marginal_span_upper": hi - lo + 1e-14,
        "any_population_span_upper": full + 1e-14,
        "enclosed_rest_stock_upper": 1,
    }


def run():
    iv.dps = 45
    caps = source_caps()
    frozen = sky.HERE / "outer_requirement.json"
    meta_path, rot_path = [
        sky.sample.measure.DATA / n for n in ("table1.dat", "table2.dat")
    ]
    meta = {
        r[:11].strip(): (i, r)
        for i, r in enumerate(meta_path.read_text().splitlines(), 1)
    }
    data = rot_path.read_text().splitlines()
    surveys = {}
    for name, entry in json.loads(frozen.read_bytes())["galaxies"].items():
        num, t = meta[name]
        _, hmax = sky.sample.slope_box(t)
        hmax = F(sky.sample.upper(hmax))
        lum, lum_err, hi = [F(t[a:b]) for a, b in ((40, 47), (48, 55), (86, 93))]
        rings = sky.rings_in(data, entry["rings"])
        result = {}
        for label, gas_factors in (
            ("catalogue", (F(1),)),
            ("expanded_HI", (F(1, 2), F(2))),
        ):
            # Maximize the enclosed baryonic velocity and the total baryonic
            # stock separately to favor an ordinary-material interpretation.
            stock = (
                F("1.4") * (lum + 2 * lum_err) + F("1.33") * max(gas_factors) * hi
            ) * 10**9
            assert stock > 0
            lower = []
            for radius, slow, _, (gas, disk, bulge), line in rings:
                baryon = max(
                    g * gas + d * disk + b * bulge
                    for g in gas_factors
                    for d, b in sky.sample.RECT
                )
                mass = radius * max(0, slow / hmax - baryon) / (G_HIGH * stock)
                lower.append((mass, line))
            worst, line = max(lower)
            result[label] = {
                "required_fraction_lower": str(worst),
                "witness_line": line,
                "stock_upper_at_catalogue_distance_msun": str(stock),
                "excludes": [
                    k
                    for k in (
                        "fixed_marginal_span_upper",
                        "any_population_span_upper",
                        "enclosed_rest_stock_upper",
                    )
                    if worst > F(caps[k])
                ],
            }
        surveys[name] = {"table1_line": num, "ring_count": len(rings), "bounds": result}
    files = [
        HERE / "ARCHIVED_OVERLAP.md",
        HERE / "reference_arrival.json",
        frozen,
        meta_path,
        rot_path,
    ]
    summary = {
        m: {
            k: sum(k in row["bounds"][m]["excludes"] for row in surveys.values())
            for k in (
                "fixed_marginal_span_upper",
                "any_population_span_upper",
                "enclosed_rest_stock_upper",
            )
        }
        for m in ("catalogue", "expanded_HI")
    }
    return {
        "source_caps": caps,
        "galaxies": surveys,
        "rings": sum(r["ring_count"] for r in surveys.values()),
        "summary": summary,
        "inputs": {relpath(q, HERE): sha256(q.read_bytes()).hexdigest() for q in files},
    }


if __name__ == "__main__":
    z = run()
    source.coherent.steady.exp.emit("modular_stock_sky", z)
    print(z["source_caps"], z["summary"])
