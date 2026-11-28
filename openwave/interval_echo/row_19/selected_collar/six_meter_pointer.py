"""Certify a common six-meter ordinary response and test its held clock geometry."""

import argparse
from hashlib import sha256
import json

import numpy as np
import sympy as sp

from . import charge_selected_flight as selected

mat, reg, lab, iv = selected.material, selected.reg, selected.lab, selected.iv
HERE = selected.HERE
SEED = 28112118


def components(rows):
    low = sp.Matrix([[sp.Rational(z[0]) for z in row] for row in rows])
    high = sp.Matrix([[sp.Rational(z[1]) for z in row] for row in rows])
    return (low + high) / 2, (high - low) / 2


def positive(matrix):
    return matrix.applyfunc(abs)


def bounds(center, error):
    return [
        float(np.nextafter(float(center - error), -np.inf)),
        float(np.nextafter(float(center + error), np.inf)),
    ]


def run():
    iv.dps = 35

    def load_record(n):
        return json.loads((HERE / (n + ".json")).read_bytes())

    prepared = load_record("charge_selected_flight")
    z = prepared["windows"]
    times = list(z)
    ordinary = load_record("ordinary_service_distribution")["windows"]
    source = load_record("source_interference")["windows"]
    costs = load_record("density_material")["windows"]
    labels = list(z[times[0]]["material"])

    def row(name):
        return [z[t]["meter_quadratures"][name][j] for t in times for j in (0, 1)]

    A, E = components([row(k) for k in labels])
    K = A.inv(method="DM")
    H = positive(K) * E
    G = (sp.eye(6) - H).inv(method="DM")
    assert all(v >= 0 for v in G)
    scaling = G * sp.ones(6, 1)
    assert all(v >= 1 for v in scaling)
    assert H * scaling == scaling - sp.ones(6, 1)
    target = [
        components([[ordinary[t]["held_ordinary"][k]] for k in labels]) for t in times
    ]
    y = sum((a for a, _ in target), sp.zeros(6, 1))
    ey = sum((b for _, b in target), sp.zeros(6, 1))
    x = K * y
    radius = G * (positive(K) * ey + H * positive(x))
    maxabs = positive(x) + radius
    beta = sp.Rational(1)
    while any(
        beta**2 * (maxabs[j] ** 2 + maxabs[j + 1] ** 2) > sp.Rational(9, 16)
        for j in (0, 2, 4)
    ):
        beta /= 2
    norm_factor = beta / 3

    # The same solution occurs in every held observable. Propagate the
    # ordinary-target and response-matrix errors before projecting them.
    def projected(prices):
        g, eg = components([prices])
        e = (positive(g * K) * (ey + E * maxabs) + eg * maxabs)[0]
        return bounds(norm_factor * (g * x)[0], norm_factor * e)

    charge = projected(
        [b for t in times for b in source[t]["measured_price_quadratures"]["connected"]]
    )
    price_sum = iv.mpf(0)
    for t in times:
        price_sum += iv.mpf(costs[t]["ordinary_energy"])
    mass = lab.six.bounds(float(beta) * price_sum / 3)
    knobs = {
        t: [bounds(beta * x[j], beta * radius[j]) for j in (i * 2, i * 2 + 1)]
        for i, t in enumerate(times)
    }
    centers = {
        t: complex(float(beta * x[2 * i]), float(beta * x[2 * i + 1]))
        for i, t in enumerate(times)
    }
    assert max(abs(v) for v in centers.values()) <= 0.75
    rng = np.random.default_rng(SEED)
    _, supplied, *_ = mat.c.trials()
    X = supplied["reference"]
    ports, bins = lab.six.survey.return_bins()
    clock = 3 * np.eye(4) + np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    routes = lab.six.survey.ROUTES
    names = [f"clock/{a}/{b}" for a, b in routes]
    anchors = [i for i, r in enumerate(routes) if max(r) < 4]
    ref = reg.src.observed(prepared["source_reference_heads"], lab.SHOTS, False)
    bases = []
    for t in times:
        a, b = (
            reg.src.observed(z[t]["phase_standard_heads"][k], lab.SHOTS, False)
            for k in ("plus", "minus")
        )
        bases.append(
            tuple(
                (1 - mat.density.DOSE) * u + mat.density.DOSE * (v + w) / 2
                for u, v, w in zip(ref, a, b, strict=True)
            )
        )
    base = tuple(sum(m[j] for m in bases) / 3 for j in (0, 1))
    start = [selected.fringe.read(R, base, False) for R in clock]
    dr = {n: iv.mpf(projected(row(n))) for n in names}
    ordinary_delta = {
        n: float(beta) * sum(iv.mpf(ordinary[t]["held_ordinary"][n]) for t in times) / 3
        for n in names
    }
    predicted_geometry = selected.geometry_change(start, list(dr.values()), anchors)
    ordinary_geometry = selected.geometry_change(
        start, list(ordinary_delta.values()), anchors
    )
    # Acquire one mixed source, rather than averaging the nonlinear areas.
    mixed, _, active = mat.density.history(
        X, ports, rng, {"held_" + t: centers[t] for t in times}
    )
    states = {}
    active_states = {}
    for key in ("held", "dense"):
        arms = ["held_" + t if key == "held" else "dense" for t in times]
        states[key] = (
            sum(mixed[int(t)][a] for t, a in zip(times, arms, strict=True)) / 3
        )
        active_states[key] = (
            sum(active[int(t)][a].sum(0) for t, a in zip(times, arms, strict=True)) / 3
        )
    heads = {k: lab.pairs(z, rng, False) for k, z in active_states.items()}
    observed = {k: reg.src.observed(v, lab.SHOTS, False) for k, v in heads.items()}
    q = mat.density.DOSE * mat.density.price_contrast(
        ref, observed["held"], observed["dense"], "connected"
    )
    radar = {
        k: [lab.counts(np.r_[p, 1 - p.sum()], rng) for p in mat.c.station.pair(bins, z)]
        for k, z in states.items()
    }
    means, _ = reg.timing_slopes(radar, lab.SHOTS, (("held", "dense", 1),))
    direct_geometry = selected.geometry_change(
        means["dense"],
        [v - w for v, w in zip(means["held"], means["dense"], strict=True)],
        anchors,
    )
    effects, scales, hits = mat.beam.flight(mat.METERS, include_edges=False)
    material = {}
    for key, Eout in effects.items():
        receipts = {}
        measured = {}
        for arm, src in active_states.items():
            receipts[arm], measured[arm] = mat.acquire(
                rng, 2**72, src, Eout, hits[key], scales[key]
            )
        n = f"{key[1]}/{key[0]}"
        delta = mat.density.DOSE * (measured["held"] - measured["dense"])
        c, ce = components([row(n)])
        center = bounds(
            norm_factor * (c * x)[0],
            norm_factor * (ce * positive(x))[0] + sp.Rational("1e-17"),
        )
        material[n] = {
            "counts": receipts,
            "observed": lab.six.bounds(delta),
            "forecast_at_center": center,
        }
    return {
        "seed": SEED,
        "times": times,
        "labels": labels,
        "common_ratio": [int(beta.p), int(beta.q)],
        "common_positive_energy": mass,
        "preconditioned_contraction_upper": float(
            np.nextafter(float(max((v - 1) / v for v in scaling)), np.inf)
        ),
        "exact_phase_enclosures": knobs,
        "central_settings": {t: [p.real, p.imag] for t, p in centers.items()},
        "connected_charge_exact": charge,
        "connected_charge_programmed_read": lab.six.bounds(q),
        "source_heads": heads,
        "geometry_exact": predicted_geometry,
        "equal_ordinary_geometry": ordinary_geometry,
        "geometry_programmed_read": direct_geometry,
        "clock_minus_equal_ordinary": {
            n: lab.six.bounds(dr[n] - ordinary_delta[n]) for n in names
        },
        "clock_counts": radar,
        "material_at_center": material,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "SIX_METER_POINTER.md",
                "charge_selected_flight.json",
                "source_interference.json",
                "ordinary_service_distribution.json",
                "density_material.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "six_meter_pointer.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    print(
        r["common_ratio"],
        r["common_positive_energy"],
        r["connected_charge_exact"],
        r["geometry_exact"],
        r["equal_ordinary_geometry"],
    )
