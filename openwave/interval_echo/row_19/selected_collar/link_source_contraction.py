"""Decode a finite entropy and area bound from ten terminal transfers and an independent volume."""

import argparse
import json
from fractions import Fraction
from hashlib import sha256
from os.path import relpath
from pathlib import Path

import numpy as np
from matched_einstein_limit.source_implication.full_algebra.enclosure import (
    geometry_jet,
)

from . import finite_dose_material as beam

c, reg, iv = beam.curve, beam.reg, beam.iv
bounds = c.lab.six.bounds


def clip(vertices, a, b, h):
    """Intersect a convex polygon with a*t+b*u <= h using exact dyadic inputs."""
    a, b, h = map(Fraction, (a, b, h))
    out = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1], strict=True):
        fs, fe = (a * t + b * u - h for t, u in (start, end))
        if fs <= 0:
            out.append(start)
        if (fs < 0 < fe) or (fe < 0 < fs):
            w = fs / (fs - fe)
            out.append(tuple(s + w * (e - s) for s, e in zip(start, end, strict=True)))
    return out


def infer(links, gains):
    polygon = [
        (Fraction(0), Fraction(0)),
        (Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(1)),
    ]
    for r in links:
        p, w, y = (r[k] for k in ("independent", "compensator", "held"))
        polygon = clip(polygon, p[0], w[0], y[1])
        polygon = clip(polygon, -p[1], -w[1], -y[0])
    if not polygon:
        raise ValueError("The measured transfers exclude the promised family")

    def outward(x, direction):
        return float(np.nextafter(float(x), direction))

    entropy = min(
        sum(Fraction(g) * v for g, v in zip(gains, vertex, strict=True))
        for vertex in polygon
    )
    return {
        "entropy_lower": outward(entropy, -np.inf),
        "mixture_intervals": [
            [outward(min(vs), -np.inf), outward(max(vs), np.inf)]
            for vs in zip(*polygon, strict=True)
        ],
    }


def clock_counts(p, rng, n):
    p = np.append(p, 1 - p.sum())
    order = np.argsort(p, kind="stable")
    p = np.maximum(p[order], 0)
    p[-1] = 1 - p[:-1].sum()
    totals = np.zeros(len(p), dtype=object)
    for _ in range(n // (2**52 * 1024)):
        totals += rng.multinomial(2**52, p, size=1024).sum(axis=0).astype(object)
    totals = totals[np.argsort(order)].tolist()
    assert sum(totals) == n
    return totals


def run():
    iv.dps = 35
    base, states, *_ = c.trials()
    cert_path = Path(c.dose.enclosure.__file__).with_name("certificate.json")
    cert = json.loads(cert_path.read_bytes())
    names = ("reference", "independent", "compensator", "held")
    gains = [base["entropy_gain_product"][0], base["entropy_gain_compensator"][0]]
    effects, scales, active = beam.flight()
    seed = 28112090
    rng = np.random.default_rng(seed)
    n = 2**72
    links = []
    for edge in beam.cuts.CHANNELS:
        score = {}
        counts = {}
        for name in names:
            y = states[name]
            scale = scales[edge]
            counts[name], score[name] = beam.port_read(
                rng,
                n,
                c.station.pair(active[edge], y) / scale,
                c.station.pair(effects[edge], y) / scale,
                scale,
                np.log(4 * len(names) * len(beam.cuts.CHANNELS) / 0.01),
            )
        record = {k: bounds(score[k] - score["reference"]) for k in names[1:]}
        links.append(
            {
                "edge": list(edge),
                "plus_minus_miss": counts,
                **record,
                "alone": infer([record], gains),
            }
        )
    decoded = infer(links, gains)
    # Decode the source before acquiring any volume or area.
    _, bins = c.lab.six.survey.return_bins()
    ids = [i for i, p in enumerate(c.lab.six.survey.ROUTES) if max(p) < 4]
    n_clock = 2**70
    time = {
        name: [
            clock_counts(p, rng, n_clock)
            for p in c.station.pair(bins[ids], states[name])
        ]
        for name in ("reference", "held")
    }
    means, _ = reg.timing_slopes(time, n_clock, (("held", "reference", 1),))
    ref, y = map(lambda k: np.array(means[k], dtype=object), ("reference", "held"))
    direction = y - ref
    path = np.array(
        [r + reg.band(0, 1) * d for r, d in zip(ref, direction, strict=True)]
    )
    area, volume = geometry_jet(path)
    delta_volume = volume.g @ direction
    # The archived global Hessian proof gives a ceiling on S+eta*A-lambda*V
    # throughout the convex preparation domain. The reference's certified
    # gap to its exact maximizer is retained when comparing finite volumes.
    drift = cert["multiplier_drift_upper"]
    lam = reg.band(51.0886428 - drift - 1e-12, 51.0886428 + drift + 1e-12)
    deduction = (
        decoded["entropy_lower"]
        - lam * delta_volume
        - base["reference_lagrangian_gap_upper"]
    ) / iv.mpf(cert["gain_interval"])
    return {
        "source_family": list(names[:3]),
        "entropy_gain_stock_lowers": gains,
        "links": links,
        "seed": seed,
        "starts_per_transfer": n,
        "blind_source_bounds": decoded,
        "starts_per_clock": n_clock,
        "time_histograms": time,
        "volume_change": bounds(delta_volume),
        "ordinary_volume_multiplier": bounds(lam),
        "minimum_area_contraction": bounds(deduction)[0],
        "independent_area_change": bounds(area.g @ direction),
        "input_sha256": {
            relpath(f, c.lab.HERE): sha256(f.read_bytes()).hexdigest()
            for f in (
                c.lab.HERE / "LINK_SOURCE_CONTRACTION.md",
                c.lab.HERE / "finite_pairing_dose.json",
                cert_path,
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    data = run()
    path = c.lab.HERE / "link_source_contraction.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    print(
        {
            k: data[k]
            for k in (
                "blind_source_bounds",
                "volume_change",
                "minimum_area_contraction",
                "independent_area_change",
            )
        }
    )
