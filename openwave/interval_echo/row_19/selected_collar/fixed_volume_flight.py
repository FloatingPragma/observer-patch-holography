"""Use charge-bearing energy-neutral controls to fix survey volume and predict its area."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import zero_flight_work as zero

HERE, lab, reg, iv, d, mat = zero.HERE, zero.lab, zero.reg, zero.iv, zero.d, zero.mat
fringe = zero.impulse.selected.fringe
TIMES = ("2", "4")
SEED = 28112128


def run():
    iv.dps = 35
    read = zero.energy.flow.dist.read
    neutral, prior, source, costs, ordinary = (
        read(k)
        for k in (
            "zero_flight_work",
            "charge_selected_flight",
            "source_interference",
            "density_material",
            "ordinary_service_distribution",
        )
    )
    certfile = Path(mat.c.dose.enclosure.__file__).with_name("certificate.json")
    gain = iv.mpf(json.loads(certfile.read_bytes())["gain_interval"])
    ports, bins = lab.six.survey.return_bins()
    routes = lab.six.survey.ROUTES
    anchors = [i for i, route in enumerate(routes) if max(route) < 4]
    labels = [f"clock/{a}/{b}" for a, b in routes]
    effects = 3 * np.eye(4) + np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    reference = reg.src.observed(source["reference_head_counts"], lab.SHOTS, False)
    dense_offsets, lengths, baselines, meters = {}, {}, {}, {}
    for t in TIMES:
        phase = list(map(iv.mpf, neutral["windows"][t]["phase"]))
        std = {
            k: reg.src.observed(v, lab.SHOTS, False)
            for k, v in source["windows"][t]["head_quadrature_counts"].items()
        }
        dense = tuple(
            (u + v) / 2 for u, v in zip(std["plus"], std["minus"], strict=True)
        )
        dense_offsets[t] = fringe.trace_bound(fringe.difference(dense, reference))
        lengths[t] = sum(
            abs(g).b * fringe.trace_bound(fringe.difference(std[a], std[b])) / 2
            for g, (a, b) in zip(
                phase, (("plus", "minus"), ("left", "right")), strict=True
            )
        )
        baselines[t] = [fringe.read(e, dense, False) for e in effects]
        meters[t] = {
            k: zero.impulse.apply_phase(v, phase)
            for k, v in prior["windows"][t]["meter_quadratures"].items()
        }
    base = [(a + b) / 2 for a, b in zip(*baselines.values(), strict=True)]

    def mix(a, b, u):
        return (u * b - a) / 2

    def geometry(u):
        delta = [mix(meters["2"][k], meters["4"][k], u) for k in labels]
        return zero.impulse.selected.geometry_change(base, delta, anchors)

    # Refine both certified sides. The guard width grows if the middle
    # point's sign is unresolved, and never invalidates a retained bracket.
    lo, hi = 0.125, 0.25
    left, right = geometry(lo)["volume"], geometry(hi)["volume"]
    assert left[0] > 0 > right[1]
    for _ in range(20):
        x, y = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        vx, vy = geometry(x)["volume"], geometry(y)["volume"]
        old = lo, hi
        for z, v in ((x, vx), (y, vy)):
            if v[0] > 0 and z > lo:
                lo, left = z, v
            if v[1] < 0 and z < hi:
                hi, right = z, v
        if old == (lo, hi):
            break
    u = iv.mpf([lo, hi])
    full, connected = (
        mix(*(iv.mpf(neutral["windows"][t]["modular"][k]) for t in TIMES), u)
        for k in ("full", "connected")
    )
    length = (lengths["2"] + hi * lengths["4"]) / 2
    radius = sum(dense_offsets.values()) / 2 + length
    err = (
        source["reference_gradient_error_upper"]
        + source["mixed_geometry_Hessian_upper"] * radius
    ) * length + 1e-12
    # Build the reported prediction from the same outward-rounded inputs
    # an independent observer receives in the JSON record.
    error_upper = lab.six.bounds(err)[1]
    area = (
        -iv.mpf(lab.six.bounds(full)) + iv.mpf([-error_upper, error_upper])
    ) / iv.mpf(lab.six.bounds(gain))
    mid = (lo + hi) / 2
    knobs = {
        t: complex(*neutral["windows"][t]["programmed_phase"]) * a
        for t, a in zip(TIMES, (-1, mid), strict=True)
    }
    rng = np.random.default_rng(SEED)
    _, stocks, *_ = mat.c.trials()
    mixed, _, active = d.history(
        stocks["reference"],
        ports,
        rng,
        {"p" + t: g for t, g in knobs.items()},
        windows=tuple(map(int, TIMES)),
    )

    def average(table, arm, local):
        pair = [table[int(t)]["dense" if arm == "dense" else "p" + t] for t in TIMES]
        return sum(x.sum(0) if local else x for x in pair) / 2

    arms = ("held", "dense")
    heads = {k: lab.pairs(average(active, k, True), rng, False) for k in arms}
    observed = [reg.src.observed(heads[k], lab.SHOTS, False) for k in arms]
    radar = {
        k: [
            lab.counts(np.r_[p, 1 - p.sum()], rng)
            for p in mat.c.station.pair(bins, average(mixed, k, False))
        ]
        for k in arms
    }
    means, _ = reg.timing_slopes(radar, lab.SHOTS, (("held", "dense", 1),))
    direct = zero.impulse.selected.geometry_change(
        means["dense"],
        [a - b for a, b in zip(means["held"], means["dense"], strict=True)],
        anchors,
    )
    beam, scales, hits = mat.beam.flight(("shift",), include_edges=False)
    ordinary_e = sum(iv.mpf(costs["windows"][t]["ordinary_energy"]) for t in TIMES) / 2
    material = {}
    for sign in (-1, 1):
        key = sign, "shift"
        n = f"shift/{sign}"
        common = (
            sum(iv.mpf(ordinary["windows"][t]["held_ordinary"][n]) for t in TIMES) / 2
        )
        assert 0 not in common
        score = {}
        tape = {}
        for k in arms:
            tape[k], score[k] = mat.acquire(
                rng, 2**76, average(active, k, True), beam[key], hits[key], scales[key]
            )
        material[n] = {
            "counts": tape,
            "equivalent_at_exact_volume": lab.six.bounds(
                ordinary_e * mix(meters["2"][n], meters["4"][n], u) / common
            ),
            "equivalent_at_programmed_center": lab.six.bounds(
                ordinary_e * d.DOSE * (score["held"] - score["dense"]) / common
            ),
        }
    return {
        "times": list(TIMES),
        "weights": [0.5, 0.5],
        "seed": SEED,
        "starts": lab.SHOTS,
        "material_starts": 2**76,
        "amplitude_bracket": [lo, hi],
        "volume_at_ends": [left, right],
        "central_phases": {t: [g.real, g.imag] for t, g in knobs.items()},
        "connected_at_exact_volume": lab.six.bounds(connected),
        "full_at_exact_volume": lab.six.bounds(full),
        "stationarity_error_upper": error_upper,
        "ordinary_gain": lab.six.bounds(gain),
        "area_at_exact_volume": lab.six.bounds(area),
        "area_from_source_routes": geometry(u)["area"],
        "heads": heads,
        "clocks": radar,
        "direct_geometry": direct,
        "material": material,
        "direct_modular": {
            k: lab.six.bounds(d.DOSE * d.price_contrast(reference, *observed, k))
            for k in ("full", "connected")
        },
        "inputs": {
            k: sha256((HERE / k).read_bytes()).hexdigest()
            for k in (
                "FIXED_VOLUME_FLIGHT.md",
                "zero_flight_work.json",
                "charge_selected_flight.json",
                "source_interference.json",
                "density_material.json",
                "ordinary_service_distribution.json",
            )
        },
        "equilibrium_certificate_sha256": sha256(certfile.read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    path = HERE / "fixed_volume_flight.json"
    if args.write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    print(
        {
            k: result[k]
            for k in (
                "amplitude_bracket",
                "volume_at_ends",
                "connected_at_exact_volume",
                "full_at_exact_volume",
                "area_at_exact_volume",
                "direct_geometry",
            )
        }
    )
