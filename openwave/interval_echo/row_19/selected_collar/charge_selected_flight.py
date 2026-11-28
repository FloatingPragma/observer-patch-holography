"""Forecast area and material motion from a source selected only by its recorded charge."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import density_material as material
from . import source_interference as fringe

d, lab, reg, iv = material.density, material.lab, material.reg, material.iv
HERE = material.HERE
SEED = 28112116
POWER = 72
MENU = [
    (5, 0),
    (0, 5),
    (-5, 0),
    (0, -5),
    *((a * i, b * j) for a, b in ((3, 4), (4, 3)) for i in (-1, 1) for j in (-1, 1)),
]


def choice(quadratures):
    r, i = map(iv.mpf, quadratures)
    candidates = [(a * r + b * i) / 5 for a, b in MENU]
    k = max(range(len(MENU)), key=lambda j: float(candidates[j].a))
    return (
        MENU[k],
        lab.six.bounds(candidates[k]),
        min(
            float((candidates[k].a - c.b).a) for j, c in enumerate(candidates) if j != k
        ),
    )


def geometry_change(start, delta, anchors):
    x = np.array([start[i] for i in anchors], dtype=object)
    dx = np.array([delta[i] for i in anchors], dtype=object)
    jet = fringe.pos.geometry_jet(x + np.array([reg.band(0, 1) * v for v in dx]))
    return {
        name: lab.six.bounds(j.g @ dx)
        for name, j in zip(("area", "volume"), jet, strict=True)
    }


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)

    def read_file(name):
        return json.loads((HERE / (name + ".json")).read_bytes())

    source = read_file("source_interference")
    ordinary = read_file("ordinary_service_distribution")
    cost = read_file("density_material")
    selections = {
        t: choice(source["windows"][str(t)]["measured_price_quadratures"]["connected"])
        for t in d.WINDOWS
    }
    pointers = {**fringe.PHASES}
    for t, ((a, b), _, _) in selections.items():
        pointers[f"held_{t}"] = complex(a / 5, b / 5)
    _, inputs, *_ = material.c.trials()
    x = inputs["reference"]
    ports, bins = lab.six.survey.return_bins()
    mixed, _, local = d.history(x, ports, rng, pointers)
    effects, scales, hits = material.beam.flight(material.METERS, include_edges=False)
    prices = {f"{k}/{q}": np.asarray(v) for (q, k), v in effects.items()}
    clock = 3 * np.eye(4) + np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    routes = lab.six.survey.ROUTES
    clocks = [f"clock/{a}/{b}" for a, b in routes]
    anchors = [j for j, route in enumerate(routes) if max(route) < 4]
    prices.update(dict(zip(clocks, clock, strict=True)))
    ref_tape = lab.pairs(x, rng, False)
    ref = reg.src.observed(ref_tape, lab.SHOTS, False)
    rows = {}
    for tau in d.WINDOWS:
        name = str(tau)
        held_name = f"held_{tau}"
        (a, b), charge, pick_gap = selections[tau]
        active = {k: v.sum(0) for k, v in local[tau].items()}
        standards = {k: lab.pairs(active[k], rng, False) for k in fringe.PHASES}
        obs = {k: reg.src.observed(v, lab.SHOTS, False) for k, v in standards.items()}
        quadratures = {}
        expected = {}
        for k, E in prices.items():
            # Differences in normalized preparations have trace zero.
            E = E - E.trace(axis1=-2, axis2=-1).real.sum() / 8 * np.eye(4)
            r = d.DOSE * material.forecast(E, obs["plus"], obs["minus"]) / 2
            j = d.DOSE * material.forecast(E, obs["left"], obs["right"]) / 2
            quadratures[k] = [lab.six.bounds(r), lab.six.bounds(j)]
            expected[k] = (a * r + b * j) / 5
        base = tuple(
            (1 - d.DOSE) * u + d.DOSE * (v + w) / 2
            for u, v, w in zip(ref, obs["plus"], obs["minus"], strict=True)
        )
        start = [fringe.read(E, base, False) for E in clock]
        forecast = geometry_change(start, [expected[k] for k in clocks], anchors)
        # Held acquisitions begin after choosing gamma and computing every forecast.
        heads = {k: lab.pairs(active[k], rng, False) for k in ("dense", held_name)}
        high, low = (
            reg.src.observed(heads[k], lab.SHOTS, False) for k in (held_name, "dense")
        )
        observed_q = d.DOSE * d.price_contrast(ref, high, low, "connected")
        radar = {
            k: [
                lab.counts(np.r_[v, 1 - v.sum()], rng)
                for v in material.c.station.pair(bins, mixed[tau][k])
            ]
            for k in ("dense", held_name)
        }
        means, _ = reg.timing_slopes(radar, lab.SHOTS, ((held_name, "dense", 1),))
        direct = geometry_change(
            means["dense"],
            [v - w for v, w in zip(means[held_name], means["dense"], strict=True)],
            anchors,
        )
        M = {}
        for key, E in effects.items():
            data = {}
            scores = {}
            for k in ("dense", held_name):
                data[k], scores[k] = material.acquire(
                    rng, 2**POWER, active[k], E, hits[key], scales[key]
                )
            label = f"{key[1]}/{key[0]}"
            denom = iv.mpf(ordinary["windows"][name]["held_ordinary"][label])
            value = d.DOSE * (scores[held_name] - scores["dense"])
            energy = iv.mpf(cost["windows"][name]["ordinary_energy"])
            assert 0 not in denom
            M[label] = {
                "counts": data,
                "direct": lab.six.bounds(value),
                "ordinary": lab.six.bounds(denom),
                "equivalent_from_source": lab.six.bounds(
                    energy * expected[label] / denom
                ),
                "equivalent_from_material": lab.six.bounds(energy * value / denom),
            }
        rows[name] = {
            "pointer_fifths": [a, b],
            "selection_charge": charge,
            "lower_minus_next_upper": pick_gap,
            "held_charge": lab.six.bounds(observed_q),
            "phase_standard_heads": standards,
            "held_heads": heads,
            "meter_quadratures": quadratures,
            "forecast": {k: lab.six.bounds(v) for k, v in expected.items()},
            "geometry_from_source": forecast,
            "geometry_from_radar": direct,
            "clock_counts": radar,
            "material": M,
        }
    return {
        "seed": SEED,
        "material_starts": 2**POWER,
        "head_and_radar_starts": lab.SHOTS,
        "source_reference_heads": ref_tape,
        "windows": rows,
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "CHARGE_SELECTED_FLIGHT.md",
                "source_interference.json",
                "ordinary_service_distribution.json",
                "density_material.json",
                "finite_pairing_dose.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "charge_selected_flight.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dst.read_bytes()), result)
    for t, r in result["windows"].items():
        print(
            t,
            r["pointer_fifths"],
            r["held_charge"],
            r["geometry_from_radar"]["area"],
            [r["material"][f"shift/{q}"]["equivalent_from_material"] for q in (-1, 1)],
        )
