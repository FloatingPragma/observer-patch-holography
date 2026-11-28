"""Construct and bound a physically adjustable source with equal impulse equivalents."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import charge_selected_flight as selected

d, lab, reg, iv = selected.d, selected.lab, selected.reg, selected.iv
mat = selected.material
HERE = selected.HERE
SEED = 28112117


def apply_phase(row, phase):
    return sum(x * iv.mpf(z) for x, z in zip(phase, row, strict=True))


def calibrated(w, control, energy):
    quad = w["meter_quadratures"]
    ordinary_impulses = [iv.mpf(control[f"shift/{sign}"]) for sign in (-1, 1)]
    assert all(0 not in x for x in ordinary_impulses)
    determinants = [
        iv.mpf(quad["shift/-1"][j]) * ordinary_impulses[1]
        - iv.mpf(quad["shift/1"][j]) * ordinary_impulses[0]
        for j in range(2)
    ]
    norm = iv.sqrt(determinants[0] ** 2 + determinants[1] ** 2)
    assert norm.a > 0
    phase = [determinants[1] / norm, -determinants[0] / norm]
    amp = apply_phase(quad["shift/-1"], phase) / ordinary_impulses[0]
    assert 0 not in amp
    sign = 1 if amp.a > 0 else -1
    phase = [sign * v for v in phase]
    target = [
        apply_phase(quad[f"shift/{q}"], phase) / op
        for q, op in zip((-1, 1), ordinary_impulses, strict=True)
    ]
    ratio = iv.mpf([max(z.a for z in target), min(z.b for z in target)])
    assert ratio.a > 0
    centers = np.array([float(v.mid) for v in determinants])
    fixed = sign * np.array([centers[1], -centers[0]]) / np.linalg.norm(centers)
    fixed *= 1 - 1e-14
    assert all(float(z.a) <= v <= float(z.b) for z, v in zip(phase, fixed, strict=True))
    return phase, ratio, fixed, lab.six.bounds(energy * ratio)


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)

    def source(name):
        return json.loads((HERE / (name + ".json")).read_bytes())

    prior = source("charge_selected_flight")
    ords = source("ordinary_service_distribution")["windows"]
    costs = source("density_material")["windows"]
    fringes = source("source_interference")["windows"]
    records = {
        k: calibrated(w, ords[k]["held_ordinary"], iv.mpf(costs[k]["ordinary_energy"]))
        for k, w in prior["windows"].items()
    }
    _, prep, *_ = mat.c.trials()
    ports, bins = lab.six.survey.return_bins()
    mixed, _, active = d.history(
        prep["reference"],
        ports,
        rng,
        {"selected_" + k: complex(*r[2]) for k, r in records.items()},
    )
    effects, scales, hits = mat.beam.flight(("shift",), include_edges=False)
    routes = lab.six.survey.ROUTES
    clocks = [f"clock/{a}/{b}" for a, b in routes]
    anchors = [j for j, r in enumerate(routes) if max(r) < 4]
    clock = 3 * np.eye(4) + np.einsum("rtbij,t->rbij", bins, reg.SCORE[:-1] - 3)
    ref = reg.src.observed(prior["source_reference_heads"], lab.SHOTS, False)
    out = {}
    for k, (enclosure, ratio, fixed, mass) in records.items():
        tau = int(k)
        name = "selected_" + k
        w = prior["windows"][k]
        q = apply_phase(
            fringes[k]["measured_price_quadratures"]["connected"], enclosure
        )
        src = {
            label: reg.src.observed(tape, lab.SHOTS, False)
            for label, tape in w["phase_standard_heads"].items()
        }
        base = tuple(
            (1 - d.DOSE) * r + d.DOSE * (u + v) / 2
            for r, u, v in zip(ref, src["plus"], src["minus"], strict=True)
        )
        predicted = {
            j: apply_phase(row, enclosure) for j, row in w["meter_quadratures"].items()
        }
        predicted_fixed = {
            j: apply_phase(row, fixed) for j, row in w["meter_quadratures"].items()
        }
        start = [selected.fringe.read(E, base, False) for E in clock]
        geom = selected.geometry_change(start, [predicted[j] for j in clocks], anchors)
        diff = {
            j: lab.six.bounds(z - ratio * iv.mpf(ords[k]["held_ordinary"][j]))
            for j, z in predicted.items()
        }
        heads = {
            j: lab.pairs(active[tau][j].sum(0), rng, False) for j in (name, "dense")
        }
        obs = {j: reg.src.observed(v, lab.SHOTS, False) for j, v in heads.items()}
        measured_q = d.DOSE * d.price_contrast(
            ref, obs[name], obs["dense"], "connected"
        )
        radar = {
            j: [
                lab.counts(np.r_[p, 1 - p.sum()], rng)
                for p in mat.c.station.pair(bins, mixed[tau][j])
            ]
            for j in (name, "dense")
        }
        means, _ = reg.timing_slopes(radar, lab.SHOTS, ((name, "dense", 1),))
        observed_geom = selected.geometry_change(
            means["dense"],
            [x - y for x, y in zip(means[name], means["dense"], strict=True)],
            anchors,
        )
        impulses = {}
        for sign in (-1, 1):
            key = (sign, "shift")
            counts = {}
            read = {}
            for j in (name, "dense"):
                counts[j], read[j] = mat.acquire(
                    rng,
                    2**76,
                    active[tau][j].sum(0),
                    effects[key],
                    hits[key],
                    scales[key],
                )
            label = f"shift/{sign}"
            scale = iv.mpf(costs[k]["ordinary_energy"]) / iv.mpf(
                ords[k]["held_ordinary"][label]
            )
            impulses[str(sign)] = {
                "counts": counts,
                "equivalent_from_meter": lab.six.bounds(
                    scale * d.DOSE * (read[name] - read["dense"])
                ),
                "equivalent_from_source_at_programmed": lab.six.bounds(
                    scale * predicted_fixed[label]
                ),
            }
        out[k] = {
            "exact_phase_enclosure": list(map(lab.six.bounds, enclosure)),
            "programmed_phase": fixed.tolist(),
            "common_positive_energy_at_exact_phase": mass,
            "ordinary_ratio_at_exact_phase": lab.six.bounds(ratio),
            "connected_charge_at_exact_phase": lab.six.bounds(q),
            "source_geometry_at_exact_phase": geom,
            "meter_minus_equal_ordinary_at_exact_phase": diff,
            "heads": heads,
            "observed_programmed_charge": lab.six.bounds(measured_q),
            "clock_counts": radar,
            "observed_programmed_geometry": observed_geom,
            "material_at_programmed_phase": impulses,
        }
    return {
        "seed": SEED,
        "material_starts": 2**76,
        "windows": out,
        "inputs": {
            n: sha256((HERE / n).read_bytes()).hexdigest()
            for n in (
                "COMMON_IMPULSE_POINTER.md",
                "charge_selected_flight.json",
                "ordinary_service_distribution.json",
                "density_material.json",
                "source_interference.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "common_impulse_pointer.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), r)
    for k, w in r["windows"].items():
        print(
            k,
            w["common_positive_energy_at_exact_phase"],
            w["connected_charge_at_exact_phase"],
            w["source_geometry_at_exact_phase"]["area"],
        )
