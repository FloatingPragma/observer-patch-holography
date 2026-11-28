"""Predict held geometry and material response with an ordinary spatial energy calibration."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import ordinary_service_distribution as service
from . import regional_density as regions

HERE, lab, reg, iv = regions.HERE, regions.lab, regions.reg, regions.iv
SEED = 28112120


def read(name):
    return json.loads((HERE / (name + ".json")).read_bytes())


def acquire():
    rng = np.random.default_rng(SEED)
    d = regions.density
    _, inputs, *_ = regions.c.trials()
    ports = list(service.lab.six.survey.return_bins()[0])
    matched = read("six_meter_pointer")
    settings = {t: complex(*v) for t, v in matched["central_settings"].items()}
    _, _, held = d.history(
        inputs["reference"], ports, rng, {"p" + t: g for t, g in settings.items()}
    )

    def tape(x):
        return regions.joint_tape(x, rng)

    targets = {
        name: tape(
            sum(
                held[int(t)][("p" + t) if name == "overlap" else "dense"]
                for t in settings
            )
            / len(settings)
        )
        for name in ("overlap", "dense")
    }
    standards = {}
    for home in ports:
        for sign in (-1, 1):
            runs = {
                arm: d.history(
                    inputs[arm], (home,), rng, signs=(sign,), read_ports=ports
                )[2]
                for arm in ("reference", "ordinary")
            }
            for t in settings:
                standards[f"{t}/{home}/{sign}"] = {
                    arm: tape(z[int(t)]["dense"]) for arm, z in runs.items()
                }
    return {
        "seed": SEED,
        "starts": lab.SHOTS,
        "ports": ports,
        "reference_heads": lab.pairs(inputs["reference"], rng, False),
        "standard_counts": standards,
        "target_counts": targets,
    }


def decode(data):
    iv.dps = 35
    matched, clocks, service_record = (
        read(n)
        for n in (
            "six_meter_pointer",
            "six_meter_clocks",
            "ordinary_service_distribution",
        )
    )
    receivers = [str(p) for p in data["ports"]] + ["complement"]
    prices = regions.costs(
        reg.src.observed(data["reference_heads"], data["starts"], False)
    )

    def charges(high, low):
        return {
            n: regions.local_contrasts(
                prices, regions.at_port(high, i), regions.at_port(low, i)
            )
            for i, n in enumerate(receivers)
        }

    held = charges(data["target_counts"]["overlap"], data["target_counts"]["dense"])
    controls = {
        k: charges(v["ordinary"], v["reference"])
        for k, v in data["standard_counts"].items()
    }

    def full(z):
        return z["connected"] + z["marginal"]

    def energy(z):
        return sum(v["energy"] for v in z.values())

    energy_per_history = {k: energy(v) for k, v in controls.items()}
    assert all(v.a > 0 for v in energy_per_history.values())
    ordinary_conversion = sum(
        full(v) for z in controls.values() for v in z.values()
    ) / sum(energy_per_history.values())
    # A conversion whose interval includes zero cannot identify an energy.
    if ordinary_conversion.a <= 0 <= ordinary_conversion.b:
        return {
            "conversion": lab.six.bounds(ordinary_conversion),
            "status": "ordinary slope unresolved",
        }
    means, _ = reg.timing_slopes(
        clocks["clock_counts"], clocks["starts"], (("held", "dense", 1),)
    )
    dose = regions.density.DOSE
    base = [
        (1 - dose) * u + dose * v
        for u, v in zip(means["reference"], means["dense"], strict=True)
    ]
    directions = [
        dose * (u - v) for u, v in zip(means["held"], means["dense"], strict=True)
    ]
    ratios = regions.c.clock_map()[1:] - regions.c.clock_map()[0]
    names = [f"ratio/{i}/{j}" for i, j in regions.c.EDGES[1:]]

    def shape(timing):
        return {
            name: lab.six.bounds(
                sum(
                    float(c) * val / m
                    for c, val, m in zip(row, timing, base, strict=True)
                    if c
                )
            )
            for name, row in zip(names, ratios, strict=True)
        }

    material = matched["labels"]
    cols = {}
    for k, local in controls.items():
        t, sender = k.split("/", 1)
        ordinary = service_record["windows"][t]["ordinary_columns"][sender]
        route_deltas = [
            iv.mpf(ordinary[f"clock/{i}/{j}"]) for i, j in lab.six.survey.ROUTES
        ]
        cols[k] = {
            **{n: ordinary[n] for n in material},
            **shape(route_deltas),
            **{"site/" + n: lab.six.bounds(v["energy"]) for n, v in local.items()},
        }
    target_shape = shape(directions)
    targets = {
        n: matched["material_at_center"][n]["forecast_at_center"] for n in material
    }
    direct = {n: matched["material_at_center"][n]["observed"] for n in material}
    sites = ["site/" + n for n in receivers]
    groups = {
        "sites": sites,
        "shape": names,
        "material": material,
        "shape_material": names + material,
        "sites_shape": sites + names,
        "sites_material": sites + material,
        "all": sites + names + material,
    }
    sources = {}
    for name, price in (("connected", lambda v: v["connected"]), ("full", full)):
        inferred = {n: price(v) / ordinary_conversion for n, v in held.items()}
        total = sum(inferred.values())
        readings = {
            **target_shape,
            **targets,
            **{"site/" + n: lab.six.bounds(v) for n, v in inferred.items()},
        }
        direct_readings = {**readings, **direct}
        envelope = {}
        for meter in [*names, *material]:
            ratios_by_history = [
                iv.mpf(z[meter]) / energy_per_history[k] for k, z in cols.items()
            ]
            if total.b < 0:
                # All ordinary histories carry strictly positive total energy.
                continue
            unit = iv.mpf(
                [
                    min(v.a for v in ratios_by_history),
                    max(v.b for v in ratios_by_history),
                ]
            )
            prediction = iv.mpf([max(total.a, 0), total.b]) * unit

            def excluded(a):
                return bool(prediction.a > a[1] or prediction.b < a[0])

            envelope[meter] = {
                "prediction": lab.six.bounds(prediction),
                "excludes_programmed_read": excluded(direct_readings[meter]),
                "excludes_source_forecast": excluded(readings[meter])
                if meter in material
                else None,
            }
        sources[name] = {
            "inferred_total_energy": lab.six.bounds(total),
            "target": readings,
            "direct_material": direct,
            "source_only_prediction": envelope,
            "comparisons": {
                g: service.separate(cols, readings, ns) for g, ns in groups.items()
            },
            "direct_comparisons": {
                g: service.separate(cols, direct_readings, groups[g])
                for g in ("material", "all")
            },
        }
    return {
        "ordinary_full_modular_per_energy": lab.six.bounds(ordinary_conversion),
        "ordinary_total_energy": {
            k: lab.six.bounds(v) for k, v in energy_per_history.items()
        },
        "held_sources": {
            n: {name: lab.six.bounds(v) for name, v in z.items()}
            for n, z in held.items()
        },
        "ordinary_regional_sources": {
            k: {
                n: {name: lab.six.bounds(v) for name, v in z.items()}
                for n, z in loc.items()
            }
            for k, loc in controls.items()
        },
        "groups": groups,
        "columns": cols,
        "sources": sources,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    receipt = acquire()
    receipt["result"] = decode(receipt)
    receipt["inputs"] = {
        n: sha256((HERE / n).read_bytes()).hexdigest()
        for n in (
            "DISTRIBUTED_SOURCE_PROTOCOL.md",
            "six_meter_pointer.json",
            "six_meter_clocks.json",
            "ordinary_service_distribution.json",
        )
    }
    dst = HERE / "distributed_source.json"
    if args.write:
        dst.write_text(json.dumps(receipt, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(dst.read_bytes()), receipt)
    out = receipt["result"]
    print(
        "ordinary modular price",
        out.get("ordinary_full_modular_per_energy", out.get("conversion")),
    )
    for name, x in out.get("sources", {}).items():
        print(
            name,
            x["inferred_total_energy"],
            {g: v["status"] for g, v in x["comparisons"].items()},
        )
