"""Weigh native record logarithms with the ordinary thermometer and insertion spectra."""

import argparse
from hashlib import sha256
import json

from mpmath import iv
import numpy as np

from . import thermal_mass as mass

rp = mass.rp
HERE = rp.HERE
LAB = mass.LAB


def ratio_ceiling(classes, bins=1024):
    """Bound absolute nats per loaded owner for every filling and marginal-preserving repair."""
    cells = [[iv.mpf(x) / r["words"] for x in r["magnitude_counts"]] for r in classes]
    assert classes[0] == classes[4] and classes[1] == classes[3]
    assert all(v == 0 for v in [*cells[0][1:], cells[1][3]])
    max_nats = iv.mpf(0)
    for i in range(bins):
        p = iv.mpf([i, i + 1]) / bins
        t = p * (1 - p)
        # Enclose p(1-p) also at its stationary point; enforce its exact range.
        t = iv.mpf([t.a, min(t.b, iv.mpf(1) / 4)])
        s = 1 - 2 * t
        aa = (
            cells[0][0] * (s**2 - 2 * t**2)
            + 4 * cells[1][0] * t * s
            + 6 * cells[2][0] * t**2
        )
        bt, ct = [4 * cells[1][j] * s + 6 * cells[2][j] * t for j in (1, 2)]
        dt2 = 6 * cells[2][3]
        # Cancel both powers of t analytically, including the limiting endpoints.
        ell = iv.log(aa * dt2 / (bt * ct))
        marginal = min((bt + dt2 * t).b, (ct + dt2 * t).b)
        max_nats = max(max_nats, (abs(ell) * (1 - p) * marginal / 4).b)
    return float(np.nextafter(float(max_nats), np.inf))


def run():
    iv.dps = 45
    prices = json.loads((LAB / "equilibrium_odds.json").read_bytes())
    ordinary = json.loads((HERE / "occupied_mass.json").read_bytes())["results"]
    maps = json.loads((HERE / "record_profile_maps.json").read_bytes())
    bank = json.loads((LAB / "record_density.json").read_bytes())
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    beta = iv.mpf(prices["beta"])
    beta0 = iv.log(iv.mpf(8) / 7)
    assert beta.a <= beta0.a <= beta0.b <= beta.b
    upper_nats = ratio_ceiling(bank["classes"])
    benchmark = json.loads((HERE / "thermal_mass.json").read_bytes())[
        "planck_cold_per_baryon"
    ]
    fields = {}
    for tracer, family in maps["enclosed"].items():
        fields[tracer] = {}
        for capacity, menus in family.items():
            # Four loaded owners at full occupancy represent the physical cell mass.
            scale = int(capacity) * maps["peak_msun_pc2"][tracer] / 4
            fields[tracer][capacity] = {
                menu: {
                    name: (np.asarray(arr) * scale).tolist()
                    for name, arr in galaxies.items()
                }
                for menu, galaxies in menus.items()
            }
    tests = {}
    for tag, row in ordinary.items():
        energy = iv.mpf(row["ordinary"]["spectroscopic_interval"])
        price = 1 / (beta * energy)
        conversion = [
            float(np.nextafter(float(price.a), -np.inf)),
            float(np.nextafter(float(price.b), np.inf)),
        ]
        limit = float(
            np.nextafter(float((iv.mpf(upper_nats) / (beta * energy)).b), np.inf)
        )
        scored = mass.forecast(fields, conversion, inputs)
        fits = {}
        for tracer, family in scored.items():
            for cap, menus in family.items():
                for menu, panel in menus.items():
                    counts = {n: d["within_two"] for n, d in panel["galaxies"].items()}
                    fits[f"{tracer}/{cap}/{menu}"] = dict(
                        within_two=panel["within_two"],
                        complete=sum(
                            v == len(inputs[n]["radius_kpc"]) for n, v in counts.items()
                        ),
                        maximum_speed_change_kms=panel[
                            "maximum_baryon_speed_change_kms"
                        ],
                        median_abs_kms=panel["median_abs_kms"],
                        by_galaxy=counts,
                    )
        tests[tag] = dict(
            ordinary_energy=row["ordinary"]["spectroscopic_interval"],
            mass_per_nat_in_ordinary_units=conversion,
            mass_ratio_ceiling=limit,
            forecasts=fits,
        )
    limit = max(v["mass_ratio_ceiling"] for v in tests.values())
    lower_factor = float(np.nextafter(float((iv.mpf(12000) / 2237 / limit).a), -np.inf))
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                LAB / "EQUILIBRIUM_ODDS.md",
                LAB / "record_density.json",
                LAB / "equilibrium_odds.json",
                HERE / "record_profile_maps.json",
                HERE / "occupied_mass.json",
                HERE / "thermal_mass.json",
            )
        },
        ordinary_owners_per_full_cell=4,
        filling_intervals=1024,
        nats_per_loaded_owner_ceiling=upper_nats,
        galaxies=len(inputs),
        rings=sum(len(v["radius_kpc"]) for v in inputs.values()),
        cosmic_benchmark=benchmark,
        minimum_cosmic_shortfall=lower_factor,
        baryons_within_two=sum(
            int(
                np.sum(
                    abs(np.asarray(v["baryons_kms"]) - v["observed_kms"])
                    <= 2 * np.asarray(v["error_kms"])
                )
            )
            for v in inputs.values()
        ),
        results=tests,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    result = run()
    path = HERE / "charge_temperature.json"
    if opt.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        rp.galaxy.check(json.loads(path.read_bytes()), result)
    for name, v in result["results"].items():
        rows = v["forecasts"].values()
        print(
            name,
            v["mass_ratio_ceiling"],
            sorted({(r["within_two"], r["complete"]) for r in rows}),
            max(r["maximum_speed_change_kms"] for r in rows),
        )
