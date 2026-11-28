"""Use stationary native motion records to predict the excess at observed galaxy radii."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np

from . import closed_inventory as inventory
from . import nested_budget
from . import passage_images
from . import thermal_mass

HERE = inventory.HERE
LAB = inventory.LAB


def branching_profiles(eligible, inputs, local):
    # The per-accepted-patch source is fixed by the paired record and common
    # to all images. Free its positive mass conversion to test the shape alone.
    G = inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    panels = {}
    for tracer, cs in eligible.items():
        for cap, menu in cs.items():
            common = []
            allowed = {}
            for name, counts in menu["closed"].items():
                tab = inputs[name]
                r, v, b, e = [
                    np.asarray(tab[k])
                    for k in ("radius_kpc", "observed_kms", "baryons_kms", "error_kms")
                ]
                x = G * local * np.array(counts) / r
                lo = np.maximum(v - 2 * e, 0) ** 2 - b * b
                hi = (v + 2 * e) ** 2 - b * b
                rows = [
                    (max(0.0, a), b) for a, b in passage_images.intervals(x, lo, hi)
                ]
                lower, upper = max(a for a, b in rows), min(b for a, b in rows)
                allowed[name] = [lower, upper] if lower <= upper else None
                common.extend(rows)
            # Endpoint counting with closed intervals, including any zero-field rows.
            candidates = {0.0} | {
                v for a, b in common if a <= b for v in (a, b) if np.isfinite(v)
            }
            best, k = max((sum(a <= v <= b for a, b in common), -v) for v in candidates)
            panels[tracer + "/" + cap] = dict(
                max_rings_one_global_gain=int(best),
                favorable_gain_double_masses_per_bit=float(-k),
                per_galaxy_complete_gain=allowed,
                whole_curves_even_with_individual_gains=sum(
                    v is not None for v in allowed.values()
                ),
            )
    return panels


def run():
    motion = json.loads((LAB / "moving_records.json").read_bytes())
    eligible = json.loads((HERE / "closed_inventory_maps.json").read_bytes())
    disks = json.loads((HERE / "nested_budget_maps.json").read_bytes())
    # The maps contain physical mass per enrolled closed attempt, using the
    # observed gas or gas plus starlight, independent of the measured orbits.
    inputs, *_ = inventory.rp.galaxy.transfer.budget.load_inputs()
    blank = float(Q(motion["blank_double_bits"]))
    loaded = motion["changes"]["ordinary"]["sweep_bits"]
    source = motion["changes"]["source"]["sweep_bits"]
    baseline = nested_budget.evaluate(disks)
    out = {}
    for name, standard in (("bare_double", blank), ("added_double", loaded)):
        transfer = source / standard
        calibration = (
            [blank, blank]
            if name == "bare_double"
            else motion["changes"]["ordinary"]["sweep_bits_interval"]
        )
        sl, sh = motion["changes"]["source"]["sweep_bits_interval"]
        limits = [
            float(np.nextafter(sl / calibration[1], -np.inf)),
            float(np.nextafter(sh / calibration[0], np.inf)),
        ]
        fitted = thermal_mass.forecast(eligible, limits, inputs)
        forecasts = {}
        for tracer, cs in fitted.items():
            for cap, rows in cs.items():
                model = rows["closed"]
                forecasts[tracer + "/" + cap] = dict(
                    within_two=model["within_two"],
                    whole_curves=sum(
                        g["within_two"] == len(inputs[n]["radius_kpc"])
                        for n, g in model["galaxies"].items()
                    ),
                    counts={n: g["within_two"] for n, g in model["galaxies"].items()},
                    maximum_baryon_speed_change_kms=model[
                        "maximum_baryon_speed_change_kms"
                    ],
                )
        # This envelope grants *any* attainable stationary mixture and any
        # distribution of fuel across patches, with a common service rate.
        spec = motion["spectrum"]
        assert spec["0"]["maximum_upper"] == 0
        lower = min(Q(v["minimum_lower"]) / int(k) for k, v in spec.items() if int(k))
        upper = max(Q(v["maximum_upper"]) / int(k) for k, v in spec.items() if int(k))
        ratio = float(np.nextafter(float((upper - lower) / Q(calibration[0])), np.inf))
        G = inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
        galaxies = {}
        for n, im in disks.items():
            tab = inputs[n]
            r, b, v, e = [
                np.array(tab[k])
                for k in ("radius_kpc", "baryons_kms", "observed_kms", "error_kms")
            ]
            lo = r / G * (np.maximum(v - 2 * e, 0) ** 2 - b * b)
            hi = r / G * ((v + 2 * e) ** 2 - b * b)
            j = list(r).index(im["start_kpc"])
            need = baseline["galaxies"][n]["nonnegative_mass_integral_lower"]
            panels = {}
            for tracer, m in im["tracers"].items():
                bound = ratio * np.array(m["enclosed_msun"])
                allowed = (lo <= bound) & (hi >= -bound)
                allmass = ratio * m["imaged_msun"] * (r[-1] ** 2 - r[j] ** 2) / 2
                panels[tracer] = dict(
                    rings=len(r),
                    compatible_rings=int(allowed.sum()),
                    whole_curve_compatible=bool(allowed.all()),
                    required_per_max_local_integral=need
                    / (ratio * m["integrated_rM_msun_kpc2"]),
                    required_per_max_relocated_integral=need / allmass,
                )
            galaxies[n] = panels
        out[name] = dict(
            standard_bits=standard,
            extra_double_masses_per_accepted=transfer,
            forecasts=forecasts,
            stationary_mass_ratio_upper=ratio,
            galaxies=galaxies,
            base_lcdm_central_to_maximum=(0.1200 / 0.02237) / ratio,
            envelope={
                t: dict(
                    compatible_rings=sum(
                        r[t]["compatible_rings"] for r in galaxies.values()
                    ),
                    whole_curves=sum(
                        r[t]["whole_curve_compatible"] for r in galaxies.values()
                    ),
                    excluded_after_relocation=sum(
                        bool(r[t]["required_per_max_relocated_integral"] > 1)
                        for r in galaxies.values()
                    ),
                )
                for t in inventory.rp.TRACERS
            },
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                LAB / "MOVING_RECORDS.md",
                LAB / "moving_records.json",
                HERE / "closed_inventory_maps.json",
                HERE / "nested_budget_maps.json",
                HERE / "source_profile.json",
            )
        },
        calibrations=out,
        branching_profile_shapes=branching_profiles(
            eligible,
            inputs,
            motion["branch_balance"]["changes"]["source"]["whole_sweep_residual_bits"],
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opts = cli.parse_args()
    result = run()
    path = HERE / "moving_budget.json"
    if opts.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        inventory.rp.galaxy.check(json.loads(path.read_bytes()), result)
    for k, v in result["calibrations"].items():
        print(
            k,
            v["extra_double_masses_per_accepted"],
            v["stationary_mass_ratio_upper"],
            v["envelope"],
            [(f["within_two"], f["whole_curves"]) for f in v["forecasts"].values()],
        )
