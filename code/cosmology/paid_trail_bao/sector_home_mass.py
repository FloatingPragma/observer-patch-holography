"""Test the ordinary rest inventory inferred from overlap records on imaged galaxies."""

import argparse
import json
from fractions import Fraction
from hashlib import sha256

import numpy as np

from . import sector_mass

inventory = sector_mass.inventory
HERE = inventory.HERE


def free_placement(inputs, galaxy_masses, ratio, grav):
    # Allow every ring to use the entire budget independently. A failure
    # of this necessary bound excludes even a galaxy-specific mass profile.
    rows = {}
    for name, data in inputs.items():
        radii = np.asarray(data["radius_kpc"])
        baryons = np.asarray(data["baryons_kms"])
        lower = np.maximum(
            0, np.asarray(data["observed_kms"]) - 2 * np.asarray(data["error_kms"])
        )
        required = np.maximum(0, (lower**2 - baryons**2) * radii / grav)
        budget = ratio * galaxy_masses[name] * 1e9
        rows[name] = {
            "ordinary_budget_msun": galaxy_masses[name] * 1e9,
            "extra_budget_msun": budget,
            "minimum_extra_by_ring_msun": required.tolist(),
            "exceeds_whole_galaxy_budget": int((required > budget).sum()),
            "largest_required_per_budget": float(required.max() / budget),
        }
    return {
        "extra_per_ordinary_upper": ratio,
        "by_galaxy": rows,
        "excluded_speeds": sum(v["exceeds_whole_galaxy_budget"] for v in rows.values()),
        "excluded_curves": sum(
            v["exceeds_whole_galaxy_budget"] > 0 for v in rows.values()
        ),
    }


def run():
    path = inventory.LAB / "sector_home.json"
    chart = HERE / "record_profile_maps.json"
    standard = inventory.LAB / "thermal_standard.json"
    lab = json.loads(path.read_text())
    peaks = json.loads(chart.read_text())["peak_msun_pc2"]
    calibration = json.loads(standard.read_text())
    single, double = (calibration[k]["energy_band"] for k in ("single", "double"))
    # Price a rest-inventory unit with the measured single; divide by the
    # independently measured double that the four-port image inventory uses.
    units_per_double = [double[0] / single[1], double[1] / single[0]]
    bins = {
        k: lab["all_rates"][k] for k in ("four_home_single_units", "home_single_units")
    }
    inputs, _, _, galaxy_masses = inventory.rp.galaxy.transfer.budget.load_inputs()
    G = inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    bbn = 3 * 0.1**2 / (8 * np.pi * G) * 0.02235
    profiles = sector_mass.radial_channels(
        bins, 1.0, units_per_double, peaks, inputs, bbn
    )
    caps = {}
    for name, panel in bins.items():
        # Even if the report probability and baryon filling are unrelated,
        # p*f_response(p) has this upper bound over all observed return rates.
        upper = (
            max(max(0, r["per_rate"][1]) * r["rate"][1] for r in panel)
            / units_per_double[0]
        )
        caps[name] = {
            "ordinary_doubles_per_prepared_return_upper": upper,
            "maximum_ratio_at_any_filling": upper * (2 / 9) * (448 / 2110),
        }
    stock = lab["ordinary_rest_budget"]
    reference = float(Fraction(stock["reference_mean_single_units"]))
    dynamic = max(
        r["rate"][1] * max(0, r["per_rate"][1]) for r in bins["four_home_single_units"]
    )
    bounds = {
        "any_geometric_reporting_mixture": dynamic / reference,
        "any_population_on_returned_support": float(
            Fraction(stock["gain_per_reference_upper"])
        ),
        "any_supported_preparation_full_board_energy": float(
            Fraction(stock["board_energy_gain_per_ordinary_rest_upper"])
        ),
    }
    return {
        "sha256": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (path, chart, standard, sector_mass.pooled.sky.MAPS)
        },
        "sites_per_register": 4,
        "native_enrollment": [448, 2110],
        "imaged_galaxies": len(inputs),
        "speeds": sum(len(v["radius_kpc"]) for v in inputs.values()),
        "singles_per_double": units_per_double,
        "all_rate_bounds": caps,
        "arbitrary_spherical_placement": {
            name: free_placement(inputs, galaxy_masses, bound, G)
            for name, bound in bounds.items()
        },
        "profiles": profiles,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    result = run()
    dst = HERE / "sector_home_mass.json"
    if opt.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        inventory.rp.galaxy.check(json.loads(dst.read_text()), result)
    for k, v in result["profiles"].items():
        print(
            k,
            [(x["possible"], x["complete"]) for x in v["scores"].values()],
            result["all_rate_bounds"][k],
            v["bbn_ratio_bounds"],
        )
