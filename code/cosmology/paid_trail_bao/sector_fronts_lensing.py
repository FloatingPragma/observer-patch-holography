"""Test independently priced single and collective record covariances in isolated-lens apertures."""

import argparse
import json
import math
from fractions import Fraction
from hashlib import sha256

import numpy as np
from cosmology.rar_deep_regime.compensation import measure as apertures

from . import sector_passage_sky as galaxy


def run():
    here = galaxy.HERE
    lab = galaxy.inventory.LAB
    sky = here / "sector_passage_sky.json"
    source = lab / "sector_passage.json"
    aggregate = lab / "sector_fronts_collective.json"
    published = apertures.HERE / "inputs.json"
    data = json.loads(sky.read_bytes())["phase_clock_spherical_transfers"]["double"]
    s = json.loads(source.read_bytes())["spatial_modular_source"]
    reg = json.loads(aggregate.read_bytes())
    observed = apertures.run(json.loads(published.read_bytes()))
    D = s["mean_coefficients"][-1]
    coeff_error = s["block_and_calibration"]["coefficient_error_upper"] + 1e-12
    # The home-20 clock coefficient is the positive ordinary delivery standard.
    # Use its independently calibrated mean and a conservative sum of the two
    # coefficient error allowances. A double grants the larger inferred mass.
    price = (
        (D + coeff_error) / (s["ordinary_mean_response"] - 2 * coeff_error) / (37 / 28)
    )
    n = reg["size"]
    base = float(Fraction(reg["ordinary_covariance"]))
    native = reg["forecasts"]["one_eighth"]
    low_price = (D - coeff_error) / (s["ordinary_mean_response"] + 2 * coeff_error)
    high_price = price * (37 / 28)
    bins = []
    for row in reg["continuous_cadence"]["bins"]:
        product = np.multiply.outer(row["per_rate"], [low_price, high_price])
        bins.append({"rate": row["rate"], "per_rate": [product.min(), product.max()]})
    max_register_mass = max(r["rate"][1] * max(0, r["per_rate"][1]) for r in bins) / (
        37 / 28
    )
    # For totals in [0,n], Var <= n^2/4. Independent ordinary returns
    # supply exactly n*base as the control covariance. The population cap
    # applies even if every measured member in the source is correlated.
    models = {
        "additive_connected_all_rates": data["connected_geometry"]["free_placement"][
            "extra_per_ordinary_upper"
        ],
        "additive_with_single_home_terms_all_rates": data["full_geometry"][
            "free_placement"
        ]["extra_per_ordinary_upper"],
        "collective_native_any_subgroup_at_eighth": price
        * max(0, native["all_subgroup_sizes_per_member_upper"]),
        "arbitrary_dependence_at_register_size": price * (n / 4 - base),
        "collective_native_any_subgroup_all_rates": max(
            max_register_mass,
            data["connected_geometry"]["free_placement"]["extra_per_ordinary_upper"],
        ),
        "any_independent_member_reader_at_same_ordinary_covariance": price
        * reg["arbitrary_member_read_covariance_gain_upper"],
    }
    out = {}
    for label, p in observed["apertures"].items():
        obs = p["cases"]["two_systematic_double_baryons"]
        limit = np.array(obs["min_exterior_debit_over_allowed_galaxy_baryons"])
        budget = np.array(obs["subtracted_baryon_cap_msun"])
        need = np.array(obs["min_net_excess_and_exterior_debit_msun"])
        out[label] = {
            "radius_kpc": p["radius_kpc"],
            "required_extra_to_baryon_lower": limit.tolist(),
            "required_extra_mass_msun": need.tolist(),
            "smallest_possible_correlated_group_by_stack": [
                max(1, math.ceil(4 * (x / price + base))) for x in limit
            ],
            "models": {
                name: {
                    "max_extra_mass_msun": (r * budget).tolist(),
                    "mass_shortfall_msun": np.maximum(0, need - r * budget).tolist(),
                    "excluded_stacks": int(sum(r < limit)),
                    "required_factor_over_maximum": (limit / r).tolist()
                    if r > 0
                    else None,
                }
                for name, r in models.items()
            },
        }
    maps = here / "record_profile_maps.json"
    peaks = json.loads(maps.read_bytes())["peak_msun_pc2"]
    inputs, _, _, masses = galaxy.inventory.rp.galaxy.transfer.budget.load_inputs()
    G = galaxy.inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    bbn = 3 * 0.1**2 / (8 * np.pi * G) * 0.02235
    radial = galaxy.sector_mass.radial_channels(
        {"whole_register": bins}, 1.0, [1, 1], peaks, inputs, bbn
    )["whole_register"]
    radial["free_placement"] = galaxy.sector_home_mass.free_placement(
        inputs, masses, models["collective_native_any_subgroup_all_rates"], G
    )
    # A mapped baryon tile has ordinary mass p times its capacity. The
    # source uses 6p^2(1-p)^2 eligible configurations times the measured
    # mass per report and the four-home area conversion.
    sparse = bins[0]
    end = sparse["rate"][1]
    k = 1.5 * 448 / 2110
    sparse_ratio = {
        "maximum_filling": end,
        "extra_over_baryons_divided_by_p_squared": [
            k * (1 - end) ** 2 * sparse["per_rate"][0],
            k * sparse["per_rate"][1],
        ],
    }
    assert (
        sparse["rate"][0] == 0
        and sparse_ratio["extra_over_baryons_divided_by_p_squared"][0] > 0
    )
    return {
        "sha256": {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                source,
                aggregate,
                sky,
                published,
                maps,
                galaxy.sector_mass.pooled.sky.MAPS,
                here / "source_profile.json",
            )
        },
        "covariance_price_per_baryon_upper": price,
        "ordinary_register_baseline": "independent reported native returns",
        "maximum_extra_to_baryon": models,
        "stacks": out,
        "predicted_imaged_orbits": radial,
        "sparse_mapped_baryon_tiles": sparse_ratio,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    s = run()
    dest = galaxy.HERE / "sector_fronts_lensing.json"
    if opt.write:
        dest.write_text(json.dumps(s, indent=2, allow_nan=False) + "\n")
    else:
        galaxy.inventory.rp.galaxy.check(json.loads(dest.read_bytes()), s)
    print(
        s["maximum_extra_to_baryon"],
        s["stacks"]["outer"]["smallest_possible_correlated_group_by_stack"],
    )
