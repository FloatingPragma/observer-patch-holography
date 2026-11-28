"""Transfer the calibrated regional modular reading to galaxy images and lensing apertures."""

import argparse
import json
from hashlib import sha256

import numpy as np
from cosmology.rar_deep_regime.compensation import measure as apertures

from . import sector_home_mass, sector_mass

inventory = sector_mass.inventory
HERE = sector_mass.HERE


def run():
    path = inventory.LAB / "sector_register_source.json"
    standard = inventory.LAB / "thermal_standard.json"
    maps = HERE / "record_profile_maps.json"
    lens = apertures.HERE / "inputs.json"
    lab = json.loads(path.read_text())
    price = lab["any_common_observation_of_product_references"]
    paired = lab["fixed_volume_pooled_home_reading"]
    double = json.loads(standard.read_text())["double"]["energy_band"]
    inputs, _, _, masses = inventory.rp.galaxy.transfer.budget.load_inputs()
    grav = inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    bbn = 3 * 0.1**2 / (8 * np.pi * grav) * 0.02235
    bins = [
        {
            "rate": r["rate"],
            "per_rate": (
                np.array(r["per_rate"]) / paired["members_in_reading"]
            ).tolist(),
        }
        for r in paired["bins"]
    ]
    profiles = sector_mass.radial_channels(
        {"distinct_pair": bins},
        price["beta"],
        double,
        json.loads(maps.read_text())["peak_msun_pc2"],
        inputs,
        bbn,
    )["distinct_pair"]
    budgets = {
        "any_pair_rate_and_radial_placement": max(
            max(abs(np.array(r["per_rate"]))) * r["rate"][1] for r in bins
        )
        / (price["beta"] * price["minimum_ordinary_inventory_energy_per_member"]),
        "any_register_any_common_readout": price["absolute_extra_to_ordinary_upper"],
        "difference_of_two_connected_contrasts": 2
        * price["absolute_extra_to_ordinary_upper"],
    }
    observed = apertures.run(json.loads(lens.read_text()))
    stacks = {}
    for key, r in observed["apertures"].items():
        need = np.array(
            r["cases"]["two_systematic_double_baryons"][
                "min_exterior_debit_over_allowed_galaxy_baryons"
            ]
        )
        stacks[key] = {
            "radius_kpc": r["radius_kpc"],
            "required_extra_to_baryon_lower": need.tolist(),
            "excluded_stacks": {k: int((v < need).sum()) for k, v in budgets.items()},
        }
    return {
        "input_sha256": {
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (
                path,
                standard,
                maps,
                lens,
                sector_mass.pooled.sky.MAPS,
                HERE / "source_profile.json",
            )
        },
        "maximum_extra_to_ordinary": budgets,
        "imaged_profiles": profiles,
        "arbitrary_spherical_placement": {
            k: sector_home_mass.free_placement(inputs, masses, v, grav)
            for k, v in budgets.items()
        },
        "isolated_lens_stacks": stacks,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    out = run()
    dst = HERE / "sector_register_sky.json"
    if cli.parse_args().write:
        dst.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        inventory.rp.galaxy.check(json.loads(dst.read_text()), out)
    print(out["maximum_extra_to_ordinary"], out["imaged_profiles"]["scores"])
