"""Forecast the imaged and cosmic mass from the independently measured sector charge."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np

from . import closed_inventory as inventory
from . import pooled_ordinary as pooled

HERE = inventory.HERE


def temperature_ceiling(classes):
    # A centered energy line is always at least -33. Check the exact
    # rational factors, including the vertex of each quadratic.
    found = False
    for shape in classes:
        for factor in shape["centered_characteristic_factors"]:
            c = list(map(Q, factor["polynomial"]))
            if len(c) == 2:
                low = -c[1] / c[0]
                assert low >= -33
                found |= low == -33
            else:
                a, b, d = c
                assert a > 0 and -b / (2 * a) >= -33 and a * 33**2 - 33 * b + d >= 0
    assert found
    # Jensen gives 0 <= beta^-1 log mean exp(-beta*(E-Ebar)) <= 33
    # for every beta>0, independently of source or temperature.
    return 33


def priced_filling(filling, bins, beta, mass):
    ends = np.array([row["rate"][1] for row in bins])
    idx = np.searchsorted(ends, filling)
    limits = (
        np.array([row["per_rate"] for row in bins])[idx]
        * np.asarray(filling)[..., None]
        / beta
    )
    values = limits[..., None] / np.asarray(mass)
    return np.stack((values.min(axis=(-1, -2)), values.max(axis=(-1, -2))), axis=-1)


def radial_prediction(history, source, double, peaks, inputs, bbn):
    channels = {
        tag: history[method]["rate_bins"]
        for tag, method in (
            ("full", "continuous_cadence"),
            ("heat_subtracted", "continuous_sector_excess"),
        )
    }
    return radial_channels(channels, source["beta"], double, peaks, inputs, bbn)


def radial_channels(channels, signal_per_energy, double, peaks, inputs, bbn):
    result = {}
    lib = inventory.rp.galaxy.transfer.budget
    G = lib.local.profiles.galaxies.G
    width = inventory.rp.galaxy.WIDTH
    for tag, bins in channels.items():
        scores = {}
        profiles = {}
        cosmic = []
        with np.load(pooled.sky.MAPS) as charts:
            for tracer, peak in peaks.items():
                for cap in inventory.rp.CAPACITY:
                    pbar = bbn * width / (1e6 * peak * cap)
                    priced = priced_filling(pbar, bins, signal_per_energy, double)
                    cosmic.extend(priced * 1.5 * pbar * (1 - pbar) ** 2 * 448 / 2110)
                    key = f"{tracer}/{cap}"
                    row = dict(possible=0, complete=[])
                    for name, tab in inputs.items():
                        fill = charts[f"{name}_{tracer}"] / cap
                        radius = np.array(tab["radius_kpc"])
                        bins_r = np.searchsorted(
                            radius, charts[f"{name}_{tracer}_radius_kpc"]
                        )
                        weight = 6 * fill**2 * (1 - fill) ** 2 * 448 / 2110
                        units = priced_filling(fill, bins, signal_per_energy, double)
                        mass = weight[:, None] * units * width**2 * 1e6 * peak * cap / 4
                        enclosed = np.array(
                            [
                                np.bincount(
                                    bins_r, weights=v, minlength=len(radius)
                                ).cumsum()
                                for v in mass.T
                            ]
                        )
                        v2 = np.array(tab["baryons_kms"]) ** 2 + G * enclosed / radius
                        assert min(v2.ravel()) >= 0
                        speed = np.sqrt(v2)
                        vel, err = (
                            np.array(tab["observed_kms"]),
                            np.array(tab["error_kms"]),
                        )
                        possible = (speed[1] >= vel - 2 * err) & (
                            speed[0] <= vel + 2 * err
                        )
                        row["possible"] += int(possible.sum())
                        if possible.all():
                            row["complete"].append(name)
                        if tracer == "mapped_baryons" and cap == 1:
                            profiles[name] = dict(
                                radius_kpc=radius.tolist(),
                                enclosed_extra_mass_msun=enclosed.T.tolist(),
                                predicted_speed_kms=speed.T.tolist(),
                            )
                    scores[key] = row
        result[tag] = dict(
            scores=scores,
            mapped_baryon_profiles=profiles,
            bbn_ratio_bounds=[min(cosmic), max(cosmic)],
        )
    return result


def run():
    native = inventory.LAB / "sector_calorimetry.json"
    clock = inventory.LAB / "thermal_standard.json"
    image = HERE / "closed_inventory_maps.json"
    profiles = HERE / "record_profile_maps.json"
    source = json.loads(native.read_bytes())
    double = json.loads(clock.read_bytes())["double"]["energy_band"]
    maps = json.loads(image.read_bytes())
    prior = json.loads((inventory.LAB / "peer_register_lab.json").read_bytes())
    incidence = float(Q(prior["passive_probability"]))
    inputs, *_ = inventory.rp.galaxy.transfer.budget.load_inputs()
    history = inventory.LAB / "sector_retention.json"
    delayed = json.loads(history.read_bytes())
    cap = temperature_ceiling(source["spectral_classes"])
    experiments = {
        "registered_source": source["full_equivalent_energy"],
        "any_age_full_read": [
            v / source["beta"]
            for v in delayed["all_word_independent_age_mixtures_nats"]
        ],
        "any_age_after_heat_subtraction": [
            v / source["beta"] for v in delayed["all_age_sector_excess_nats"]
        ],
        "one_passive_register": [
            incidence * x for x in source["full_equivalent_energy"]
        ],
        "any_matched_mixture_at_thermometer": source["matched_compact_energy_envelope"],
        "any_matched_mixture_at_any_positive_temperature": [-cap, cap],
    }
    out = {}
    for name, energy in experiments.items():
        scaled = [u / v for u in energy for v in double]
        limits = [
            float(np.nextafter(min(scaled), -np.inf)),
            float(np.nextafter(max(scaled), np.inf)),
        ]
        raw = inventory.mass.forecast(maps, limits, inputs)
        scores = {}
        for tag, v in raw.items():
            for c, s in v.items():
                row = s["closed"]
                key = f"{tag}/{c}"
                scores[key] = dict(
                    possible=row["within_two"],
                    complete=[
                        g
                        for g, d in row["galaxies"].items()
                        if d["within_two"] == len(inputs[g]["radius_kpc"])
                    ],
                    maximum_speed_change_kms=row["maximum_baryon_speed_change_kms"],
                )
        out[name] = dict(ordinary_double_units=limits, scores=scores)
    peak = json.loads(profiles.read_bytes())["peak_msun_pc2"]
    lib = inventory.rp.galaxy.transfer.budget
    G = lib.local.profiles.galaxies.G
    width = inventory.rp.galaxy.WIDTH
    bbn = 3 * 0.1**2 / (8 * np.pi * G) * 0.02235
    # 6p²(1-p)² source opportunities per four-port cell, containing
    # 4p ordinary doubles on average, followed by native enrollment.
    for panel in out.values():
        limits = panel["ordinary_double_units"]
        panel["maximum_positive_ratio_any_filling"] = (
            max(0, limits[1]) * (2 / 9) * (448 / 2110)
        )
        panel["maximum_positive_ratio_one_source_per_two_doubles"] = (
            max(0, limits[1]) / 2
        )
        bbn_range = []
        for dens in peak.values():
            for c in inventory.rp.CAPACITY:
                p = bbn * width / (1e6 * dens * c)
                frac = 1.5 * p * (1 - p) ** 2 * 448 / 2110
                bbn_range.extend(frac * np.array(limits))
        panel["bbn_ratio_envelope"] = [min(bbn_range), max(bbn_range)]
    radial = radial_prediction(delayed, source, double, peak, inputs, bbn)
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                native,
                clock,
                image,
                profiles,
                history,
                pooled.sky.MAPS,
                inventory.LAB / "peer_register_lab.json",
            )
        },
        reporting_probability_equals_imaged_filling=radial,
        adopted_base_lcdm_cold_to_baryons=0.1200 / 0.02237,
        imaged_galaxies=len(inputs),
        speeds=sum(len(v["radius_kpc"]) for v in inputs.values()),
        transfers=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "sector_mass.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        inventory.rp.galaxy.check(json.loads(path.read_text()), out)
    print(
        {
            k: (
                max(y["possible"] for y in v["scores"].values()),
                v["maximum_positive_ratio_any_filling"],
            )
            for k, v in out["transfers"].items()
        }
    )
