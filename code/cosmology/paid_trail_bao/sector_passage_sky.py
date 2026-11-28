"""Test uniform mean-clock rescaling and the calibrated exterior phase-clock source on SPARC."""

import argparse
import json
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256

import numpy as np

from . import release_rotation, sector_home_mass, sector_mass

inventory = sector_mass.inventory
HERE = inventory.HERE


def interval(baryon2, speed, error):
    lo = max(0, speed - 2 * error) ** 2
    hi = (speed + 2 * error) ** 2
    if baryon2 > 0:
        return lo / baryon2, hi / baryon2
    if lo == 0:
        return Q(0), None if baryon2 == 0 else Q(0)
    return None


def most(rows):
    events = defaultdict(lambda: [0, 0])
    for row in rows:
        if row is not None:
            lo, hi = row
            events[lo][0] += 1
            if hi is not None:
                events[hi][1] += 1
    live = best = 0
    gain = Q(0)
    for x, (begins, ends) in sorted(events.items()):
        live += begins
        if live > best:
            best, gain = live, x
        live -= ends
    return best, gain


def shape(groups):
    allrows = []
    result = {}
    for name, rows in groups.items():
        intervals = [interval(*r) for r in rows]
        best, gain = most(intervals)
        nominal = sum(
            z is not None and z[0] <= 1 and (z[1] is None or z[1] >= 1)
            for z in intervals
        )
        result[name] = {
            "speeds": len(rows),
            "nominal": nominal,
            "maximum": best,
            "squared_gain": str(gain),
        }
        if best < len(rows):
            invalid = next((i for i, r in enumerate(intervals) if r is None), None)
            if invalid is not None:
                result[name]["impossible_ring"] = invalid
            else:
                a = max(range(len(rows)), key=lambda i: intervals[i][0])
                b = min(
                    (i for i, z in enumerate(intervals) if z[1] is not None),
                    key=lambda i: intervals[i][1],
                )
                assert intervals[a][0] > intervals[b][1]
                result[name]["incompatible_rings"] = [a, b]
                result[name]["squared_gain_lower_upper"] = [
                    str(intervals[a][0]),
                    str(intervals[b][1]),
                ]
        allrows.extend(intervals)
    n, g = most(allrows)
    return {
        "galaxies": len(groups),
        "speeds": len(allrows),
        "at_ordinary_scale": sum(r["nominal"] for r in result.values()),
        "one_scale_per_galaxy": sum(r["maximum"] for r in result.values()),
        "complete": sum(r["maximum"] == r["speeds"] for r in result.values()),
        "one_scale_over_all": n,
        "common_squared_gain": str(g),
        "by_galaxy": result,
    }


def run():
    lab = inventory.LAB / "sector_passage.json"
    chart = HERE / "record_profile_maps.json"
    standard = inventory.LAB / "thermal_standard.json"
    receipt = json.loads(lab.read_text())
    source = receipt["spatial_modular_source"]
    inputs, _, _, galaxy_masses = inventory.rp.galaxy.transfer.budget.load_inputs()
    _, catalogue = release_rotation.galaxies.catalogue()
    groups = {}
    for name, rows in catalogue.items():
        valid = [
            r
            for r in rows
            if r["rad_kpc"] > 0 and r["vobs_kms"] > 0 and r["e_vobs_kms"] > 0
        ]
        groups[name] = []
        for row in valid:
            gas, disk, bulge = [
                Q(str(row[k])) for k in ("vgas_kms", "vdisk_kms", "vbul_kms")
            ]
            bary = gas * abs(gas) + disk * abs(disk) / 2 + Q(7, 10) * bulge * abs(bulge)
            groups[name].append(
                (bary, Q(str(row["vobs_kms"])), Q(str(row["e_vobs_kms"])))
            )
    imaged = {
        name: [
            (Q(str(b)) ** 2, Q(str(v)), Q(str(e)))
            for b, v, e in zip(
                d["baryons_kms"], d["observed_kms"], d["error_kms"], strict=True
            )
        ]
        for name, d in inputs.items()
    }
    nominal = shape(groups)
    assert (nominal["speeds"], nominal["at_ordinary_scale"]) == (3391, 452)
    raw = json.loads(standard.read_text())
    single, double = (raw[k]["energy_band"] for k in ("single", "double"))
    masses = {
        "single": [double[0] / single[1], double[1] / single[0]],
        "double": [1, 1],
    }
    G = inventory.rp.galaxy.transfer.budget.local.profiles.galaxies.G
    peaks = json.loads(chart.read_text())["peak_msun_pc2"]
    bbn = 3 * 0.1**2 / (8 * np.pi * G) * 0.02235
    bins = source["all_rates"]["bins"]
    profiles = {}
    for load, units in masses.items():
        panels = sector_mass.radial_channels(bins, 1.0, units, peaks, inputs, bbn)
        for method, curve in panels.items():
            cap = (
                max(row["rate"][1] * max(0, row["per_rate"][1]) for row in bins[method])
                / units[0]
            )
            # Grant the entire catalogued baryon budget to disjoint sources.
            # Their native reference has 37/7 ordinary single units per packet.
            ratio = cap / (37 / 28)
            curve["free_placement"] = sector_home_mass.free_placement(
                inputs, galaxy_masses, ratio, G
            )
            curve["maximum_extra_per_baryon_at_native_enrollment"] = (
                cap * (2 / 9) * (448 / 2110)
            )
        profiles[load] = panels
    paths = [
        lab,
        chart,
        standard,
        sector_mass.pooled.sky.MAPS,
        HERE / "source_profile.json",
        release_rotation.galaxies.sparc.DATA / "table2.dat",
    ]
    return {
        "sha256": {p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
        "serial_all_sparc": nominal,
        "serial_imaged": shape(imaged),
        "phase_clock_spherical_transfers": profiles,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "sector_passage_sky.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        inventory.rp.galaxy.check(json.loads(dst.read_text()), result)
    for name in ("serial_all_sparc", "serial_imaged"):
        print(name, {k: v for k, v in result[name].items() if k != "by_galaxy"})
