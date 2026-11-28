"""Test the all-route radial acceleration and differential bearing on galaxy disks."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import collective_rotation as collective
from . import release_rotation as rotation

scales = collective.scales
HERE = scales.HERE


def run():
    lab = scales.budget.LAB / "landmark_release.json"
    record = json.loads(lab.read_bytes())["panels"]["1"]
    _, mu = scales.source_data()
    lib = collective.peers.inventory.rp.galaxy.transfer.budget
    obs, _, meta, _ = lib.load_inputs()
    masses = scales.pixel_maps()
    menus = list(masses)
    result, shape_tests, example = {}, {}, {}
    arrival = scales.budget.LAB / "landmark_arrival.json"
    profiles = record["profiles"] | {
        f"parked-{k}": v
        for k, v in json.loads(arrival.read_bytes())["parked_profiles"].items()
    }
    best_shape = {}
    for anchor, profile in profiles.items():
        shells = [s for s in profile["shells"] if s["radius"] > 0]
        r = np.array([s["radius"] for s in shells])
        # The annular integral projects a vector of magnitude r times the kernel.
        shape = np.r_[0.0, [s["inward"] / s["radius"] for s in shells], 0.0]
        knots = np.r_[0.0, r, r[-1] + profile["bin_width"]]
        gains = {
            s["bin"]: float(mu) / (s["radius"] ** 2 * s["ordinary"])
            for s in shells
            if s.get("ordinary", 0) > 0
        }
        strongest = max(shells, key=lambda s: s.get("ordinary", 0))["bin"]
        bounds = {menu: {} for menu in menus}
        for name, table in obs.items():
            radius, v, e, baryons = [
                np.array(table[k])
                for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
            ]
            low = np.maximum(0, v - 2 * e) ** 2 - baryons**2
            high = (v + 2 * e) ** 2 - baryons**2
            annuli = np.column_stack(
                [np.diff(np.r_[0, masses[m][name]]) for m in menus]
            )
            length = meta[name]["disk_scale_kpc"]
            projection = [
                lib.local.profiles.galaxies.G
                * radius[:, None]
                * (
                    collective.annular_force(radius, length, knots, shape, order)
                    @ annuli
                )
                for order in (12, 24)
            ]
            for j, menu in enumerate(menus):
                raw = scales.budget.passage_images.intervals(
                    projection[1][:, j], low, high
                )
                bounds[menu][name] = [(max(0.0, a), b) for a, b in raw]
            for k, gain in gains.items():
                coarse, fine = [v2 * gain for v2 in projection]
                ok = (low[:, None] <= fine) & (fine <= high[:, None])
                rough = (low[:, None] <= coarse) & (coarse <= high[:, None])
                for j, menu in enumerate(menus):
                    key = f"{anchor}/{k}/{menu}"
                    row = result.setdefault(
                        key, dict(within_two=0, complete=[], refinement_flips=0)
                    )
                    row["within_two"] += int(ok[:, j].sum())
                    row["refinement_flips"] += int((ok[:, j] != rough[:, j]).sum())
                    if ok[:, j].all():
                        row["complete"].append(name)
                    if k == strongest and menu == "mapped_baryons/1":
                        example.setdefault(anchor, {})[name] = dict(
                            radius_kpc=radius.tolist(),
                            excess_squared_speed=fine[:, j].tolist(),
                            within_two=int(ok[:, j].sum()),
                        )
        for menu, disks in bounds.items():
            best, _ = scales.budget.passage_images.best_gain(
                [band for bands in disks.values() for band in bands]
            )
            best_shape[f"{anchor}/{menu}"] = int(best)
            shape_tests[f"{anchor}/{menu}"] = {
                name: max(a for a, _ in bands) <= min(b for _, b in bands)
                for name, bands in disks.items()
            }
    metadata, catalogue = rotation.galaxies.catalogue()
    fixed = {}
    for name, rows in catalogue.items():
        good = [
            r for r in rows if min(r["rad_kpc"], r["vobs_kms"], r["e_vobs_kms"]) > 0
        ]
        if good:
            fixed[name] = rotation.zero_excess_fit(metadata[name], good)
    paths = [
        lab,
        arrival,
        HERE / "closed_inventory_maps.json",
        rotation.galaxies.sparc.DATA / "table2.dat",
        rotation.galaxies.sparc.DATA / "table1.dat",
    ]
    return dict(
        sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
        number_of_imaged_speeds=sum(len(v["radius_kpc"]) for v in obs.values()),
        strongest_ordinary_examples=example,
        ordinary_calibrations=result,
        whole_curve_possible_free_gain=shape_tests,
        most_bands_one_fitted_gain=best_shape,
        differential_landmarks={
            k: sum(t[k] for t in fixed.values())
            for k in ("rings", "within_two", "nominal_whole", "shared_parameters_whole")
        }
        | {"galaxies": len(fixed)},
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    raw = json.dumps(out, indent=2, allow_nan=False) + "\n"
    path = HERE / "landmark_rotation.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    fits = out["ordinary_calibrations"]
    print(
        len(fits),
        max(v["within_two"] for v in fits.values()),
        max(len(v["complete"]) for v in fits.values()),
        out["differential_landmarks"],
    )
