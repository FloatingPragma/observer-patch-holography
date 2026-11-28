"""Apply the measured multi-source tides and ordinary controls to the enrolled disks."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np

from . import tidal_neighbors as peers
from . import tidal_scales as scales
from .tidal_scale_kernel import annular_force

HERE = scales.HERE


def profiles(lab):
    out = {}
    for name, panel in lab["carriers"].items():
        if not panel["enrolled"]:
            continue
        m = int(name.split("/")[-1])
        if m == 1 and "dispersed" in name:
            continue
        shells = {int(k): v for k, v in panel["shells"].items() if v["chords"]}
        last = max(shells)
        xs = [0, *sorted(shells), last + 1]
        shape = np.interp(
            np.arange(last + 2),
            xs,
            [
                0.0,
                *[-float(Q(shells[k]["average_change"])) for k in sorted(shells)],
                0.0,
            ],
        )
        ordinary = {
            k: -Q(v["average_ordinary"])
            for k, v in shells.items()
            if Q(v["average_ordinary"]) < 0
        }
        assert np.all(shape >= 0)
        out[name] = (
            shape,
            ordinary,
            dict(
                missing=[k for k in range(1, last + 1) if k not in shells],
                measured_shells=sorted(shells),
                nearest_source_ranges=panel["nearest_source_ranges"],
                ordinary_tide={str(k): str(v) for k, v in ordinary.items()},
                peak_ordinary_shell=max(ordinary, key=ordinary.get),
            ),
        )
    return out


def image_terms(name, r, disk_scale, raw, stock, menus):
    # Each entry specifies a hop and the imaged source mass residing in
    # each annulus for every menu. The two instruments have different scales.
    rc = np.column_stack([np.diff(np.r_[0, stock[menu][name]]) for menu in menus])
    terms = {"stellar_scale": [(disk_scale, rc)], "six_hop_peer": []}
    all_bins = {}
    for j, menu in enumerate(menus):
        tag, cap = menu.split("/")
        bins = raw[f"{name}_{tag}_bins"]
        values = raw[f"{name}_{tag}_{cap}"]
        for k, row in zip(bins, values[:-1], strict=True):
            if k not in all_bins:
                all_bins[k] = np.zeros_like(rc)
            all_bins[k][:, j] = row
    for k, arr in sorted(all_bins.items()):
        terms["six_hop_peer"].append((peers.WIDTH * np.sqrt(k) / 6, arr))
    return terms


def hat_fields(r, terms, active, order):
    # Kernel interpolation is linear. Integrate each triangular basis shape
    # once per physical length and reuse it for all measured carrier profiles.
    result = {}
    for mode, rows in terms.items():
        moments = []
        for center in active:
            total = np.zeros_like(rows[0][1])
            for hop, stock in rows:
                if (center - 1) * hop >= 2 * r[-1]:
                    continue
                area = annular_force(
                    r, hop, [center - 1, center, center + 1], [0.0, 1.0, 0.0], order
                )
                total += area @ stock
            moments.append(total)
        result[mode] = np.array(moments)
    return result


def median_matched(r, shape, source_ranges, raw, name, menus):
    native = np.median(source_ranges)
    knots = np.arange(len(shape))
    fine = np.zeros((len(r), len(menus)))
    coarse = np.zeros_like(fine)
    metadata = {}
    areas = {}
    for j, menu in enumerate(menus):
        tag, cap = menu.split("/")
        bins = raw[f"{name}_{tag}_bins"]
        masses = raw[f"{name}_{tag}_{cap}"][:-1]
        by_size = masses.sum(1)
        assert by_size.sum() > 0
        k = int(bins[np.searchsorted(by_size.cumsum(), by_size.sum() / 2)])
        hop = peers.WIDTH * np.sqrt(k) / native
        if k not in areas:
            areas[k] = [annular_force(r, hop, knots, shape, order) for order in (8, 16)]
        annular = masses.sum(0)
        coarse[:, j], fine[:, j] = [area @ annular for area in areas[k]]
        metadata[menu] = float(hop)
    return (coarse, fine), metadata


def run():
    path = scales.budget.LAB / "collective_tide.json"
    survey = profiles(json.loads(path.read_bytes()))
    _, mu = scales.source_data()
    lib = peers.inventory.rp.galaxy.transfer.budget
    obs, _, meta, _ = lib.load_inputs()
    G = lib.local.profiles.galaxies.G
    stock = scales.pixel_maps()
    menus = list(stock)
    active = sorted(
        set.union(*(set(map(int, np.flatnonzero(v[0]))) for v in survey.values()))
    )
    with np.load(peers.MAPS) as data:
        raw = dict(data)
    result = {}
    cases = {}
    measured_lengths = {}
    shape_bands = {}
    for name, tab in obs.items():
        r, v, e, b = [
            np.asarray(tab[k])
            for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
        ]
        low, hi = np.maximum(v - 2 * e, 0) ** 2 - b**2, (v + 2 * e) ** 2 - b**2
        terms = image_terms(name, r, meta[name]["disk_scale_kpc"], raw, stock, menus)
        coarse, fine = (hat_fields(r, terms, active, order) for order in (8, 16))
        for label, (shape, refs, info) in survey.items():
            values = np.array([shape[k] if k < len(shape) else 0.0 for k in active])
            pair_fields = {
                mode: tuple(
                    np.tensordot(values, tables[mode], axes=(0, 0))
                    for tables in (coarse, fine)
                )
                for mode in terms
            }
            if info["nearest_source_ranges"]:
                field, lengths = median_matched(
                    r, shape, info["nearest_source_ranges"], raw, name, menus
                )
                pair_fields["matched_peer_median"] = field
                measured_lengths.setdefault(label, {})[name] = lengths
            for mode, (co, project) in pair_fields.items():
                project = G * r[:, None] * project
                co = G * r[:, None] * co
                for j, menu in enumerate(menus):
                    bands = scales.budget.passage_images.intervals(
                        project[:, j], low, hi
                    )
                    shape_bands.setdefault((label, mode, menu), {})[name] = [
                        (max(0, a), b) for a, b in bands
                    ]
                for k, ref in refs.items():
                    gain = float(mu / (k**3 * ref))
                    prediction = project * gain
                    hit = (low[:, None] <= prediction) & (prediction <= hi[:, None])
                    other = (low[:, None] <= co * gain) & (co * gain <= hi[:, None])
                    for j, menu in enumerate(menus):
                        key = f"{label}/{mode}/{k}/{menu}"
                        stats = result.setdefault(
                            key,
                            dict(within_two=0, complete=[], quadrature_disagreements=0),
                        )
                        stats["within_two"] += int(hit[:, j].sum())
                        if hit[:, j].all():
                            stats["complete"].append(name)
                        stats["quadrature_disagreements"] += int(
                            (hit[:, j] != other[:, j]).sum()
                        )
                        if (
                            menu == "mapped_baryons/1"
                            and k == info["peak_ordinary_shell"]
                        ):
                            cases.setdefault(f"{label}/{mode}", {})[name] = dict(
                                radius_kpc=r.tolist(),
                                excess_squared_speed=prediction[:, j].tolist(),
                                within_two=int(hit[:, j].sum()),
                            )
    shape_result = {}
    for (label, mode, menu), panels in shape_bands.items():
        bands = [v for row in panels.values() for v in row]
        count, _ = scales.budget.passage_images.best_gain(bands)
        shape_result[f"{label}/{mode}/{menu}"] = dict(
            common_amplitude_rings=int(count),
            complete_with_galaxy_amplitude=[
                n
                for n, row in panels.items()
                if max(a for a, b in row) <= min(b for a, b in row)
            ],
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                path,
                HERE / "COLLECTIVE_ROTATION.md",
                peers.MAPS,
                HERE / "source_profile.json",
                HERE / "closed_inventory_maps.json",
                scales.budget.LAB / "moving_records.json",
            )
        },
        profiles={k: v[2] for k, v in survey.items()},
        measured_kpc_per_hop=measured_lengths,
        scores=result,
        shape_only=shape_result,
        mapped_baryons_first_capacity_peak_ordinary=cases,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "collective_rotation.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        peers.inventory.rp.galaxy.check(json.loads(dst.read_bytes()), result)
    print(
        {
            mode: max(
                v["within_two"] for key, v in result["scores"].items() if mode in key
            )
            for mode in ("stellar_scale", "six_hop_peer", "matched_peer_median")
        }
    )
