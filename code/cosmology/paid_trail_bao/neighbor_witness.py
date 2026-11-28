"""Predict local recovery access from independent neighboring native head records."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from . import witness_access as access

rp = access.rp
HERE = rp.HERE


def neighbors(field):
    pad = np.pad(field, 1)
    return (pad[:-2, 1:-1] + pad[2:, 1:-1] + pad[1:-1, :-2] + pad[1:-1, 2:]) / 4


def dilute_limit(bank):
    Q = rp.source.Q
    blind, informed = (bank["dilute"][k] for k in rp.MENUS)
    assert (blind["order"], informed["order"]) == (2, 1)
    flags, sizes = bank["third_magnitude_two_words"], bank["classes"]
    assert flags[0] == 0
    flag_slope = Q(4 * flags[1], sizes[1]["words"])
    local = -Q(blind["coefficient"])
    neighbor = Q(informed["coefficient"]) * flag_slope
    assert local > 0 and neighbor > 0
    return dict(
        local=str(local), neighbor=str(neighbor), threshold=str(local / neighbor)
    )


def acquire(cache):
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    prior = json.loads((HERE / "record_profile_maps.json").read_bytes())
    peak = prior["peak_msun_pc2"]
    bank = json.loads((rp.source.HERE / "record_density.json").read_bytes())
    rates = [
        a / b["words"]
        for a, b in zip(bank["third_magnitude_two_words"], bank["classes"], strict=True)
    ]
    assert rates[0] == 0
    profiles = {tag: {str(c): {} for c in rp.CAPACITY} for tag in rp.TRACERS}
    for name, tag, img, shell, selected in rp.observed_cells(cache, inputs):
        density = img[selected]
        nearness = float((density @ neighbors(img)[selected]) / (density @ density))
        for c in rp.CAPACITY:
            cap = c * peak[tag]
            raw = img / cap
            oversized = raw > 1
            flags = sum(
                w * rate
                for w, rate in zip(
                    rp.source.input_weights(np.minimum(raw, 1)), rates, strict=True
                )
            )
            fraction = neighbors(flags)[selected]
            assert min(fraction) >= 0 and max(fraction) <= 1
            endpoint = rp.source_laws(raw[selected], bank)
            change = endpoint["none"] + fraction * (
                endpoint["third_head"] - endpoint["none"]
            )
            count = len(inputs[name]["radius_kpc"])
            deposited = np.bincount(shell[selected], change, minlength=count)
            affected = neighbors(oversized.astype(float))[selected] > 0
            profiles[tag][str(c)][name] = dict(
                enclosed=(
                    deposited.cumsum() * rp.galaxy.WIDTH**2 * rp.maps.AREA_UNITS
                ).tolist(),
                access_range=[float(min(fraction)), float(max(fraction))],
                neighbor_to_self_product=nearness,
                saturated_neighbors_recipient_cells=int(affected.sum()),
                maximum_image_filling_excess=float(max(0, raw.max() - 1)),
            )
    return profiles


def run(profiles):
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    challenge = json.loads((HERE / "witness_access.json").read_bytes())
    radii = challenge["radius_kpc"]
    rational = access.rational
    indices = [inputs["NGC2841"]["radius_kpc"].index(r) for r in radii]
    bounds = [
        [rational(v) for v in challenge[k]]
        for k in ("mass_lower_msun", "mass_upper_msun")
    ]
    G = rp.galaxy.transfer.budget.local.profiles.galaxies.G
    tests = {}
    for tag, scales in profiles.items():
        tests[tag] = {}
        for c, rows in scales.items():
            f = {
                n: G * np.asarray(row["enclosed"]) / inputs[n]["radius_kpc"]
                for n, row in rows.items()
            }
            comparison = rp.scale_scores(inputs, f)
            enclosed = [rational(rows["NGC2841"]["enclosed"][i]) for i in indices]
            increments = [rational(0)] + [
                b - a for a, b in zip(enclosed[:-1], enclosed[1:], strict=True)
            ]
            shell = access.variable_access(radii, *bounds, increments, increments)
            comparison["ngc2841_with_free_interior"] = shell
            tests[tag][c] = comparison
    limit = dilute_limit(
        json.loads((rp.source.HERE / "record_density.json").read_bytes())
    )
    threshold = float(rp.source.Q(limit["threshold"]))
    limit["apertures"] = {
        tag: {
            name: {
                "neighbor_to_self_product": row["neighbor_to_self_product"],
                "leading_source_sign": int(
                    np.sign(threshold - row["neighbor_to_self_product"])
                ),
            }
            for name, row in scale["1"].items()
        }
        for tag, scale in profiles.items()
    }
    return dict(
        dilute_limit=limit,
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "NEIGHBOR_WITNESS.md",
                HERE / "neighbor_witness_maps.json",
                HERE / "record_profile_maps.json",
                HERE / "witness_access.json",
                rp.source.HERE / "record_density.json",
            )
        },
        tests=tests,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--image-cache", type=Path)
    cli.add_argument("--write", action="store_true")
    opts = cli.parse_args()
    src = HERE / "neighbor_witness_maps.json"
    if opts.image_cache:
        maps = acquire(opts.image_cache)
        if opts.write:
            src.write_text(json.dumps(maps, indent=2, allow_nan=False) + "\n")
        else:
            rp.galaxy.check(json.loads(src.read_bytes()), maps)
    result = run(json.loads(src.read_bytes()))
    dst = HERE / "neighbor_witness.json"
    if opts.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        rp.galaxy.check(json.loads(dst.read_bytes()), result)
    for tag, scales in result["tests"].items():
        print(
            tag,
            {
                c: (
                    row["within_two"],
                    row["held_full_curves"],
                    row["ngc2841_with_free_interior"]["possible"],
                )
                for c, row in scales.items()
            },
        )
