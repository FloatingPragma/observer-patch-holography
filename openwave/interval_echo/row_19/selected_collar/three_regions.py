"""Try a common positive local stress on three overlapping observer regions."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np

from . import midpoint_bulk as bulk

HERE = bulk.HERE
SEED = 28122833


def run():
    iv = bulk.shells.iv
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    src = bulk.shells.dist.read("resolved_shells")
    price = iv.mpf(src["conditional_entropy_slope"])
    word = [F(src["word_event_probabilities"][name]) for name in ("fresh", "shared")]
    specs = (
        ("half", 0.5, iv.mpf(3) / 4),
        ("boundary", 1.0, iv.sqrt(11) / 4),
        ("central", 0.0, iv.sqrt(15) / 4),
    )
    vertices, faces = bulk.base_chart()
    levels = []
    spreads = []
    for level in range(1, 5):
        vertices, faces = bulk.geometry._subdivide(vertices, faces)
        x, area, _ = bulk.sites(vertices, faces)
        spreads.append(float(np.max(np.abs(np.linalg.norm(x, axis=1) - 1))))
        weights = np.array([np.full(len(x), 1 / len(x)), area])
        panels = [{}, {}]
        for name, pos, radius in specs:
            rad = float(radius.mid)
            center = pos * x[0]
            pops = bulk.moments(x, center, weights, np.array([rad]))
            for k, law in enumerate(weights):
                trials, hits = bulk.direct_addresses(
                    x, law, center, np.array([rad]), rng
                )
                rows = {
                    "distinct_address_starts": trials,
                    "endpoint_hits": hits[0][0],
                    "weighted_midpoint_hits": hits[0][2],
                }
                for col, instrument in ((0, "endpoint"), (2, "weighted_midpoint")):
                    rows[instrument] = bulk.native_read(
                        F(float(pops[k, 0, col])), price, word, rng
                    )
                panels[k][name] = rows
        for k, row in enumerate(panels):
            for inst in ("endpoint", "weighted_midpoint"):
                residual = iv.mpf(0)
                for name, _, radius in specs:
                    q = iv.mpf(row[name][inst]["charge"])
                    if inst == "weighted_midpoint":
                        q *= 15 * radius / 32
                        row[name]["stationary_gain_times_area_forecast"] = (
                            bulk.shells.bounds(-q)
                        )
                    residual += (2 if name == "half" else -1) * radius * q
                row[inst + "_residual"] = bulk.shells.bounds(residual)
                if inst == "endpoint":
                    width = 2 * specs[0][2] + specs[1][2] + specs[2][2]
                    row["missing_regional_term_lower_bound"] = bulk.shells.bounds(
                        residual / width
                    )[0]
                    row["negative_mass_lower_bound"] = bulk.shells.bounds(
                        8 * residual / (7 * iv.pi)
                    )[0]
            levels.append(
                {
                    "level": level,
                    "stations": len(x),
                    "law": ("uniform_station", "solid_angle")[k],
                    **row,
                }
            )
    sweep_rng = np.random.default_rng(SEED + 1)
    signed = []
    for denom in (16, 256, 4096):
        assert max(spreads) + bulk.TOL < 1 / denom
        below = bulk.native_read(F(0), price, word, sweep_rng)
        above = bulk.native_read(F(1), price, word, sweep_rng)
        jump = iv.mpf(above["charge"]) - iv.mpf(below["charge"])
        tv = jump * denom / (4 * iv.pi)
        signed.append(
            {
                "halfwidth": str(F(1, denom)),
                "below": below,
                "above": above,
                "difference": bulk.shells.bounds(jump),
                "minimum_total_variation": bulk.shells.bounds(tv)[0],
                "minimum_in_midpoint_masses": bulk.shells.bounds(
                    tv / (15 * price / (32 * iv.pi))
                )[0],
            }
        )
    # Round isotropic endpoint probabilities. All central words are missed,
    # even the individual words. The boundary and intermediate caps differ.
    analytic = (
        2 * (iv.mpf(3) / 4) * 25 / 1024 - (iv.sqrt(11) / 4) * 121 / 4096
    ) * price
    return {
        "seed": SEED,
        "radial_normalization_errors": spreads,
        "signed_stress_sweep": signed,
        "starts_per_word_arm": bulk.shells.STARTS,
        "radius_squared": {"half": "9/16", "boundary": "11/16", "central": "15/16"},
        "round_endpoint_charge_residual": bulk.shells.bounds(analytic),
        "levels": levels,
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("THREE_REGIONS.md", "midpoint_bulk.json", "resolved_shells.json")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    z = run()
    target = HERE / "three_regions.json"
    if opt.write:
        target.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(target.read_bytes()), z)
    print(z["round_endpoint_charge_residual"])
    for row in z["levels"]:
        print(
            row["stations"],
            row["law"],
            row["endpoint_residual"],
            row["weighted_midpoint_residual"],
        )
