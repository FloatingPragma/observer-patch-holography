"""Read the same surface source through the distance obtained by two angular deliveries."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import cos, pi, tan

import numpy as np

from . import midpoint_bulk as bulk

HERE = bulk.HERE
SEED = 28122863


def run():
    iv = bulk.shells.iv
    iv.dps = 35
    source = bulk.shells.dist.read("resolved_shells")
    price = iv.mpf(source["conditional_entropy_slope"])
    word = [F(source["word_event_probabilities"][a]) for a in ("fresh", "shared")]
    rng = np.random.default_rng(SEED)
    vertices, faces = bulk.base_chart()
    data = []
    for level in range(1, 5):
        vertices, faces = bulk.geometry._subdivide(vertices, faces)
        x, weights, mesh = bulk.sites(vertices, faces)
        projection = x @ x[0]
        distinct = 1 - weights @ weights
        for R in (pi / 6, pi / 4, pi / 3):
            c = cos(R)
            mask = projection >= c - 1e-13
            end = ((weights @ mask) ** 2 - weights[mask] @ weights[mask]) / distinct
            moments = np.zeros(3)
            for a in range(0, len(x), 64):
                z = (projection[a : a + 64, None] + projection[None, :]) / 2
                d = np.arccos(np.clip(z, -1, 1))
                block = weights[a : a + 64, None] * weights[None, :]
                block[np.arange(len(z)), a + np.arange(len(z))] = 0
                moments += (
                    np.array(
                        [
                            np.sum(block * v)
                            for v in (
                                z >= c - 1e-13,
                                np.maximum(1 - (d / R) ** 2, 0),
                                np.maximum((z - c) / (1 - c), 0),
                            )
                        ]
                    )
                    / distinct
                )
            row = {}
            for key, v in zip(
                ("endpoint", "enclosed", "flat_weight", "cosine_weight"),
                [end, *moments],
                strict=True,
            ):
                row[key] = bulk.native_read(F(float(v)), price, word, rng)
            q = iv.mpf(row["endpoint"]["charge"])
            rr = iv.mpf(R)
            row["mass_over_price_flat"] = bulk.shells.bounds(
                q / (iv.pi * rr * iv.mpf(row["flat_weight"]["charge"]))
            )
            row["mass_over_price_cosine"] = bulk.shells.bounds(
                q
                / (2 * iv.pi * iv.tan(rr / 2) * iv.mpf(row["cosine_weight"]["charge"]))
            )
            data.append(
                {
                    "level": level,
                    "stations": len(x),
                    "radius": R,
                    "mesh_bound": mesh,
                    "rows": row,
                    "isotropic_endpoint": (1 - c) ** 2 / 4,
                    "isotropic_enclosed": (1 - c) ** 2 / 2,
                    "isotropic_cosine_mass_over_price": 3 / (4 * pi * tan(R / 2)),
                }
            )
    return {
        "seed": SEED,
        "measurements": data,
        "inputs": {
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in ("PAIRED_RANGING.md", "resolved_shells.json", "midpoint_bulk.json")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    z = run()
    target = HERE / "hemisphere_source.json"
    if args.write:
        target.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(target.read_bytes()), z)
    print(
        [
            (
                r["radius"],
                r["rows"]["mass_over_price_flat"],
                r["rows"]["mass_over_price_cosine"],
            )
            for r in z["measurements"][-3:]
        ]
    )
