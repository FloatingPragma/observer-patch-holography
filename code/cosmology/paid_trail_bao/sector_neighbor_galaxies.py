"""Compare galaxy potential differences with mapped, azimuth-averaged round-trip clock responses."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import local_well
from . import sector_return_lapse as lapse

HERE = lapse.HERE
RATIOS = (1 / 16, 1 / 8, 1 / 4, 1 / 2, 1, 2, 4)
SCALES = ("pixel", *(str(r) for r in RATIOS))


def within(radius, source_radii, cutoff):
    """Fraction of a source ring within cutoff of a receiver at radius."""
    if radius == 0:
        return (source_radii < cutoff).astype(float)
    den = 2 * radius * source_radii
    q = np.divide(
        radius**2 + source_radii**2 - cutoff**2,
        den,
        out=np.zeros_like(den),
        where=den != 0,
    )
    return np.where(
        den > 0, np.arccos(np.clip(q, -1, 1)) / np.pi, float(radius < cutoff)
    )


def convolve(radii, source_radii, length, shape, amplitudes):
    # Integral of a shell step function via its jumps at half-hop boundaries.
    shape = np.asarray(shape)
    jumps = shape - np.r_[shape[1:], 0.0]
    grid = np.empty((len(radii), len(source_radii)))
    for i, r in enumerate(radii):
        grid[i] = sum(
            jump * within(r, source_radii, (k + 0.5) * length)
            for k, jump in enumerate(jumps)
            if jump
        )
    return (grid @ amplitudes).T


def differences(speeds_squared, radius):
    return (
        (speeds_squared[:-1] + speeds_squared[1:])
        / 2
        * np.diff(np.log(radius))
        / lapse.sky.C**2
    )


def feasible(diff, error, low, high):
    lo, hi = 0.0, np.inf
    for x, e, a, b in zip(diff, error, low, high, strict=True):
        l, h = lapse.gain_interval([x - e, x + e], [a, b])
        lo, hi = max(lo, l), min(hi, h)
    return bool(lo <= hi and hi > 0)


def run():
    inventory = lapse.sector_mass.inventory
    paths = [
        inventory.LAB / "sector_neighbor_ticks.json",
        HERE / "source_profile.json",
        lapse.sector_mass.pooled.sky.MAPS,
    ]
    shells = json.loads(paths[0].read_text())["shells"]
    data = json.loads(paths[1].read_text())["resolved_maps"]["inputs"]
    meta, _ = local_well.profiles.galaxies.catalogue()
    caps = np.array(inventory.rp.CAPACITY)
    targets, bounds = {}, {}
    fields = {}
    with np.load(paths[2]) as maps:
        for name, rows in data.items():
            r, v, e, b = (
                np.array(rows[k])
                for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
            )
            targets[name] = differences(b**2, r)
            bounds[name] = (
                differences(np.maximum(0, v - 2 * e) ** 2 - b**2, r),
                differences((v + 2 * e) ** 2 - b**2, r),
            )
            for tag in ("mapped_baryons", "atomic"):
                rp, idx = np.unique(
                    maps[f"{name}_{tag}_radius_kpc"], return_inverse=True
                )
                p = maps[f"{name}_{tag}"][:, None] / caps
                assert p.min() >= 0 and p.max() <= 1 + 1e-12
                strengths = {
                    "ordinary": 4 * p,
                    "cross": 3 * p**2 * (1 - p) ** 2 * 448 / 2110,
                }
                for key, strength in strengths.items():
                    amp = np.array(
                        [
                            np.bincount(idx, weights=col, minlength=len(rp))
                            for col in strength.T
                        ]
                    ).T
                    lengths = [
                        inventory.rp.galaxy.WIDTH,
                        *(meta[name]["disk_scale_kpc"] * x for x in RATIOS),
                    ]
                    for si, length in enumerate(lengths):
                        sh = shells[key]
                        val = convolve(r, rp, length, sh["log_rate"], amp)
                        # At most the largest calibration error per contributing tile,
                        # even across both endpoints of a potential difference.
                        band = 2 * np.max(sh["error"]) * amp.sum(0)
                        for ci, cap in enumerate(caps):
                            fields[name, tag, int(cap), si, key] = (
                                np.diff(val[ci]),
                                float(band[ci]),
                            )
    settings = {}
    for tag in ("mapped_baryons", "atomic"):
        for cap in caps:
            cap = int(cap)
            sub = {
                name: np.array(
                    [
                        fields[name, tag, cap, i, "ordinary"][0]
                        for i in range(len(SCALES))
                    ]
                )
                for name in data
            }
            normal = {
                n: float(np.linalg.norm(t) / np.sqrt(len(t)))
                for n, t in targets.items()
            }
            scores = {n: f / normal[n] for n, f in sub.items()}
            unit = {n: t / normal[n] for n, t in targets.items()}
            possible = np.zeros(len(SCALES), int)
            results = {}
            for n in data:
                lo, hi = bounds[n]
                train = sorted(set(data) - {n})
                i, g = local_well.calibrate(scores, unit, train)
                q, err = fields[n, tag, cap, i, "cross"]
                err = np.full_like(q, err)
                admitted = (g * (q + err) >= lo) & (g * (q - err) <= hi)
                for si in range(len(SCALES)):
                    x, w = fields[n, tag, cap, si, "cross"]
                    possible[si] += feasible(x, np.full_like(x, w), lo, hi)
                results[n] = {
                    "held_scale": SCALES[i],
                    "ordinary_gain": g,
                    "ordinary_relative_rms": float(
                        np.sqrt(np.mean((g * sub[n][i] - targets[n]) ** 2)) / normal[n]
                    ),
                    "excess_annuli": len(q),
                    "within_two": int(admitted.sum()),
                    "model_extra": (g * q).tolist(),
                    "model_numeric_halfwidth": float(g * err[0]),
                    "target_extra_lo": lo.tolist(),
                    "target_extra_hi": hi.tolist(),
                }
            settings[f"{tag}/{cap}"] = {
                "fully_admitted": sum(
                    s["within_two"] == s["excess_annuli"] for s in results.values()
                ),
                "all_annuli": sum(s["excess_annuli"] for s in results.values()),
                "admitted_annuli": sum(s["within_two"] for s in results.values()),
                "feasible_individual_gains": possible.tolist(),
                "galaxies": results,
            }
    return {
        "input_sha256": {p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
        "scales": list(SCALES),
        "settings": settings,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    dst = HERE / "sector_neighbor_galaxies.json"
    data = run()
    if cli.parse_args().write:
        dst.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        lapse.sector_mass.inventory.rp.galaxy.check(json.loads(dst.read_text()), data)
    print(
        {
            k: (
                v["fully_admitted"],
                v["admitted_annuli"],
                v["feasible_individual_gains"],
            )
            for k, v in data["settings"].items()
        }
    )
