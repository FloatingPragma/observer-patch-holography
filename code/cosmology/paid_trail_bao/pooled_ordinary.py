"""Predict ordinary-equivalent mass distributions from conditional native relay matching."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np

from cosmology.rar_deep_regime.boundary_source.stellar_maps import ordinary_disk_v2
from . import collective_occupancy as sky
from . import tidal_scales as scales

HERE = scales.HERE
PORTS = 48
WIDTH = sky.neighbor.WIDTH


def expected_per_tile(w, group, coefficients):
    sizes = np.bincount(group)
    live = sizes > 0
    lam = np.bincount(group, weights=w)
    log_miss = np.bincount(group, weights=np.log1p(-w))
    a, b, c = coefficients
    mean = (a * PORTS * sizes + c) * (-np.expm1(log_miss)) + b * lam
    assert np.all(mean >= 0)
    result = np.zeros(len(sizes))
    result[live] = mean[live] / sizes[live]
    return result[group]


def field_from_tiles(r, rings, rr, unique, mass, G):
    shell = np.bincount(rings, weights=mass, minlength=len(r))
    ring_mass = np.bincount(rr, weights=mass)
    fields = {"spherical": G * shell.cumsum() / r}
    for h in (WIDTH / 2, WIDTH):
        fields[f"disk_height_{h:g}"] = np.array(
            [ordinary_disk_v2(radius, ring_mass, unique, h) for radius in r]
        )
    return shell, fields


def homogeneous(n, p, coefficients, mu, acceptance=1.0):
    a, b, c = coefficients
    w = 6 * p**2 * (1 - p) ** 2 * (448 / 2110) * acceptance
    if n == math.inf:
        return mu * (a * PORTS + b * w) / (4 * p)
    any_source = -math.expm1(n * math.log1p(-w))
    return mu * ((a * PORTS + c / n) * any_source + b * w) / (4 * p)


def cosmic(peaks, coeffs, mu, G, acceptance=1.0):
    # H100 is 0.1 km/s/kpc. BBN fixes physical density independently
    # of the value assigned to H0. The surface density is in Msun/pc².
    omega_b = 0.02235
    rho = 3 * 0.1**2 / (8 * np.pi * G) * omega_b
    target = 0.1200 / 0.02237
    out = {}
    for tag, peak in peaks.items():
        for cap in sky.neighbor.inventory.rp.CAPACITY:
            p = rho * WIDTH / (1e6 * peak * cap)
            ratios = {
                str(n): homogeneous(n, p, coeffs, mu, acceptance)
                for n in (1, 5, 20, 80)
            }
            limit = homogeneous(math.inf, p, coeffs, mu, acceptance)
            # On this input the curve rises across the observed crossing.
            # Verify the bracket; the inferred scale is a sensitivity diagnostic.
            lo, hi = 1, 2
            while homogeneous(hi, p, coeffs, mu, acceptance) < target:
                lo, hi = hi, 2 * hi
            while hi - lo > 1:
                mid = (hi + lo) // 2
                if homogeneous(mid, p, coeffs, mu, acceptance) >= target:
                    hi = mid
                else:
                    lo = mid
            assert (
                homogeneous(lo, p, coeffs, mu, acceptance)
                < target
                <= homogeneous(hi, p, coeffs, mu, acceptance)
            )
            out[f"{tag}/{cap}"] = dict(
                filling=p,
                ratios_by_tiles=ratios,
                infinite_pool_ratio=limit,
                inferred_tiles=hi,
                inferred_side_kpc=WIDTH * hi ** (1 / 3),
            )
    return dict(
        bbn_omega_b=omega_b,
        present_baryons_msun_kpc3=rho,
        planck_base_lcdm_cold_to_baryon=target,
        tests=out,
    )


def run(receipt=None, acceptance=1.0):
    if receipt is None:
        receipt = scales.budget.LAB / "collective_relay_lab.json"
    coeff = list(
        map(
            lambda x: float(Q(x)),
            json.loads(receipt.read_bytes())["ordinary_coefficients"],
        )
    )
    _, mu = scales.source_data()
    lib = sky.neighbor.inventory.rp.galaxy.transfer.budget
    G = lib.local.profiles.galaxies.G
    inputs, *_ = lib.load_inputs()
    peaks = json.loads((HERE / "record_profile_maps.json").read_bytes())[
        "peak_msun_pc2"
    ]
    scores = {}
    selected = {}
    ratios = {}
    with np.load(sky.MAPS) as raw:
        for name, tab in inputs.items():
            r, v, e, b = [
                np.array(tab[k])
                for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
            ]
            for tag, peak in peaks.items():
                density = raw[f"{name}_{tag}"]
                actual_radius = raw[f"{name}_{tag}_radius_kpc"]
                unique, rr = np.unique(actual_radius, return_inverse=True)
                rings = np.searchsorted(r, actual_radius)
                for cap in sky.neighbor.inventory.rp.CAPACITY:
                    p = density / cap
                    w = 6 * p**2 * (1 - p) ** 2 * (448 / 2110) * acceptance
                    per_double = WIDTH**2 * 1e6 * peak * cap / 4
                    for pool, group in (
                        ("tile", np.arange(len(p))),
                        ("annulus", rings),
                        ("disk", np.zeros(len(p), int)),
                    ):
                        donors = expected_per_tile(w, group, coeff)
                        mass = donors * float(mu) * per_double
                        ratio = float(mass.sum() / (4 * per_double * p.sum()))
                        tag2 = f"{pool}/{tag}/{cap}"
                        ratios.setdefault(tag2, {})[name] = ratio
                        shell, fields = field_from_tiles(r, rings, rr, unique, mass, G)
                        for mode, extra in fields.items():
                            total = b * b + extra
                            hit = (np.maximum(v - 2 * e, 0) ** 2 <= total) & (
                                total <= (v + 2 * e) ** 2
                            )
                            key = f"{tag2}/{mode}"
                            row = scores.setdefault(
                                key, dict(within_two=0, complete=[], undefined=0)
                            )
                            row["within_two"] += int(hit.sum())
                            row["undefined"] += int((total < 0).sum())
                            if hit.all():
                                row["complete"].append(name)
                            if tag == "mapped_baryons" and cap == 1:
                                selected.setdefault(f"{pool}/{mode}", {})[name] = dict(
                                    radius_kpc=r.tolist(),
                                    predicted_excess_v2=extra.tolist(),
                                    equivalent_shell_msun=shell.tolist(),
                                    native_equivalent_to_enrolled_baryons=ratio,
                                    within_two=int(hit.sum()),
                                )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                receipt,
                sky.MAPS,
                HERE / "source_profile.json",
                HERE / "record_profile_maps.json",
                HERE / "POOLED_ORDINARY.md",
            )
        },
        homogeneous=cosmic(peaks, coeff, float(mu), G, acceptance),
        scores=scores,
        ratios=ratios,
        mapped_baryons_first_capacity=selected,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "pooled_ordinary.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        sky.neighbor.inventory.rp.galaxy.check(json.loads(dst.read_bytes()), result)
    print(
        {
            pool: max(
                v["within_two"]
                for k, v in result["scores"].items()
                if k.startswith(pool)
            )
            for pool in ("tile", "annulus", "disk")
        }
    )
