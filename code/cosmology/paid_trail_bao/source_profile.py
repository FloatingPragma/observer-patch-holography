"""Forward galaxy profiles from measured modular shell increments."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import measure as bao
from cosmology.rar_deep_regime.boundary_source import source as maps
from cosmology.rar_deep_regime.individual_trail_budget import measure as galaxies
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    source_reader,
)

RATIOS = 2.0 ** (np.arange(-64, 33) / 8)
CAP = galaxies.native_caps()["32"]["cap"]


def annulus_fraction(radius, shell, a, b):
    """Spherical enclosure of an isotropic write about a seed uniform in annular area."""
    # For a seed at s the allowed direction cosine occupies the fraction
    # clip((R²-(s-h)²)/(4sh),0,1). Integrate against 2s ds/(b²-a²).
    if np.any(np.asarray(shell) <= 0) or np.any(b <= a):
        raise ValueError("A positive shell radius and nonempty annulus are required")
    full = np.maximum(0, np.minimum(b, np.maximum(0, radius - shell)) ** 2 - a**2)
    lower = np.maximum(a, abs(radius - shell))
    upper = np.minimum(b, radius + shell)
    width = np.maximum(0, upper - lower)
    mid = (lower + upper) / 2
    partial = width * (radius**2 - (mid - shell) ** 2 - width**2 / 12) / (2 * shell)
    value = (full + partial) / (b * b - a * a)
    assert np.min(value) > -1e-9 and np.max(value) < 1 + 1e-9
    return np.clip(value, 0, 1)


def observed_shapes(path):
    labs = json.loads(path.read_bytes())["cases"]
    shapes = {}
    for name in ("outward", "folded"):
        tape = np.array(labs["8_A"]["decoded"][name]["shells"])
        active = tape[:, 1] > 0
        radii, charge = tape[active, :2].T
        shapes[name] = radii, charge / charge.sum()
    return shapes


def mapped_annuli(kept):
    out = {}
    for name, (meta, survey, rows, _, _) in kept.items():
        s = (
            np.array([float(row["r"]) for row in rows])
            * meta["dist_mpc"]
            / float(survey["Dist"])
        )
        edges = maps.edges(s)
        surface = sum(
            maps.filled(s, rows, key)[0] for key in ("Sigma*", "SigmaHI", "SigmaH2")
        )
        mass = surface * np.diff(edges**2)
        assert min(mass) >= 0 and sum(mass) > 0
        out[name] = edges, mass / sum(mass)
    return out


def bare_speed2(rows):
    return np.array(
        [
            sum(
                weight * galaxies.signed_square(row[k])
                for k, weight in (
                    ("vdisk_kms", 0.5),
                    ("vbul_kms", 0.7),
                    ("vgas_kms", 1),
                )
            )
            for row in rows
        ]
    )


def forecasts(meta, rows, shape, annuli=None):
    radius = np.array([x["rad_kpc"] for x in rows])
    bar = bare_speed2(rows)
    mass = 1e9 * (0.5 * meta["L36_gsun"] + galaxies.sparc.GAS_HELIUM * meta["MHI_gsun"])
    hops, share = shape
    fields = []
    for scale in RATIOS:
        shells = hops * meta["disk_scale_kpc"] * scale
        if annuli is None:
            inside = (radius[:, None] >= shells) @ share
        else:
            edges, seed = annuli
            frac = annulus_fraction(
                radius[:, None, None], shells[None, :, None], edges[:-1], edges[1:]
            )
            inside = (frac @ seed) @ share
        fields.append(inside)
    enclosed = np.array(fields)
    assert np.all(np.diff(enclosed, axis=1) >= -1e-10)
    return np.sqrt(bar + galaxies.G * mass * CAP * enclosed / radius), enclosed


def choose(errors, train):
    if not train:
        raise ValueError("No calibration galaxies")
    return int(np.argmin(np.mean([errors[k] for k in sorted(train)], axis=0)))


def compare(sample, meta, shapes, annuli=None):
    truth = {k: np.array([r["vobs_kms"] for r in rows]) for k, rows in sample.items()}
    inputs = {
        k: dict(
            radius_kpc=[r["rad_kpc"] for r in rows],
            observed_kms=truth[k].tolist(),
            error_kms=[r["e_vobs_kms"] for r in rows],
            baryons_kms=np.sqrt(bare_speed2(rows)).tolist(),
        )
        for k, rows in sample.items()
    }
    predictions = {
        "baryons": {k: np.array(v["baryons_kms"]) for k, v in inputs.items()}
    }
    results = {}
    for route, shape in shapes.items():
        for seeded in [False] if annuli is None else [False, True]:
            name = route + ("_mapped" if seeded else "_central")
            banks = {
                k: forecasts(meta[k], rows, shape, annuli[k] if seeded else None)
                for k, rows in sample.items()
            }
            loss = {
                k: np.mean((curves - truth[k]) ** 2, axis=1)
                for k, (curves, _) in banks.items()
            }
            held = {}
            predictions[name] = {}
            for k, (curves, profiles) in banks.items():
                train = {
                    j
                    for j in sample
                    if (
                        j != k
                        if annuli is not None
                        else galaxies.fold(j) != galaxies.fold(k)
                    )
                }
                ix = choose(loss, train)
                held[k] = dict(
                    hop_per_disk_scale=RATIOS[ix],
                    speed_kms=curves[ix].tolist(),
                    enclosed_fraction=profiles[ix].tolist(),
                )
                predictions[name][k] = curves[ix]
            results[name] = held
    metrics = {}
    for name, rows in predictions.items():
        residuals = {k: v - truth[k] for k, v in rows.items()}
        metrics[name] = dict(
            pooled_median_abs_kms=np.median(
                abs(np.concatenate(list(residuals.values())))
            ),
            galaxy_equal_rmse_kms=np.sqrt(
                np.mean([np.mean(r * r) for r in residuals.values()])
            ),
            within_two_speed_errors=sum(
                int(sum(abs(r) <= 2 * np.array(inputs[k]["error_kms"])))
                for k, r in residuals.items()
            ),
            whole_curves_within_two=sum(
                bool(np.all(abs(r) <= 2 * np.array(inputs[k]["error_kms"])))
                for k, r in residuals.items()
            ),
        )
    return dict(
        galaxies=len(sample),
        rings=sum(map(len, sample.values())),
        inputs=inputs,
        forecasts=results,
        scores=metrics,
    )


def run():
    path = source_reader.HERE / "source_reader.json"
    shapes = observed_shapes(path)
    meta, groups = galaxies.catalogue()
    kept, _ = galaxies.select(meta, groups)
    charts = mapped_annuli(maps.enroll()[0])
    mapped = {
        k: sorted(
            [
                r
                for r in groups[k]
                if r["rad_kpc"] > 0
                and r["vobs_kms"] > 0
                and 0 < r["e_vobs_kms"] <= 0.1 * r["vobs_kms"]
                and bare_speed2([r])[0] > 0
            ],
            key=lambda r: r["rad_kpc"],
        )
        for k in charts
    }
    assert all(mapped.values())
    files = [
        path,
        bao.HERE / "SOURCE_PROFILE_PROTOCOL.md",
        *[galaxies.sparc.DATA / f for f in ("table1.dat", "table2.dat")],
        *[maps.HERE / f for f in ("leroy_sample.dat", "leroy_table7.dat")],
    ]
    return dict(
        sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in files},
        granted_inventory_ratio=CAP,
        scale_grid=RATIOS.tolist(),
        bulgeless=compare(kept, meta, shapes),
        resolved_maps=compare(mapped, meta, shapes, charts),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.loads(json.dumps(run(), allow_nan=False))
    path = bao.HERE / "source_profile.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), out)
    for k in ("bulgeless", "resolved_maps"):
        print(k, out[k]["galaxies"], out[k]["rings"], out[k]["scores"])
