"""Predict held rotation curves from the gradient of measured trail incidence."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import source_profile as profiles
from . import well_certificate
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    cold_release as lab,
)

HERE = profiles.bao.HERE
WINDOW = 0.02


def shapes(receipt):
    board = lab.federation(1)
    origin = lab.distances(board, [0])
    result = {}
    for route in ("outward", "folded"):
        tape = receipt["rows"][route + "_1_8"]["tape"]
        mask, _ = lab.read_map(board, tape)
        counts = np.bincount(origin[np.asarray(mask)].ravel())
        hops = np.flatnonzero(counts)
        result[route] = hops, counts[hops] / counts.sum()
    return result


def density_times_length(r, lengths, shape, annuli):
    """s times the spherically averaged density at one radius, s in kpc."""
    lengths = np.asarray(lengths)
    edges, seed = annuli
    a, b = edges[:-1], edges[1:]
    weight = seed / (b**2 - a**2)
    hops, shares = shape
    result = np.zeros_like(lengths, dtype=float)
    for hop, share in zip(hops, shares, strict=True):
        if hop == 0:
            interior = (a < r) & (r < b)
            boundary = (a == r) | (b == r)
            result += (
                lengths
                * (share / (2 * np.pi * r))
                * (weight @ (interior + 0.5 * boundary))
            )
            continue
        if len(lengths) <= 128:
            h = (lengths * hop)[:, None]
            overlap = np.maximum(0, np.minimum(b, r + h) - np.maximum(a, abs(r - h)))
            result += share * (overlap @ weight) / (4 * np.pi * r * hop)
            continue
        # Each annulus overlap is affine between these breakpoints. Evaluate
        # there, then interpolate; this also permits the continuous shape test.
        for lo, hi, w in zip(a, b, weight, strict=True):
            if w == 0:
                continue
            knots = np.unique(
                np.array([0, abs(r - lo), abs(r - hi), r, r + lo, r + hi])
            )
            overlap = np.maximum(
                0, np.minimum(hi, r + knots) - np.maximum(lo, abs(r - knots))
            )
            overlap[0] = overlap[-1] = 0
            result += (w * share / (4 * np.pi * r * hop)) * np.interp(
                lengths, knots / hop, overlap
            )
    return result


def field_times_length(r, lengths, shape, annuli):
    low = density_times_length((1 - WINDOW) * r, lengths, shape, annuli)
    high = density_times_length((1 + WINDOW) * r, lengths, shape, annuli)
    return (low - high) / np.log((1 + WINDOW) / (1 - WINDOW))


def shape_test(rows, annuli, shape):
    r = np.asarray(rows["radius_kpc"])
    obs, bar, err = [
        np.asarray(rows[k]) for k in ("observed_kms", "baryons_kms", "error_kms")
    ]
    lower = np.maximum(0, obs - 2 * err) ** 2 - bar**2 - 0.001
    upper = (obs + 2 * err) ** 2 - bar**2 + 0.001
    ix = np.flatnonzero(lower > 0)
    if not len(ix):
        return dict(ring_indices=[], excluded=False, intervals=0)
    ix = ix[np.rint(np.linspace(0, len(ix) - 1, min(7, len(ix)))).astype(int)]
    hops = shape[0][shape[0] > 0]
    R = r[ix, None] * np.array([1 - WINDOW, 1 + WINDOW])
    rr = R.ravel()[:, None]
    e = annuli[0]
    offsets = np.unique(np.r_[abs(rr - e).ravel(), (rr + e).ravel(), rr.ravel()])
    knots = np.unique((offsets[:, None] / hops).ravel())
    knots = np.r_[0, knots[knots > 0]]
    F = np.array([field_times_length(x, knots, shape, annuli) for x in r[ix]])
    # On the two open tails the shape is constant up to normalization.
    # On each interior interval positive rescalings sweep the cone of its
    # endpoint vectors, so normalizing each endpoint preserves feasibility.
    F[:, 0] = F[:, 1]
    norms = np.max(abs(F), axis=0)
    F = F / np.where(norms > 0, norms, 1)
    a, b = F[:, :-1], F[:, 1:]
    start = np.zeros(a.shape[1])
    end = np.ones_like(start)
    possible = (np.max(abs(a), axis=0) + np.max(abs(b), axis=0)) > 0

    def impose(left, right):
        # Round in favor of feasibility, after endpoint normalization.
        tol = 2e-10
        left, right = left + tol, right + tol
        possible[:] &= ~((left < 0) & (right < 0))
        x = (left < 0) & (right >= 0)
        y = (left >= 0) & (right < 0)
        start[x] = np.maximum(start[x], left[x] / (left[x] - right[x]))
        end[y] = np.minimum(end[y], left[y] / (left[y] - right[y]))

    for i, row in enumerate(ix):
        impose(a[i], b[i])
        for j, other in enumerate(ix):
            # Eliminate the common positive amplitude from lower_i/F_i <= upper_j/F_j.
            ratio = lower[row] / upper[other]
            impose(a[i] - ratio * a[j], b[i] - ratio * b[j])
    possible &= start <= end + 2e-10
    # If a point with an all-zero profile survives the relaxation, this
    # returns inconclusive. It cannot turn a feasible curve into an exclusion.
    return dict(
        ring_indices=ix.tolist(),
        excluded=not bool(np.any(possible)),
        intervals=len(start),
    )


def calibrate(fields, targets, names):
    cross = np.mean([np.mean(fields[n] * targets[n], axis=1) for n in names], axis=0)
    power = np.mean([np.mean(fields[n] ** 2, axis=1) for n in names], axis=0)
    gain = np.maximum(
        0, np.divide(cross, power, out=np.zeros_like(power), where=power > 0)
    )
    loss = np.mean(
        [np.mean((gain[:, None] * fields[n] - targets[n]) ** 2, axis=1) for n in names],
        axis=0,
    )
    index = int(np.argmin(loss))
    return index, float(gain[index])


def evaluate(inputs, charts, meta, shape):
    fields, targets = {}, {}
    for name, rows in inputs.items():
        scale = profiles.RATIOS * meta[name]["disk_scale_kpc"]
        mass = (
            0.5 * meta[name]["L36_gsun"]
            + profiles.galaxies.sparc.GAS_HELIUM * meta[name]["MHI_gsun"]
        )
        fields[name] = (
            mass
            * np.array(
                [
                    field_times_length(r, scale, shape, charts[name]) / scale
                    for r in rows["radius_kpc"]
                ]
            ).T
        )
        targets[name] = (
            np.asarray(rows["observed_kms"]) ** 2 - np.asarray(rows["baryons_kms"]) ** 2
        )
    out, errors, within, no_orbits = {}, [], 0, 0
    for name, rows in inputs.items():
        index, gain = calibrate(fields, targets, sorted(set(inputs) - {name}))
        v2 = np.asarray(rows["baryons_kms"]) ** 2 + gain * fields[name][index]
        speed = np.sqrt(np.maximum(0, v2))
        residual = abs(speed - rows["observed_kms"])
        inside = (v2 >= 0) & (residual <= 2 * np.asarray(rows["error_kms"]))
        out[name] = dict(
            hop_per_disk_scale=profiles.RATIOS[index],
            potential_gain=gain,
            squared_speed=v2.tolist(),
            within_two=int(inside.sum()),
            no_orbits=int(sum(v2 < 0)),
            continuous_shape=shape_test(rows, charts[name], shape),
        )
        errors.extend(residual)
        within += int(inside.sum())
        no_orbits += int(sum(v2 < 0))
    return dict(
        galaxies=out,
        rings=len(errors),
        median_abs_kms=float(np.median(errors)),
        within_two=within,
        no_orbits=no_orbits,
        excluded_shapes=sum(x["continuous_shape"]["excluded"] for x in out.values()),
    )


def run():
    source = lab.HERE / "cold_release.json"
    prior = HERE / "source_profile.json"
    models = shapes(json.loads(source.read_bytes()))
    inputs = json.loads(prior.read_bytes())["resolved_maps"]["inputs"]
    charts = profiles.mapped_annuli(profiles.maps.enroll()[0])
    meta, _ = profiles.galaxies.catalogue()
    assert set(charts) == set(inputs)
    files = [
        source,
        prior,
        HERE / "LOCAL_WELL_PROTOCOL.md",
        profiles.maps.HERE / "leroy_table7.dat",
        profiles.maps.HERE / "leroy_sample.dat",
        profiles.galaxies.sparc.DATA / "table1.dat",
    ]
    return dict(
        sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in files},
        ddo154_exact={
            route: well_certificate.certify(
                inputs["DDO154"], charts["DDO154"], models[route], indices
            )
            for route, indices in (("outward", [0, 2]), ("folded", [0, 2, 5]))
        },
        window_fraction=WINDOW,
        mass_unit_solar=10**9,
        scale_grid=profiles.RATIOS.tolist(),
        shapes={
            n: dict(hops=h.tolist(), share=f.tolist()) for n, (h, f) in models.items()
        },
        routes={
            n: evaluate(inputs, charts, meta, shape) for n, shape in models.items()
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    record = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "local_well.json"
    if args.write:
        path.write_text(json.dumps(record, indent=2) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), record)
    for name, row in record["routes"].items():
        print(name, {k: v for k, v in row.items() if k != "galaxies"})
