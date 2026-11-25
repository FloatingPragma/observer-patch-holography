"""Strict second-variation checks on the three-dimensional collar mixture family."""

import argparse
from itertools import combinations, product
from .verify import precision_check
import json
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import brentq
from .. import selected_station as prior
from ..geometry_gain import pooled_ranges


def geometry(r):
    # Same survey as simplex_area_volume, batched over reference points and derivatives.
    off = ~np.eye(4, dtype=bool)
    logs = np.zeros_like(r)
    logs[..., off] = np.log(r[..., off])
    gains = (logs - logs.swapaxes(-1, -2)).mean(axis=-1)
    corrected = logs - gains[..., None]
    sq = np.exp(corrected + corrected.swapaxes(-1, -2))
    ij = np.arange(4)
    sq[..., ij, ij] = 0

    def gram(inds):
        i, *other = inds
        sides = sq[..., i, other]
        g = (sides[..., None] + sides[..., None, :] - sq[..., other, :][..., other]) / 2
        return np.sqrt(np.linalg.det(g)) / (6 if len(other) == 3 else 2)

    return np.stack(
        [sum(gram(face) for face in combinations(range(4), 3)), gram(range(4))], axis=-1
    )


def affine(c):
    return np.array(
        [c[1, 0] - c[0, 0], c[0, 1] - c[0, 0], c[1, 1] + c[0, 0] - c[1, 0] - c[0, 1]]
    )


def gradients(corners, x):
    derivatives = affine(corners)
    base = corners[0, 0] + np.einsum("...d,dij->...ij", x, derivatives)
    return (
        geometry(base[..., None, :, :] + 1e-18j * derivatives).imag.swapaxes(-1, -2)
        / 1e-18
    )


def coords(first, second, s):
    return np.stack([first * (1 - s) + second * s, s + 0 * first, second * s], axis=-1)


def details(corners, click, x):
    a, v = gradients(corners, x)
    z = np.array([affine(c) for c in np.asarray(click).reshape(-1, 2, 2)])
    mean = z @ x
    first = (-np.arctanh(mean)) @ z / len(z)
    # Baseline directions du and ds at fixed measured volume, adjusted via dv.
    tangent = np.array([[1, 0, -v[0] / v[2]], [0, 1, -v[1] / v[2]]])
    da = tangent @ a
    ds = tangent @ first
    eta = -ds[0] / da[0]
    lam = (v @ (eta * a + first)) / (v @ v)
    h = 2e-4
    plus = gradients(corners, x + h * np.eye(3))
    minus = gradients(corners, x - h * np.eye(3))
    ha, hv = np.moveaxis((plus - minus) / (2 * h), 1, 0)
    hess = eta * ha - lam * hv - (z.T / (1 - mean**2)) @ z / len(z)
    basis = null_space(v[None, :])
    eig = np.linalg.eigvalsh(basis.T @ hess @ basis)
    return dict(
        x=x.tolist(),
        eta=float(eta),
        lagrange=float(lam),
        eigenvalues=eig.tolist(),
        residual=float(np.linalg.norm(tangent @ (eta * a + first))),
        area_slopes=da.tolist(),
        entropy_slopes=ds.tolist(),
    )


def click_bank(charge):
    board = prior.pairing.formation.federation(1)
    route = prior.pairing.formation.region(board)
    p, d = prior.head_source()
    rows = {}
    for odd in (False, True):
        h = prior.ranging.generator(board, route, charge, 8, odd)
        launch = np.concatenate(
            [prior.pairing.origin_columns(h, port) for port in route], axis=1
        )
        waves = prior.pairing.formation.flow(h, launch, prior.pairing.TIMES)
        for leg, wave in enumerate(waves):
            for left, right in product(range(4), repeat=2):
                src, dst = route[left], route[right]
                out = wave[4 * dst : 4 * dst + 4, 4 * left : 4 * left + 4]
                back = wave[4 * src : 4 * src + 4, 4 * right : 4 * right + 4]
                trip = back @ out
                effect = trip.conj().T @ trip
                key = (charge, leg, src, dst)
                entries = rows.setdefault(key, np.empty((2, 2)))
                entries[int(odd)] = [np.trace(effect @ rho).real for rho in (p, p + d)]
    return rows


def radar_corners(banks, index, horizon, radar):
    selected = {"mix": (-1, 1), "plus": (1,), "minus": (-1,)}[radar]
    return np.array(
        [
            [
                pooled_ranges([banks[q, odd][index] for q in selected], name, horizon)
                for name in ("independent", "formed")
            ]
            for odd in (False, True)
        ]
    )


def station_normal(x, z):
    if len(z) == 1:
        # Exclude entropy extrema with zero ordinary gain.
        return z[0]
    return -np.arctanh(x @ z.T) @ z / len(z)


def search(banks, probes):
    central = [k for k in probes if k[2:] == (1, 38)]
    menus = [(k,) for k in probes] + list(combinations(central, 2))
    rates = list(product((0.03125, 0.125, 0.25, 0.375, 0.46875), repeat=2))
    strength = np.linspace(0.001, 0.999, 11)
    grid = np.array([coords(a, b, strength) for a, b in rates])
    rows = []
    for index, horizon, radar in product((0, 1), (2, 4, 6), ("mix", "plus", "minus")):
        corners = radar_corners(banks, index, horizon, radar)
        gradients_at_grid = gradients(corners, grid)
        normals = np.cross(gradients_at_grid[..., 0, :], gradients_at_grid[..., 1, :])
        for menu in menus:
            records = np.array([probes[k] for k in menu])
            z = np.array([affine(c) for c in records])
            det = (normals * station_normal(grid, z)).sum(axis=-1)
            for rate_index, interval in np.argwhere(det[:, :-1] * det[:, 1:] < 0):
                first, second = rates[rate_index]

                def mismatch(s):
                    x = coords(first, second, s)
                    ga, gv = gradients(corners, x)
                    return np.cross(ga, gv) @ station_normal(x, z)

                lo, hi = strength[interval : interval + 2]
                star = brentq(mismatch, lo, hi, xtol=1e-13)
                row = details(corners, records, coords(first, second, star))
                row.update(
                    index=index,
                    horizon=horizon,
                    radar=radar,
                    menu=list(menu),
                    conditional_odd=[first, second],
                    s=star,
                    bracket=[lo, hi],
                    bracket_determinants=[mismatch(lo), mismatch(hi)],
                )
                rows.append(row)
    return rows


def run():
    banks = prior.radar_banks()
    probes = click_bank(-1) | click_bank(1)
    roots = search(banks, probes)
    summaries = []
    for length in (1, 2):
        subset = [r for r in roots if len(r["menu"]) == length]
        positive = [r for r in subset if r["eta"] > 0]
        maxima = [r for r in positive if r["eigenvalues"][1] < 0]
        summaries.append(
            dict(
                settings_per_record=length,
                bracketed=len(subset),
                positive_eta=len(positive),
                strict_negative_hessian=len(maxima),
                minimum_larger_eigenvalue_per_eta=min(
                    r["eigenvalues"][1] / r["eta"] for r in positive
                ),
            )
        )
    best = min(
        (r for r in roots if len(r["menu"]) == 2 and r["eta"] > 0),
        key=lambda r: r["eigenvalues"][1],
    )
    records = np.array([probes[k] for k in best["menu"]])
    corners = radar_corners(banks, best["index"], best["horizon"], best["radar"])
    eta = best["eta"]
    x = np.array(best["x"])
    out = []
    for name, geometry_source, detectors, where in (
        ("erased_probe_sign", corners, records.mean(axis=0), x),
        ("half_pairing", corners, records, coords(1 / 32, 1 / 32, best["s"] / 2)),
        ("doubled_pairing", corners, records, coords(1 / 32, 1 / 32, 2 * best["s"])),
        ("quarter_tick", radar_corners(banks, 0, 6, "plus"), records, x),
        ("horizon_four", radar_corners(banks, 1, 4, "plus"), records, x),
    ):
        r = details(geometry_source, detectors, where)
        ratios = -eta * np.array(r["area_slopes"]) / r["entropy_slopes"]
        out.append(
            dict(change=name, ordinary_transfer=ratios[0], paired_transfer=ratios[1])
        )
    return dict(
        search=summaries,
        selected=best,
        checked=precision_check(corners, records, x, eta),
        held_at_frozen_gain=out,
        raw_ranges=corners.tolist(),
        bridge_corners=records.tolist(),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    body = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_name("results.json").write_text(body)
    else:
        print(body, end="")
