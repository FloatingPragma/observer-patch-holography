"""Read the current of both separately commissioned cells on one admission record."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from . import measure
from .. import measure as stock
from ..scattering.measure import verify

JOINT_STARTS = 2**48


def response(source, coordinate, u):
    """One-axis quadrature after the signed kick, grouped by absolute centroid shift."""
    swaps, rho, _ = source
    dx = coordinate[swaps] - coordinate
    span, label = np.unique(np.round(np.abs(dx), 14), return_inverse=True)
    label = label.ravel()
    populations = np.broadcast_to(np.diag(rho), swaps.shape)
    coherence = rho[np.arange(len(rho)), swaps]
    imag = coherence * np.sin(u * dx) * np.sign(dx)
    probabilities = [
        np.bincount(label, weights=((populations + sign * imag) / 60).ravel())
        for sign in (1, -1)
    ]
    return (
        np.column_stack((15 * span, -15 * span)).ravel(),
        np.column_stack(probabilities).ravel(),
        float((dx**4).sum(axis=0).max()),
    )


def joint_bank(rng, first, second, coordinates, u, f):
    x, y = coordinates
    a, pa, da = response(first, x, f * u)
    b, pb, db = response(second, y, (1 - f) * u)
    levels = np.concatenate((2 * f * a, 2 * (1 - f) * b))
    p = np.concatenate((pa, pb)) / 2
    ao, bo = first[2], second[2]
    adm = rng.multinomial(
        JOINT_STARTS, [ao * bo, ao * (1 - bo), (1 - ao) * bo, (1 - ao) * (1 - bo)]
    )
    counts = rng.multinomial(adm[0], np.maximum(p, 0) / np.maximum(p, 0).sum())
    return {
        "levels": levels.tolist(),
        "counts": counts.tolist(),
        "admissions": adm.tolist(),
        "remainder_const": f**4 * da + (1 - f) ** 4 * db,
    }


def decode(axes, u):
    data = []
    for axis in axes:
        plus, minus = axis[str(u)], axis[str(-u)]
        vp, ep = measure.read(plus, plus["levels"])
        vm, em = measure.read(minus, minus["levels"])
        data.append(
            [
                (vp - vm) / (2 * u),
                (ep + em) / (2 * u) + u * u * plus["remainder_const"] / 12,
            ]
        )
    return np.mean(data, axis=0).tolist()


def run():
    result = {"joint_starts": JOINT_STARTS, "seed": 1401625, "settings": {}}
    rng = np.random.default_rng(result["seed"])
    xyz = measure.vertices()
    mask, _ = stock.native_geometry()
    for tag, k, enable in stock.PANELS:
        prepared = {}
        coords = {}
        for name in ("neutral", "positive", "vacancy"):
            words, sw = stock.sector(*stock.NAMES[name])
            kin, co, ref = stock.operators(words, sw, mask)
            rho, spec = stock.filtered(kin, co, ref, k, enable)
            prepared[name] = (sw, rho, spec["acceptance"])
            coords[name] = (words != 0) @ xyz / (words[0] != 0).sum()
        filled = (np.zeros((30, 1), dtype=int), np.ones((1, 1)), 1.0)
        vacancy = prepared["vacancy"]
        cal = [
            measure.inertia(rng, vacancy[1], vacancy[0], x, vacancy[2])
            for x in coords["vacancy"].T
        ]
        info = {
            "vacancy": cal,
            "vacancy_recoil": {
                str(u): measure.recoil_decode(cal, u) for u in (0.25, 0.5)
            },
            "joint": {},
        }
        for label, first, second, q, frac in (
            ("neutral_filled", "neutral", filled, np.zeros((1, 3)), 1.0),
            (
                "positive_vacancy_equal_load",
                "positive",
                vacancy,
                coords["vacancy"],
                0.75,
            ),
            ("positive_vacancy_data_only", "positive", vacancy, coords["vacancy"], 1.0),
        ):
            axes = [
                {
                    str(u): joint_bank(
                        rng,
                        prepared[first],
                        second,
                        (x, y),
                        u,
                        frac,
                    )
                    for u in measure.KICKS
                }
                for x, y in zip(coords[first].T, q.T, strict=True)
            ]
            info["joint"][label] = {
                "axes": axes,
                "inverse_inertia": {str(u): decode(axes, u) for u in (0.25, 0.5)},
            }
        result["settings"][tag] = info
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    f = parser.parse_args().record
    out = run()
    if f:
        with f.open("xb") as stream:
            stream.write(
                gzip.compress(json.dumps(out, separators=(",", ":")).encode(), mtime=0)
            )
    else:
        verify(
            out,
            json.loads(
                gzip.decompress((measure.HERE / "storage.json.gz").read_bytes())
            ),
        )
    for tag, panel in out["settings"].items():
        print(tag, "vacancy", panel["vacancy_recoil"]["0.25"])
        for name, pair in panel["joint"].items():
            print(name, pair["inverse_inertia"])
