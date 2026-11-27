"""Conditional swap and local-term readouts for a unit-free mass comparison."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

from ..row_14_profile import measure as stock
from ..row_14_profile.recoil_radius import measure as impulse
from ..row_14_profile.scattering.measure import verify

HERE = Path(__file__).resolve().parent
STARTS = 2**44
TARGET = 9396 / 9383
SEED = 150141


def record(rng, levels, weights, acceptance):
    p = np.asarray(weights)
    if abs(p.sum() - 1) > 1e-10 or p.min() < -1e-12:
        raise ArithmeticError("Incomplete outcome partition")
    p = np.maximum(p, 0)
    n = rng.binomial(STARTS, acceptance)
    counts = rng.multinomial(n, p / p.sum())
    return dict(
        levels=list(map(float, levels)), counts=counts.tolist(), refused=STARTS - n
    )


def interval(bank):
    levels = np.array(bank["levels"])
    n = sum(bank["counts"])
    mean = bank["counts"] @ levels / n
    rad = np.ptp(levels) * np.sqrt(np.log(2 * 256 / 0.01) / (2 * n))
    return [mean - rad, mean + rad]


def zero_kick(rng, words, swaps, rho, acceptance, d, edge2):
    i = np.arange(len(words))
    diag = np.diag(rho)
    mobile = (words[:, stock.EDGES[:, 0]] == 0) ^ (words[:, stock.EDGES[:, 1]] == 0)
    off = rho[i, swaps].T
    p = np.sum(mobile * diag[:, None]) / 30
    coherence = np.sum(mobile * off) / 30
    scale = 30 * edge2 / (6 * d * d)
    return record(
        rng,
        [scale, -scale, 0],
        [(p + coherence) / 2, (p - coherence) / 2, 1 - p],
        acceptance,
    )


def energy(rng, rho, maps, contrast, reference, k, acceptance):
    weights = {}
    diag = np.diag(rho)
    i = np.arange(len(rho))
    # Each selected seam exposes its commuting mismatch and exchange parity.
    for j, perm in enumerate(maps):
        for parity in (0, 1):
            prob = (diag + (-1) ** parity * rho[i, perm]) / 2
            for v in (0, 1, 4):
                e = parity + k * v
                weights[e] = weights.get(e, 0) + np.sum(prob[contrast[:, j] == v]) / 39
    for v in (-1, 0, 3):
        e = k * v
        weights[e] = weights.get(e, 0) + np.sum(diag[:, None] * (reference == v)) / 39
    labels = sorted(weights)
    return record(rng, 39 * np.array(labels), [weights[e] for e in labels], acceptance)


def ratio(num, den):
    if min(*num, *den) <= 0:
        return None
    return [num[0] / den[1], num[1] / den[0]]


def scale_score(bounds):
    if bounds is None:
        return None
    nearest = min(max(TARGET, bounds[0]), bounds[1])
    return {
        "neutral_to_positive": bounds,
        "best_common_unit_relative_error_at_least": abs(nearest - TARGET)
        / (nearest + TARGET),
        "neutral_MeV_if_positive_anchored_to_938p3": [938.3 * x for x in bounds],
    }


def compare(stocks):
    en, ep, ev = (
        interval(stocks[name]["energy"]) for name in ("neutral", "positive", "vacancy")
    )
    cn, cp, cv = (
        interval(stocks[name]["curvature"])
        for name in ("neutral", "positive", "vacancy")
    )
    sums = [x + y for x, y in zip(ep, ev, strict=True)]
    equal = [(9 * x + y) / 16 for x, y in zip(cp, cv, strict=True)]
    if min(*cp, *cv) <= 0:
        optimum = None
    else:
        optimum = ratio([x * y / (x + y) for x, y in zip(cp, cv, strict=True)], cn)
    return {
        name: scale_score(bounds)
        for name, bounds in {
            "formation_data": ratio(en, ep),
            "formation_fixed_inventory": ratio(en, sums),
            "inertia_data": ratio(cp, cn),
            "inertia_fixed_inventory_equal": ratio(equal, cn),
            "inertia_fixed_inventory_minimum_work": optimum,
        }.items()
    }


def run():
    rng = np.random.default_rng(SEED)
    result = {
        "starts": STARTS,
        "seed": SEED,
        "proton_MeV": 938.3,
        "neutron_MeV": 939.6,
        "materials": {},
    }
    mask, _ = stock.native_geometry()
    v = impulse.vertices()
    a, b = stock.EDGES.T
    lengths = ((v[a] - v[b]) ** 2).sum(axis=1)
    assert np.max(np.abs(lengths - lengths[0])) < 1e-13
    result["edge_length_squared"] = float(lengths.mean())
    for tag, k in (("native", 8), ("half_contrast", 4)):
        sources = {}
        for name, (plus, minus) in stock.NAMES.items():
            words, maps = stock.sector(plus, minus)
            kin, co, ref = stock.operators(words, maps, mask)
            rho, spec = stock.filtered(kin, co, ref, k, True)
            acc = spec["acceptance"]
            sources[name] = dict(
                clock=spec["target"],
                acceptance=acc,
                curvature=zero_kick(
                    rng, words, maps, rho, acc, plus + minus, lengths.mean()
                ),
                energy=energy(rng, rho, maps, co, ref, k, acc),
            )
        sources["filled"] = {
            "energy": record(rng, [0], [1], 1),
            "curvature": record(rng, [0], [1], 1),
        }
        result["materials"][tag] = {"sources": sources, "comparison": compare(sources)}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    path = parser.parse_args().record
    result = run()
    if path:
        with path.open("xb") as stream:
            stream.write(
                gzip.compress(
                    json.dumps(result, separators=(",", ":")).encode(), mtime=0
                )
            )
    else:
        verify(
            result, json.loads(gzip.decompress((HERE / "counts.json.gz").read_bytes()))
        )
    for tag, panel in result["materials"].items():
        print(tag)
        for name, entry in panel["comparison"].items():
            print(name, entry)
