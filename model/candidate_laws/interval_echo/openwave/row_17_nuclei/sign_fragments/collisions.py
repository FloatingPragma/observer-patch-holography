"""Compare sharp words at matched addresses on two independently cold fragment cells."""

import argparse
import json
from functools import cache
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.stats import beta

from model.candidate_laws.interval_echo.openwave.row_14_profile.measure import (
    native_geometry,
    operators,
    sector,
)
from .measure import HERE, SHOTS, STOCKS

SEED = 170127
STRENGTHS = (8, 4)


@cache
def fragment(plus, minus, k):
    words, swaps = sector(plus, minus)
    kinetic, contrast, private = operators(words, swaps, native_geometry()[0])
    h = kinetic.copy()
    h.flat[:: len(h) + 1] += k * (contrast.sum(axis=1) + private.sum(axis=1))
    spectrum, vectors = eigh(h, driver="evd")
    # Bound the numerical errors using the full orthonormal decomposition.
    residual = np.linalg.norm(h @ vectors - vectors * spectrum)
    defect = np.linalg.norm(vectors.T @ vectors - np.eye(len(h)))
    if residual > 1e-8 or defect > 1e-9:
        raise ValueError("fragment eigensystem outside checked accuracy")
    return words != 0, abs(vectors[:, 0]) ** 2, spectrum[:2]


def acquire():
    rng = np.random.default_rng(SEED)
    calibrations = []
    for overlap in (False, True):
        first = [1] + [0] * 11
        second = [0] * 12
        second[0 if overlap else 1] = 1
        measured = int(any(a and b for a, b in zip(first, second, strict=True)))
        calibrations.append([SHOTS * (1 - measured), SHOTS * measured])
    banks = []
    for plus, minus in STOCKS:
        for k in STRENGTHS:
            excess, a, ea = fragment(plus, 0, k)
            vacancy, b, eb = fragment(0, minus, k)
            common = excess.astype(int) @ vacancy.T > 0
            collision = float(a @ common @ b)
            counts = rng.multinomial(SHOTS, [1 - collision, collision]).tolist()
            banks.append(
                dict(
                    stock=[plus, minus],
                    k=k,
                    excess_levels=ea.tolist(),
                    vacancy_levels=eb.tolist(),
                    exact_model_probability=collision,
                    counts=counts,
                )
            )
    return dict(seed=SEED, shots=SHOTS, calibrations=calibrations, banks=banks)


def decode(data):
    settings = {(tuple(r["stock"]), r["k"]) for r in data["banks"]}
    if (
        data["shots"] != SHOTS
        or len(data["banks"]) != 8
        or settings != {(s, k) for s in STOCKS for k in STRENGTHS}
        or data["calibrations"] != [[SHOTS, 0], [0, SHOTS]]
    ):
        raise ValueError("incomplete fragment records")
    rows = []
    for r in data["banks"]:
        counts = r["counts"]
        if len(counts) != 2 or min(counts) < 0 or sum(counts) != SHOTS:
            raise ValueError("missing collision trials")
        c = counts[1]
        lower = float(beta.ppf(0.01 / 8, c, SHOTS - c + 1)) if c else 0.0
        gap = min(
            levels[1] - levels[0]
            for levels in (r["excess_levels"], r["vacancy_levels"])
        )
        rows.append(
            dict(
                stock=r["stock"],
                k=r["k"],
                collisions=c,
                gap=gap,
                probability_lower=lower,
                penalty_at_lower=gap * lower,
            )
        )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    if args.acquire:
        args.acquire.write_text(
            json.dumps(acquire(), separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        print(
            json.dumps(
                decode(json.loads((HERE / "collisions.json").read_text())), indent=2
            )
        )
