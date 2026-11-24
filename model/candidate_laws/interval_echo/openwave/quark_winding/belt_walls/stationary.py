"""Optimize the three-gap records over stationary states on each shared comparison ring."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate, currents
from . import belt, composite, multiplets


def ring_blocks(length, n):
    words = [composite.mask(s) for s in combinations(range(length), n)]
    all_bits = (1 << length) - 1
    orbital = {}
    orbits = []
    for w in words:
        if w in orbital:
            continue
        orbit = []
        while w not in orbital:
            orbital[w] = (len(orbits), len(orbit))
            orbit.append(w)
            w = ((w << 1) & all_bits) | (w >> (length - 1))
        orbits.append(orbit)
    return orbital, orbits


def scan(length, ones):
    n = length // 3 + ones - 1
    orbital, orbits = ring_blocks(length, n)
    moves = composite.hop_list(belt.ring(length).seams)

    def sites(w):
        return tuple(i for i in range(length) if w & (1 << i))

    arcs = np.array(
        [composite.gap_arc(sites(o[0]), length, ones, 3 - ones) for o in orbits],
        dtype=float,
    )
    low = np.nanmin(arcs)
    masks = {
        "three_walls": np.isfinite(arcs),
        "minimum_arc": arcs == low,
        "minimum_plus_three": arcs <= low + 3,
    }
    best = {k: {"probability": -1} for k in masks}
    dimensions = []
    for k in range(length // 2 + 1):  # conjugation relates momentum k and -k
        subset = [i for i, o in enumerate(orbits) if k * len(o) % length == 0]
        index = {o: j for j, o in enumerate(subset)}
        h = np.zeros((len(subset), len(subset)), complex)
        for j, o in enumerate(subset):
            orbit = orbits[o]
            h[j, j] = composite.diagonal(orbit[0], moves, candidate.KAPPA)
            for other in composite.neighbours(orbit[0], moves):
                b, t = orbital[other]
                if b in index:
                    h[index[b], j] -= (
                        0.5
                        * np.sqrt(len(orbit) / len(orbits[b]))
                        * np.exp(-2j * np.pi * k * t / length)
                    )
        np.testing.assert_allclose(h, h.conj().T, atol=2e-13)
        e, v = np.linalg.eigh(h)
        dimensions.append(len(e))
        blocks = multiplets.groups(e)
        for key, mask in masks.items():
            chosen = mask[subset]
            restrict = v[chosen]
            weights = (abs(restrict) ** 2).sum(axis=0)
            ix = np.array([int(np.argmax(weights))])
            weight = weights[ix[0]]
            coeff = np.ones(1)
            for block in blocks:
                if len(block) == 1:
                    continue
                sub = restrict[:, block]
                eig, states = np.linalg.eigh(sub.conj().T @ sub)
                if eig[-1] > weight:
                    ix, weight, coeff = block, eig[-1], states[:, -1]
            if weight > best[key]["probability"]:
                probabilities = abs(v[:, ix] @ coeff) ** 2
                best[key] = dict(
                    probability=weight,
                    k=k,
                    energy=np.mean(e[ix]),
                    block_dimension=len(ix),
                    arc_weights={
                        str(a): probabilities[arcs[subset] == a].sum()
                        for a in np.unique(arcs[np.isfinite(arcs)])
                    },
                )
    return dict(
        length=length,
        loads=n,
        minimum_arc=low,
        momentum_dimensions=dimensions,
        maxima=best,
    )


def run():
    return currents._portable(
        [scan(length, ones) for length in (12, 15, 18) for ones in (1, 2)]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(encoded)
