"""Gaussian field completion using the positive frequencies of the routed massless scalar."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from .clocked_field import eigenphase

# Euclidean time followed by the three physical spatial coordinates.
POINTS = np.array([[0, 0, 0, 0], [2, 1, 0, 0], [4, 0, 1, 0], [6, 0, 0, 1]])
PAIRS = tuple(combinations(range(4), 2))
SPEED = 0.5


def wick(pairs):
    return pairs[0] * pairs[5] + pairs[1] * pairs[4] + pairs[2] * pairs[3]


def flat_pairs(points):
    separations = np.array([points[i] - points[j] for i, j in PAIRS])
    squared = (SPEED * separations[:, 0]) ** 2 + np.sum(separations[:, 1:] ** 2, axis=1)
    return 1 / (4 * np.pi**2 * SPEED * squared)


def periodic_pairs(points, length, sites):
    if type(sites) is not int or sites < 2 or sites % 2 or length <= 0:
        raise ValueError("Expected an even periodic lattice and positive box length")
    points = np.asarray(points, dtype=float)
    h = length / sites
    if points.shape != (4, 4) or not np.isfinite(points).all():
        raise ValueError("Four finite spacetime points are required")
    if not np.allclose(
        points[:, 1:] / h, np.rint(points[:, 1:] / h), rtol=0, atol=1e-12
    ):
        raise ValueError("Spatial probes must lie on lattice sites")
    deltas = np.array([points[i] - points[j] for i, j in PAIRS])
    if np.any(deltas[:, 0] == 0):
        raise ValueError("The comparison uses positive Euclidean time separations")
    k = 2 * np.pi * np.arange(sites // 2 + 1) / length
    multiplicity = np.full(len(k), 2)
    multiplicity[[0, -1]] = 1
    weights_yz = multiplicity[:, None] * multiplicity[None, :]
    angles = np.zeros((len(k), len(k), 4))
    angles[:, :, 1] = h * k[:, None]
    angles[:, :, 2] = h * k[None, :]
    cos_yz = [
        np.cos(k[:, None] * dy) * np.cos(k[None, :] * dz) for _, _, dy, dz in deltas
    ]
    result = np.zeros(len(PAIRS))
    for i, kx in enumerate(k):
        angles[:, :, 0] = h * kx
        omega = eigenphase(angles) / h
        if i == 0:
            omega[0, 0] = (
                np.inf
            )  # Fixed spatial average in the finite-volume vacuum prescription.
        measure = multiplicity[i] * weights_yz / (2 * omega * length**3)
        for n, (dt, dx, _, _) in enumerate(deltas):
            result[n] += np.cos(kx * dx) * np.sum(
                measure * np.exp(-abs(dt) * omega) * cos_yz[n]
            )
    return result


def run():
    target = flat_pairs(POINTS)
    rows = []
    for length in (24, 48):
        for mesh in (0.25, 0.125):
            values = periodic_pairs(POINTS, length, int(length / mesh))
            rows.append(
                dict(
                    box_length=length,
                    mesh=mesh,
                    pairs=values.tolist(),
                    wick_fourth=wick(values),
                    max_pair_relative_error=max(abs(values / target - 1)),
                    fourth_relative_error=abs(wick(values) / wick(target) - 1),
                )
            )
    return dict(
        quantization="canonical free oscillator vacuum, spatial zero mode omitted",
        points=POINTS.tolist(),
        pair_order=PAIRS,
        flat_pairs=target.tolist(),
        flat_wick_fourth=wick(target),
        regulators=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
