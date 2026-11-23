"""Use two native defect clocks to test the common-lapse interpretation."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

from .. import neutral_defect as nd
from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import (
    defect_generator,
    free_seam,
)
from model.candidate_laws.contrast_memory_spherical_v1 import law
from . import force, measure

SHOTS = 2**32
TIMES = np.arange(1, 13) / 20


def response(w):
    return (np.sin(w * TIMES / 2) / w) ** 2


def fit(counts):
    observation = np.array(counts) / SHOTS

    def deviance(w):
        p = np.clip(response(w), 1e-18, 1 - 1e-15)
        return -np.sum(observation * np.log(p) + (1 - observation) * np.log1p(-p))

    grid = np.arange(0.5, 50.501, 0.1)
    best = grid[np.argmin([deviance(w) for w in grid])]
    freq = minimize_scalar(
        deviance,
        bounds=(best - 0.1, best + 0.1),
        method="bounded",
        options=dict(xatol=1e-11),
    ).x
    p = response(freq)
    slope = 0.5 * TIMES * np.sin(freq * TIMES) / freq**2 - 2 * p / freq
    uncertainty = (SHOTS * np.sum(slope**2 / np.maximum(p * (1 - p), 1e-18))) ** (-0.5)
    return float(freq), float(uncertainty)


def run():
    board = federation(0)
    edges, paths = force.cell_paths(board)
    seam = free_seam(nd.KAPPA)

    def local(l, r):
        i = law.SRC.index((l, r))
        return seam[i, i]

    records = []
    for distance, path in sorted(paths.items()):
        _, gateway = edges[path[-2], path[-1]]
        partner = min(p for p in board.neighbours[gateway] if board.degree[p] == 5)
        for ci, sign in enumerate((1, -1)):
            base = defect_generator(board, sign, nd.KAPPA, 0)[[gateway, partner]][
                :, [gateway, partner]
            ].toarray()
            for oi, source in enumerate((0, -1, 1)):
                applied = base.copy()
                if distance == 1:
                    applied[0, 0] += (
                        local(1 + source, 1 + sign)
                        - local(1 + source, 1)
                        - local(1, 1 + sign)
                        + local(1, 1)
                    )
                frequencies = np.linalg.eigvalsh(applied)
                omega = frequencies[1] - frequencies[0]
                p = np.array([abs(expm(-1j * t * applied)[1, 0]) ** 2 for t in TIMES])
                reads = [
                    measure.read_ramsey(
                        prob, SHOTS, [measure.SEED, 7, distance, oi, ci, j]
                    )
                    for j, prob in enumerate(p)
                ]
                sampled, se = fit([r["click"] for r in reads])
                records.append(
                    dict(
                        distance=distance,
                        source_gateway_deviation=source,
                        clock_load=sign,
                        rails=[gateway, partner],
                        frequency=float(omega),
                        sampled_frequency=sampled,
                        standard_error=se,
                        counts=reads,
                    )
                )
    return dict(
        protocol_sha256=hashlib.sha256(
            (measure.HERE / "CLOCK_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        times=TIMES.tolist(),
        records=records,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    for r in result["records"][:6]:
        print(
            r["clock_load"],
            r["source_gateway_deviation"],
            r["frequency"],
            r["sampled_frequency"],
            r["standard_error"],
        )
