"""Paired and separated responses to a cell-wide signed native load pulse."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import defect_generator
from .. import neutral_defect as defect
from . import measure


TIMES = (0.25, 1.0, 4.0)
SHOTS = 2**38


def partition(weights, bins, n):
    return np.bincount(bins, weights=weights, minlength=n)


def counts(prob, seed):
    return np.random.default_rng(seed).multinomial(SHOTS, prob / prob.sum()).tolist()


def contrast(readings, values):
    probabilities = np.array([s["probabilities"] for s in readings])
    observed = np.array([s["counts"] for s in readings]) / SHOTS
    variance = np.sum(probabilities @ values**2 - (probabilities @ values) ** 2) / SHOTS
    return dict(
        exact=float((probabilities[0] - probabilities[1]) @ values),
        sampled=float((observed[0] - observed[1]) @ values),
        standard_error=math.sqrt(variance),
    )


def singles(board, edge, sign, index, source_idx):
    h = defect_generator(board, sign, defect.KAPPA, 0)
    right = (np.arange(board.n) // 12 == edge[1] // 12).astype(int)
    source = np.zeros(board.n)
    source[list(edge)] = 1 / math.sqrt(2)
    initial = np.column_stack(
        [np.exp(p * 1j * math.pi / 2 * sign * right) * source for p in (1, -1)]
    )
    records = []
    for j, t in enumerate(TIMES):
        final = expm_multiply(
            (-1j * t) * h, initial, traceA=-1j * t * h.diagonal().sum()
        )
        readings = []
        for p, col in enumerate(final.T):
            prob = partition(abs(col) ** 2, right, 2)
            readings.append(
                dict(
                    probabilities=prob.tolist(),
                    counts=counts(prob, [measure.SEED, index, 5, source_idx, j, p]),
                )
            )
        records.append(
            dict(
                time=t,
                readings=readings,
                cell_deflection=contrast(readings, np.array([0.0, 1.0])),
            )
        )
    return records


def pairs(board, edge, index, h, x, y, high):
    _, sharp, filtered, _, _ = measure.source_columns(board, edge, h)
    right = edge[1] // 12
    left = edge[0] // 12
    charge = (x // 12 == right).astype(int) - (y // 12 == right).astype(int)
    bins = np.full(len(x), 2, dtype=int)
    bins[high & (x // 12 == left) & (y // 12 == left)] = 0
    bins[high & (x // 12 == right) & (y // 12 == right)] = 1
    keys = 3 * bins + charge + 1
    # Compress to real inputs, evolving the real and imaginary parts separately.
    sources = []
    for states in (sharp, filtered):
        source = states.sum(axis=1) / math.sqrt(2)
        rotation = np.exp(1j * math.pi / 2 * charge) * source
        sources.extend([rotation.real, rotation.imag])
    flights = measure.propagation(h, np.column_stack(sources), TIMES)
    sources_result = {}
    for si, name in enumerate(("sharp", "filtered")):
        data = []
        for j, t in enumerate(TIMES):
            real, imag = flights[t][:, 2 * si : 2 * si + 2].T
            readings = []
            for p, phase in enumerate((1, -1)):
                final = real + phase * 1j * imag
                prob = partition(abs(final) ** 2, keys, 9)
                readings.append(
                    dict(
                        probabilities=prob.tolist(),
                        counts=counts(prob, [measure.SEED, index, 5, si + 2, j, p]),
                    )
                )
            data.append(
                dict(
                    time=t,
                    readings=readings,
                    signed_number_change=contrast(
                        readings, np.tile([-1.0, 0.0, 1.0], 3)
                    ),
                    intact_imbalance=contrast(readings, np.repeat([-1.0, 1.0, 0.0], 3)),
                )
            )
        sources_result[name] = data
    return sources_result


def run():
    board = federation(0)
    x, y, _, _, h, cost = defect.generators(board)
    output = []
    for index in measure.PLACEMENTS:
        bridge = board.external[index]
        output.append(
            dict(
                bridge=list(bridge),
                excess=singles(board, bridge, 1, index, 0),
                vacancy=singles(board, bridge, -1, index, 1),
                **pairs(board, bridge, index, h, x, y, cost == 4 * defect.KAPPA),
            )
        )
    return dict(
        protocol_sha256=hashlib.sha256(
            (measure.HERE / "PULSE_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        placements=output,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as file:
        json.dump(result, file, indent=2, allow_nan=False)
        file.write("\n")
    for survey in result["placements"]:
        print(
            survey["bridge"],
            [
                (s, survey[s][1]["cell_deflection"]["sampled"])
                for s in ("excess", "vacancy")
            ],
            survey["filtered"][1]["signed_number_change"],
            survey["filtered"][1]["intact_imbalance"],
        )
