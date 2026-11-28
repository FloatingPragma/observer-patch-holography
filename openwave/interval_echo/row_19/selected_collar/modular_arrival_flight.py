"""Test source-only modular predictions on the restored 240-port receiver."""

from math import log, sqrt

from mpmath import iv
import numpy as np

from . import full_board_polarization as full
from . import reference_arrival as reference

local = full.local
N = 2**76
BATCH = 2**60
SEED = 28122922


def counts(weights, rng):
    p = np.maximum(0, weights)
    assert abs(p.sum() - 1) < 1e-12
    chunks = rng.multinomial(BATCH, p / p.sum(), size=N // BATCH)
    return [sum(map(int, col)) for col in chunks.T]


def band(k, scores):
    val = sum(a * b for a, b in zip(k, scores, strict=True)) / N
    half = (max(scores) - min(scores)) * (sqrt(log(2 * 10000 / 0.01) / (2 * N)) + 2e-15)
    return iv.mpf([val - half, val + half])


def run():
    iv.dps = 45
    rng = np.random.default_rng(SEED)
    initial = local.inputs("full_board_polarization")
    charge = local.inputs("reference_arrival")
    board = local.survey.bulk.federation(0)
    setup, _ = full.steady.symmetric_instrument(board)
    _, cells, ports, ab, _, degrees, edges, _ = setup
    gates = [ports[i] for i in ab]
    active, layouts, _, cert = full.native_layouts(board, setup)
    assert cert == initial["certificates"]
    sign = np.array(
        [
            int(j not in gates and tuple(sorted((gates[0], j))) in active)
            - int(j not in gates and tuple(sorted((gates[1], j))) in active)
            for j in range(board.n)
        ]
    )
    dt = full.motion.FLIGHTS[-1]
    panel = np.tile(full.empty_table(setup)[-1][None], (board.n, 1, 1))
    panel[ports] = full.motion.second_instruments(ports, ab, degrees, edges)[-1]
    pos = np.array([0, 1, 2, 0])
    quad = np.array([1, -1, 0, 0])
    responses = [(panel[:, k] - panel[:, k + 1]) @ pos / dt**2 for k in (0, 2)]
    responses.append((panel[:, 5] - panel[:, 7]) @ quad / dt)
    forward, reverse, _ = full.steady.parity(board, cells)
    results = {}
    for direction, law in (("forward", forward), ("reversed", reverse)):
        accepted = sum(float(w) * layouts[f] for f, w in law.items())
        old = initial["outcomes"][direction]
        acquired = {}
        for key, arm in charge["arms"].items():
            p = arm["collision_acceptance"]
            z = (1 - p) * layouts[0, 0] + p * accepted
            z = z / z.sum()
            # Both occupied launch gateways keep their own blocked category.
            bins = [z[sign == -1].sum(), z[gates].sum(), z[sign == 1].sum()]
            census = counts([bins[0], bins[1], 1 - sum(bins), bins[2]], rng)
            samples = [counts(row, rng) for row in np.einsum("j,jmk->mk", z, panel)]
            read = [
                (band(samples[k], pos.tolist()) - band(samples[k + 1], pos.tolist()))
                / dt**2
                for k in (0, 2)
            ]
            read.append(
                (band(samples[5], quad.tolist()) - band(samples[7], quad.tolist())) / dt
            )
            acquired[key] = (band(census, [-1, 0, 0, 1]), read, census, samples)
        base = acquired["fresh"]
        for key, (meter, read, k, shots) in acquired.items():
            arm = charge["arms"][key]
            predicted_ratio = 144 * iv.mpf(arm["arrival_from_word_measurements"])
            # The old receiver calibrated all response slopes on other trials.
            predicted = [
                predicted_ratio * iv.mpf(v)
                for v in [
                    old["available_from_counts"],
                    *old["held_flight_enclosed"][-1],
                ]
            ]
            delta = (
                arm["collision_acceptance"]
                - charge["arms"]["fresh"]["collision_acceptance"]
            )
            uncertainty = (
                144 * abs(delta) * initial["source_contrast_score_error_upper"] + 3e-14
            )
            held = []
            for interval, norm in zip(
                [meter - base[0], *[v - u for v, u in zip(read, base[1], strict=True)]],
                [1, *[float(max(abs(v))) for v in responses]],
                strict=True,
            ):
                width = uncertainty * norm + 3e-12
                held.append(reference.bounds(interval + iv.mpf([-width, width])))
            results[direction + "_" + key] = {
                "available_counts": k,
                "held_counts": shots,
                "source_only_prediction": [reference.bounds(v) for v in predicted],
                "observed_position_and_flight_change": held,
            }
    return {
        "seed": SEED,
        "starts_per_setting": N,
        "time": dt,
        "cells": list(cells),
        "gates": gates,
        "outcomes": results,
        "inputs": full.exp.save_inputs(
            "ARCHIVED_OVERLAP.md",
            "reference_arrival.json",
            "full_board_polarization.json",
        ),
    }


if __name__ == "__main__":
    result = run()
    full.exp.emit("modular_arrival_flight", result)
    for key, row in result["outcomes"].items():
        print(key, row["observed_position_and_flight_change"])
