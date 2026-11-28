"""Repeat only whole-cell energy readings and compare both finite-model populations."""

import argparse
from copy import deepcopy
from fractions import Fraction
import json
import random

import numpy as np

from matched_einstein_limit.source_implication.full_algebra.correlation_rank.spatial_cuts.chart_transport.adaptive_optics.counts import (
    random_binomial,
)
from . import regional_rate as region

RECEIPT = region.HERE / "regional_precision.json"
SHOTS = 2**94
SEED = 190823


def population(states, steps, board, route, rings, operators, strength):
    """Use the full receiver laws with the same finite differences as the counted experiment."""
    camera = region.camera
    energy = {k: [] for k in states}
    images = {k: {"forward": [], "reverse": []} for k in states}
    for offset in (-camera.OFFSET, 0, camera.OFFSET):
        t = region.CENTER + offset
        if offset:
            effect = region.energy_effect(board, route, rings, t, strength, operators)
            for key, state in states.items():
                mean = region.measure.station.pair(effect, state)
                energy[key].append(float((1 + mean / region.SCALE) / 2))
        for tag, screen in (("forward", rings), ("reverse", rings[::-1])):
            table = camera.BEAM.effects(board, route, screen, camera.Q, t, strength)
            for key, state in states.items():
                pixels = camera.BEAM.probabilities(table, state)
                images[key][tag].append(camera.BEAM.image(pixels)[0])

    rates, focus = {}, {}
    for key, panel in images.items():
        # A population clock inversion on the calibrated branch gives the scheduled times.
        rates[key] = region.SCALE * (energy[key][1] - energy[key][0]) / camera.OFFSET
        f = {}
        for tag, (left, center, right) in panel.items():
            acceleration = (left - 2 * center + right) / camera.OFFSET**2
            f[tag] = -np.trace(np.linalg.solve(center, acceleration))
        focus[key] = f["forward"] - f["reverse"]

    def contrast(data, key):
        return float((data[key + "_plus"] - data[key + "_minus"]) / (2 * steps[key]))

    ordinary_rate = contrast(rates, "ordinary")
    ordinary_f = contrast(focus, "ordinary")
    comparisons = {}
    for key in region.STOCKS[1:]:
        observed = contrast(focus, key)
        rate = contrast(rates, key)
        forecast = ordinary_f / ordinary_rate * rate
        comparisons[key] = dict(
            regional_rate=rate,
            predicted_focusing=forecast,
            measured_focusing=observed,
            measured_minus_predicted=observed - forecast,
        )
    return dict(
        energy_plus_probabilities=energy,
        ordinary_rate=ordinary_rate,
        ordinary_focusing=ordinary_f,
        comparisons=comparisons,
    )


def decode(primary, data):
    if data["shots_per_energy_setting"] != SHOTS or data["seed"] != SEED:
        raise ValueError("Unregistered precision acquisition")
    bank = deepcopy(primary["banks"]["8"]["counts"])
    if set(data["energy_plus_counts"]) != set(bank["stocks"]):
        raise ValueError("Incomplete set of energy settings")
    for key, counts in data["energy_plus_counts"].items():
        bank["stocks"][key]["energy_plus"] = counts
    return region.decode(bank, energy_shots=SHOTS)


def acquire(primary):
    states, steps, recipes = region.preparations()
    if (
        recipes != primary["endpoint_recipes"]
        or steps != primary["banks"]["8"]["counts"]["half_steps"]
    ):
        raise ValueError("Precision bank would use different source recipes")
    formation = region.camera.BEAM.front.optics.pairing.formation
    board = formation.federation(1)
    route = formation.region(board)
    rings = region.camera.BEAM.screens(board, route)
    operators = [region.regional_difference(board, route, odd) for odd in (False, True)]
    laws = {
        str(k): population(states, steps, board, route, rings, operators, k)
        for k in (8, 0)
    }
    rng = random.Random(SEED)
    counts = {
        name: [random_binomial(SHOTS, Fraction.from_float(prob), rng) for prob in probs]
        for name, probs in laws["8"]["energy_plus_probabilities"].items()
    }
    data = {
        "primary_sha256": region.digest(region.RECEIPT),
        "precision_protocol_sha256": region.digest(region.HERE / "PROTOCOL.md"),
        "shots_per_energy_setting": SHOTS,
        "seed": SEED,
        "energy_plus_counts": counts,
        "population": laws,
    }
    data["decoded"] = decode(primary, data)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--counts-only", action="store_true")
    args = parser.parse_args()
    primary = json.loads(region.RECEIPT.read_text())
    if args.counts_only:
        data = json.loads(RECEIPT.read_text())
        if data["primary_sha256"] != region.digest(region.RECEIPT) or data[
            "precision_protocol_sha256"
        ] != region.digest(region.HERE / "PROTOCOL.md"):
            raise ValueError("Wrong primary acquisition or precision protocol")
        if decode(primary, data) != data["decoded"]:
            raise ValueError("Precision intervals differ from receipt")
    else:
        data = acquire(primary)
        if args.write:
            with RECEIPT.open("x") as handle:
                json.dump(data, handle, indent=2, allow_nan=False)
                handle.write("\n")
        elif data != json.loads(RECEIPT.read_text()):
            raise ValueError("Recomputed acquisition differs from receipt")
    print(
        json.dumps(
            {
                "decoded": data["decoded"],
                "population": {
                    k: p["comparisons"] for k, p in data["population"].items()
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
