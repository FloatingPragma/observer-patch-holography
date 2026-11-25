"""Forecast a separate material-force reading from native radar timestamps alone."""

import argparse
import json
from pathlib import Path

import numpy as np

from .. import station
from ...directional_stress.connection.force_gradients import force_effects
from . import measure

HERE = Path(__file__).resolve().parent


def preparations():
    plus = np.array([1, 0])
    minus = np.array([0, 1])
    x = np.array([1, 1]) / np.sqrt(2)
    y = np.array([1, 1j]) / np.sqrt(2)

    def product(left, right):
        vector = np.kron(left, right)
        return np.outer(vector, vector.conj())

    kept = product(plus, minus)
    dephased = product(x, y)
    other = product(y.conj(), x)
    return {
        "kept_Zplus_Zminus": np.array([kept, np.zeros_like(kept)]),
        "dephased_Xplus_Yplus": np.array([dephased / 2, dephased / 2]),
        "mixed_lottery": np.array([.75 * kept + .125 * other, .125 * other]),
    }


def event_effects(bins):
    board = station.iso.pairing.formation.federation(1)
    ports = station.iso.pairing.formation.region(board)
    effects = np.zeros((len(bins), 2, 4, 4), dtype=complex)
    for odd in (False, True):
        for charge in (-1, 1):
            h = station.iso.ranging.generator(board, ports, charge, 8, odd)
            kernels = station.iso.ranging.kernels(h, ports, .5, 12)
            hits = {r: station.iso.ranging.first_hits(kernels, *r) for r in measure.ROUTES}
            for i, (start, end, returned) in enumerate(bins):
                outbound, inbound = hits[start, end], hits[end, start]
                for tick in range(1, returned):
                    amplitude = inbound[returned - tick] @ outbound[tick]
                    effects[i, int(odd)] += (amplitude.conj().T @ amplitude) / 2
    return effects


def run():
    calibration = json.loads((HERE / "clock_scores.json").read_text())["whole_algebra_absolute_force"]
    bins = calibration["event_order"]
    clock = event_effects(bins)
    forces = np.array(force_effects(), dtype=complex)
    settings = {(i, j): np.array([k for k, row in enumerate(bins) if row[:2] == [i, j]], dtype=int)
                for i, j in measure.ROUTES}
    source_states = preparations()
    trials = 2**48
    rng = np.random.default_rng(191326)
    observations = {}
    for name, state in source_states.items():
        probs = np.einsum("abij,bji->a", clock, state).real
        observed_count = np.zeros(len(probs), dtype=np.int64)
        other_counts = {}
        # Each prepared trial has one category per selected route. All outcomes
        # outside the specified completed-return bins count in the final category.
        for (i, j), rows in settings.items():
            selected = probs[rows]
            rates = np.r_[selected, 1 - selected.sum()]
            assert min(rates) >= 0
            counts = rng.multinomial(trials, rates)
            observed_count[rows] = counts[:-1]
            other_counts[f"{i}-{j}"] = int(counts[-1])
        readings = {}
        for k, charge in enumerate(("vacancy", "excess")):
            frozen = calibration["signed_force"][charge]
            weights = np.array([int(v) for v in frozen["weight_numerators"]]) / int(frozen["denominator"])
            intercept = int(frozen["intercept_numerator"]) / int(frozen["denominator"])
            predicted = intercept + weights @ probs
            estimate = intercept + weights @ (observed_count / trials)
            independent = np.einsum("bij,bji", forces[k], state).real
            spread = sum(np.ptp(np.r_[0., weights[rows]])**2 for rows in settings.values())
            width = np.sqrt(spread * np.log(2 * (2 * len(source_states)) / .01) / (2 * trials))
            assert abs(independent - predicted) < frozen["absolute_error_all_normalized_heads_upper"] + 1e-13
            readings[charge] = {
                "force_from_separate_material_operator": float(independent),
                "forecast_from_population_clock_probabilities": float(predicted),
                "forecast_from_simulated_clock_counts": float(estimate),
                "simultaneous_99pct_half_width_including_calibration": float(
                    width + frozen["absolute_error_all_normalized_heads_upper"]),
            }
        observations[name] = {"kept_bin_counts": observed_count.tolist(), "all_other_counts": other_counts,
                              "force": readings}
    return {"source": "80-cell board, independently prepared heads; equal mixture of charge",
            "synthetic_draws_per_preparation_per_route": trials, "simultaneous_coverage": .99,
            "event_order": bins, "preparations": observations}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(result)
    else:
        print(result, end="")
