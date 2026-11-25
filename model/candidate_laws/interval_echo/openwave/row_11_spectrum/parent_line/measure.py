"""Exact local moments of the archived sharp variational parents on 320 cells."""

import argparse
from collections import deque
from fractions import Fraction
import json
from math import lcm, sqrt
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_04_compact.gateway_bound import (
    dark_projector,
)
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
TRIALS = HERE.parent / "excited_bounds.json"
GROUND = HERE.parent / "lightest/results.json"
RATIOS = (
    HERE.parents[1] / "electron_field/lepton_family/arrival_spin/anchor_spectrum.json"
)


def owner_order(board):
    """Recover the source's distance ordering with deterministic ties."""
    neighbours = [[] for _ in range(board.cells)]
    for a, b in board.external:
        u, v = a // 12, b // 12
        neighbours[u].append(v)
        neighbours[v].append(u)
    distances = {0: 0}
    queue = deque([0])
    while queue:
        owner = queue.popleft()
        for target in neighbours[owner]:
            if target not in distances:
                distances[target] = distances[owner] + 1
                queue.append(target)
    assert len(distances) == board.cells
    layers = [[] for _ in range(max(distances.values()) + 1)]
    for owner in range(board.cells):
        layers[distances[owner]].append(owner)
    return [owner for layer in layers for owner in layer]


def native_word(board, positive, owners):
    x = np.ones(board.n, dtype=int)
    cell_word = np.zeros(12, dtype=int)
    cell_word[positive] = 2
    for owner in owners:
        x[12 * owner : 12 * (owner + 1)] = cell_word
    x[12 * owners[0] + positive[0]] = 1
    return x


def sharp_moments(board, word, exchange=1):
    """Formation energy in half-units and squared residual in quarter-units."""
    k = int(candidate.KAPPA)
    assert k == candidate.KAPPA
    edges = np.array(board.seams)
    a, b = edges.T
    assert len(set(board.seams)) == len(edges) and np.all(a < b)
    diff = word[a] - word[b]
    n_swap = int(np.count_nonzero(diff))
    potential = int(k * (diff @ diff + (6 - board.degree) @ (word**2 - 1)))
    return {
        "different_load_edges": n_swap,
        "external_different_load_edges": int(
            np.count_nonzero((a // 12 != b // 12) & (diff != 0))
        ),
        "formation_potential": potential,
        "twice_mean": 2 * potential + exchange * n_swap,
        "four_times_variance": exchange**2 * n_swap,
    }


def compact_light(board):
    # An ordinary one-vacancy Hamiltonian in the host-relative convention.
    # Construct a compact vector on a *320-cell* chart and verify every row.
    chart = tuple(map(int, np.flatnonzero(board.degree[:12] == 6)))
    projector = dark_projector(chart)
    index = next(i for i in range(12) if projector[i, i])
    column = projector[:, index]
    scale = lcm(*(int(v.q) for v in column))
    v = np.zeros(board.n, dtype=int)
    v[:12] = np.array(scale * column, dtype=int).ravel()
    a, b = np.array(board.seams).T
    twice_hv = (board.degree + 2 * int(candidate.KAPPA) * (2 * board.degree - 6)) * v
    np.add.at(twice_hv, a, -v[b])
    np.add.at(twice_hv, b, -v[a])
    twice_energy = int(twice_hv[index] // v[index])
    np.testing.assert_array_equal(twice_hv, twice_energy * v)
    assert twice_energy == 70 and not v[list(chart)].any()
    return {
        "twice_energy": twice_energy,
        "owner": 0,
        "gateways": list(chart),
        "local_integer_amplitudes": v[:12].tolist(),
        "rows_verified": board.n,
    }


def window_bound(mean, variance, lo, hi):
    """Cantelli bound for a closed energy interval, valid on either side."""
    if lo > hi:
        raise ValueError("Energy interval is reversed")
    distance = max(lo - mean, mean - hi, 0)
    if not distance:
        return Fraction(1)
    return variance / (variance + distance**2)


def run():
    board = federation(2)
    record = json.loads(TRIALS.read_text())
    recipes = record["high_line_trial_witness_by_board_cells"][str(board.cells)]
    ratios = json.loads(RATIOS.read_text())["scale_audit"]["closed_form_ratios"]
    light = compact_light(board)
    anchor = Fraction(light["twice_energy"], 2)
    floor = record["verified_negative_sector_floor"]
    ceiling = Fraction(recipes["twice_light_trial_energy"], 2)
    ground = json.loads(GROUND.read_text())
    ticks = Fraction(ground["native_units_per_tick"])
    ground_row = next(row for row in ground["boards"] if row["cells"] == board.cells)
    low, high = [v * ticks for v in ground_row["ground_interval_ticks"]]
    assert floor <= low <= high <= ceiling
    order = owner_order(board)
    single = np.ones(board.n, dtype=int)
    single[np.flatnonzero(board.degree == 5)[0]] = 0
    sharp_light = sharp_moments(board, single)
    assert sharp_light["twice_mean"] == recipes["twice_light_trial_energy"]
    panels = {}
    for key, ratio_digits in zip(("207", "3477"), ratios[1:], strict=True):
        recipe = recipes["first_crossing_by_rounded_ratio"][key]
        packed = order[: recipe["occupied_cells_in_trial"]]
        word = native_word(board, recipes["balanced_owner_positive_ports"], packed)
        assert np.sum(word - 1) == -1
        assert np.sum(word != 1) == recipe["defects"]
        panel = sharp_moments(board, word)
        assert panel["twice_mean"] == recipe["twice_trial_energy"]
        mean = Fraction(panel["twice_mean"], 2)
        var = Fraction(panel["four_times_variance"], 4)
        ratio = Fraction(f"{float(ratio_digits):.6f}")
        # This interval also covers the unrounded ratio in the source receipt.
        rounding = Fraction(1, 2_000_000)
        assert abs(ratio - Fraction(ratio_digits)) < rounding
        target = anchor * ratio
        point = window_bound(mean, var, target, target)
        window = window_bound(mean, var, target - 1, target + 1)
        lightest = window_bound(
            mean, var, low * (ratio - rounding), high * (ratio + rounding)
        )
        panels[key] = {
            "occupied_cells": len(packed),
            "defects": recipe["defects"],
            **panel,
            "raw_second_moment": str(mean**2 + var),
            "standard_deviation": sqrt(var),
            "width_over_centroid": sqrt(var) / float(mean),
            "return_probability_second_derivative_at_zero": str(-2 * var),
            "comparison_ratio_six_decimals": str(float(ratio)),
            "conditional_target": str(target),
            "target_minus_mean": str(target - mean),
            "exact_rounded_target_weight_upper": str(point),
            "one_unit_half_width_weight_upper": str(window),
            "ground_target_interval": [
                str(low * (ratio - rounding)),
                str(high * (ratio + rounding)),
            ],
            "ground_target_interval_weight_upper": str(lightest),
            "ground_one_unit_half_width_weight_upper": str(
                window_bound(
                    mean,
                    var,
                    low * (ratio - rounding) - 1,
                    high * (ratio + rounding) + 1,
                )
            ),
            "double_exchange_control": sharp_moments(board, word, exchange=2),
        }
    return {
        "cells": board.cells,
        "ports": board.n,
        "seams": len(board.seams),
        "energy_reference": "blank unit-filled host; native zero-load degree references",
        "compact_reference": light,
        "sharp_vacancy": sharp_light,
        "filled_sham": sharp_moments(board, np.ones(board.n, dtype=int)),
        "certified_lightest_interval": [str(low), str(high)],
        "trials": panels,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE / "results.json"
    result = run()
    if parser.parse_args().write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif result != json.loads(path.read_text()):
        raise RuntimeError("Sharp-parent moment receipt differs")
    print(json.dumps(result["trials"], indent=2))


if __name__ == "__main__":
    main()
