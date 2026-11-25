"""All-level formation bounds for the registered, blank unit host."""

import argparse
from fractions import Fraction
import json
from itertools import combinations
from math import ceil

import networkx as nx
import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from .bounds import layout_type, small_bound
from .measure import HERE, word_table


def verify_native_seam():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    blank = (exchange + 8 * contrast)[: law.N_SRC, : law.N_SRC]
    for i, (left, right) in enumerate(law.SRC):
        cap = blank[i, i] + np.sum(abs(np.delete(blank[i], i)))
        np.testing.assert_allclose(
            cap, 8 * (left - right) ** 2 + int(left != right), rtol=0, atol=1e-13
        )


def local_caps():
    """Gershgorin row costs in twice-energy units, split over cell owners."""
    _, _, words = word_table(12)
    s = words.astype(np.int16) - 1
    d = np.sum(s**2, axis=1)
    q = np.sum(s, axis=1)
    cap = np.zeros(len(s), dtype=np.int16)
    internal = np.zeros(len(s), dtype=np.int16)
    for a, b in cell().seams:
        delta = s[:, a] - s[:, b]
        internal += delta**2
        cap += 16 * delta**2 + 2 * (delta != 0)
    cap += 16 * np.sum(s**2 + 2 * s, axis=1).astype(np.int16)
    # Turning a degree-five port into a gateway removes its reference term
    # and adds its half of the external seam bound.
    advantage = 33 * s**2 - 16 * (s**2 + 2 * s)
    # Maximise over all triples, even those absent on a given board.
    n_minus = (d - q) // 2
    n_zero = 12 - d
    gate_minus = np.minimum(n_minus, 3)
    gate_plus = np.maximum(0, 3 - gate_minus - n_zero)
    arbitrary = cap + 49 * gate_minus - 15 * gate_plus
    # An independent lower envelope drops the positive external contrasts and
    # removes the three largest reference corrections. The blank exchange is PSD.
    n_plus = (d + q) // 2
    cheapest = (
        internal
        + 3 * n_plus
        - n_minus
        - (3 * np.minimum(3, n_plus) - np.maximum(0, n_minus - 9))
    )
    assert np.min(4 * cheapest - 6 * q - 3 * d) == 0
    costs = np.full(169, 10000, dtype=int)
    np.minimum.at(costs, 13 * n_minus + n_plus, cheapest)
    few = small_bound(costs.reshape(13, 13))
    floor = min(8 * min(few.values()), 6 * 9 - 12)
    boards = [federation(level) for level in range(3)]
    layouts = {
        layout_type(tuple(np.flatnonzero(board.degree[i : i + 12] == 6)))
        for board in boards
        for i in range(0, board.n, 12)
    }
    actual = np.maximum.reduce(
        [cap + advantage[:, gate].sum(axis=1) for gate in sorted(layouts)]
    )
    assert np.all(actual <= arbitrary)
    classes = 25 * d + q + 12

    def group(cost):
        maxima = np.full(325, -10000, dtype=int)
        np.maximum.at(maxima, classes, cost)
        return [
            (i // 25, i % 25 - 12, Fraction(int(value), 2))
            for i, value in enumerate(maxima)
            if value > -10000 and i // 25
        ]

    return boards, sorted(layouts), group(actual), group(arbitrary), floor, few


def certificate(rows, per_defect):
    """Find the exact lower envelope of valid tilts E <= A*count + B*charge."""
    # At fixed signed-charge/count slope, only the largest cost/count matters.
    lines = {}
    for defects, charge, cap in rows:
        count = defects if per_defect else 1
        slope = Fraction(charge, count)
        lines[slope] = max(lines.get(slope, cap / count), cap / count)
    # The convex upper envelope is minimized at an intersection (or on a flat).
    tilts = {Fraction(0)}
    for m, u in lines.items():
        for n, v in lines.items():
            if m < n:
                tilts.add((u - v) / (m - n))

    def objective(tilt):
        return max(v - tilt * m for m, v in lines.items()), abs(tilt), tilt

    cap, _, tilt = min(map(objective, tilts))
    assert all(
        value <= cap * (d if per_defect else 1) + tilt * q for d, q, value in rows
    )
    return cap, tilt


def high_line_witness(board):
    """A Rayleigh lower bound on the top level in one specified composition."""
    internal = cell().seams

    def cut(positive):
        return sum((u in positive) != (v in positive) for u, v in internal)

    # A neutral six-plus, six-minus word maximizing the internal cut.
    word = max(combinations(range(12), 6), key=cut)
    owner = -np.ones(12, dtype=np.int64)
    owner[list(word)] = 1
    owners = nx.Graph()
    owners.add_edges_from((a // 12, b // 12) for a, b in board.external)
    distances = nx.single_source_shortest_path_length(owners, 0)
    order = sorted(distances, key=lambda i: (distances[i], i))
    signed = np.zeros(board.n, dtype=np.int64)
    edges = np.array(board.seams)
    # A vacancy at a degree-five port gives true E_light <= 32 + 5/2.
    den = 69
    rows = {}
    for occupied, i in enumerate(order, start=1):
        signed[12 * i : 12 * i + 12] = owner
        if occupied == 1:
            signed[12 * i + word[0]] = 0
        delta = signed[edges[:, 0]] - signed[edges[:, 1]]
        contrast = delta @ delta + (6 - board.degree) @ (signed**2 + 2 * signed)
        twice_energy = int(16 * contrast + np.count_nonzero(delta))
        assert signed.sum() == -1 and signed @ signed == 12 * occupied - 1
        for ratio in (207, 3477):
            if twice_energy >= ratio * den and str(ratio) not in rows:
                rows[str(ratio)] = {
                    "occupied_cells_in_trial": occupied,
                    "defects": 12 * occupied - 1,
                    "twice_trial_energy": twice_energy,
                    "certified_top_line_over_lightest": str(
                        Fraction(twice_energy, den)
                    ),
                }
    return {
        "balanced_owner_positive_ports": list(word),
        "internal_cut": cut(word),
        "twice_light_trial_energy": den,
        "first_crossing_by_rounded_ratio": rows,
    }


def run():
    verify_native_seam()
    boards, layouts, actual, arbitrary, floor, few = local_caps()
    out = {
        "energy_reference": "registered blank unit-filled host, total excess -1",
        "owner_words_per_layout": 3**12,
        "actual_gateway_classes": layouts,
    }
    assert (
        floor
        == json.loads((HERE / "bounds.json").read_text())["arbitrary_three_gateways"][
            "universal_energy_floor"
        ]
        == 32
    )
    out["verified_negative_sector_floor"] = floor
    # Both sides are exact integers for every pair of adjacent sharp port records.
    assert all(
        2 * (8 * (s - t) ** 2 + int(s != t)) <= 33 * (s**2 + t**2)
        for s in (-1, 0, 1)
        for t in (-1, 0, 1)
    )
    for domain, rows in (
        ("actual_gateways", actual),
        ("any_three_gateways", arbitrary),
    ):
        a, b = certificate(rows, True)
        g, h = certificate(rows, False)
        low_fixed = {d: 8 * c for d, c in few.items()}
        low_fixed[9] = 6 * 9 - 12
        fixed_ratios = {d: (a * d - b) / lower for d, lower in low_fixed.items()}
        # For d >= 9 the derivative of (a*d-b)/(6*d-12) has numerator 6*b-12*a.
        assert 6 * b < 12 * a
        fixed_cap = max(fixed_ratios.values())
        witness = {}
        # A declared target ratio only sets the comparison threshold.
        for ratio in (Fraction("206.7682827"), Fraction(3477)):
            min_d = (floor * ratio + b) / a
            min_c = (floor * ratio + h) / g
            d_fixed = ceil(min_d)
            d_fixed += 1 - d_fixed % 2
            witness[str(float(ratio))] = {
                "minimum_mean_defect_ports": float(min_d),
                "minimum_defects_in_fixed_composition": d_fixed,
                "minimum_mean_occupied_cells": float(min_c),
                "minimum_cells_for_a_strict_dressing_support": ceil(min_c),
            }
        out[domain] = {
            "per_defect_and_charge_coefficients": [str(a), str(b)],
            "per_occupied_owner_and_charge_coefficients": [str(g), str(h)],
            "every_single_composition_excited_to_ground_ratio_ceiling": str(fixed_cap),
            "small_defect_and_tail_ratio_envelopes": {
                str(d): str(r) for d, r in fixed_ratios.items()
            },
            "necessary_dressing_sizes": witness,
            "all_level_ratio_ceiling_by_board_cells": {
                str(board.cells): float((g * board.cells - h) / floor)
                for board in boards
            },
        }
    out["high_line_trial_witness_by_board_cells"] = {
        str(b.cells): high_line_witness(b) for b in boards
    }
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE / "excited_bounds.json"
    result = json.loads(json.dumps(run()))
    if parser.parse_args().write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif result != json.loads(path.read_text()):
        raise ValueError("All-level bounds differ from the recorded certificate")
    print(json.dumps(result["actual_gateways"], indent=2))


if __name__ == "__main__":
    main()
