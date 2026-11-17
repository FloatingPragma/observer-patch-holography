#!/usr/bin/env python3
"""Inspection-sized verification of the frozen native triangle protocol.

Only Python standard library. No file, network, subprocess or package access.
The initial optical ray, coherent eraser, reference and polar ticket coupling
are supplied. This checks existing finite arithmetic; physical row stays open.
"""
from fractions import Fraction as F
from math import log, sqrt
import json

# Actual recorded polar program: integer ticket counts; no libm regeneration.
WEIGHTS = (
    (32, 17, 15), (32, 17, 15), (32, 18, 14), (32, 18, 14),
    (32, 19, 13), (31, 20, 13), (31, 21, 12), (31, 21, 12),
    (30, 22, 12), (30, 23, 11), (29, 24, 11), (29, 24, 11),
    (28, 25, 11), (28, 25, 11), (27, 26, 11), (27, 27, 10),
    (26, 27, 11), (25, 28, 11), (25, 28, 11), (24, 29, 11),
    (24, 29, 11), (23, 30, 11), (22, 30, 12), (21, 31, 12),
    (21, 31, 12), (20, 31, 13), (19, 32, 13), (18, 32, 14),
    (18, 32, 14), (17, 32, 15), (17, 32, 15), (16, 32, 16),
    (15, 32, 17), (15, 32, 17), (14, 32, 18), (14, 32, 18),
    (13, 32, 19), (13, 31, 20), (12, 31, 21), (12, 31, 21),
    (12, 30, 22), (11, 30, 23), (11, 29, 24), (11, 29, 24),
    (11, 28, 25), (11, 28, 25), (11, 27, 26), (10, 27, 27),
    (11, 26, 27), (11, 25, 28), (11, 25, 28), (11, 24, 29),
    (11, 24, 29), (11, 23, 30), (12, 22, 30), (12, 21, 31),
    (12, 21, 31), (13, 20, 31), (13, 19, 32), (14, 18, 32),
    (14, 18, 32), (15, 17, 32), (15, 17, 32), (16, 16, 32),
    (17, 15, 32), (17, 15, 32), (18, 14, 32), (18, 14, 32),
    (19, 13, 32), (20, 13, 31), (21, 12, 31), (21, 12, 31),
    (22, 12, 30), (23, 11, 30), (24, 11, 29), (24, 11, 29),
    (25, 11, 28), (25, 11, 28), (26, 11, 27), (27, 11, 26),
    (27, 11, 26), (28, 11, 25), (28, 11, 25), (29, 11, 24),
    (29, 11, 24), (30, 11, 23), (30, 12, 22), (31, 12, 21),
    (31, 12, 21), (31, 13, 20), (32, 13, 19), (32, 14, 18),
    (32, 14, 18), (32, 15, 17), (32, 15, 17), (32, 16, 16),
)
EDGES = ((0, 1), (1, 2), (2, 0))
BASIS = ((1, -1, 0), (1, 1, -2))


def native_transfer(weights):
    """Independent endpoint averaging; count64 tickets and2 fair coins."""
    assert len(weights) == 3 and sum(weights) == 64 and min(weights) > 0
    counts = [[0] * 3 for _ in range(3)]
    flags = {}
    for initial in range(3):
        loads = [int(k == initial) for k in range(3)]
        for edge, (a, b) in enumerate(EDGES):
            for coin in (0, 1):
                total = loads[a] + loads[b]
                low, high = total // 2, -((-total) // 2)
                out = loads.copy()
                out[a], out[b] = (low, high) if coin == 0 else (high, low)
                current = loads[a] - out[a]
                final = out.index(1)
                # Fine flag distinguishes inputs sharing a source endpoint.
                outputs = flags.setdefault((edge, coin, current), set())
                assert final not in outputs
                outputs.add(final)
                counts[final][initial] += weights[edge]
    assert all(sum(counts[j][i] for j in range(3)) == 128 for i in range(3))
    assert counts == [list(row) for row in zip(*counts)]
    # Coordinate chart of an EXACT invariant plane. No physical projection.
    return tuple(tuple(F(sum(BASIS[i][a] * counts[a][b] * BASIS[j][b]
                             for a in range(3) for b in range(3)),
                         128 * (2 if i == 0 else 6))
                       for j in range(2)) for i in range(2))


def apply_product(matrices):
    """Initial f=(1,1,-2)/sqrt6; coefficients use unnormalized basis."""
    a, b = F(0), F(1)
    for matrix in matrices:
        a, b = (matrix[0][0] * a + matrix[0][1] * b,
                matrix[1][0] * a + matrix[1][1] * b)
    return a, b


def compare(matrices, reference):
    a, b = apply_product(matrices)
    c, d = apply_product([reference] * len(matrices))
    assert c == 0 and d > 0
    loop_norm = a * a / 3 + b * b
    total_norm = loop_norm + d * d
    visibility = 2 * b * d / total_norm
    minus = (1 - visibility) / 2
    assert 0 <= minus <= 1
    # Complete declared384-flag uniform eraser has success mapP/sqrt3.
    # All3^-N and optical attenuation factors remain in original success.
    success = total_norm / (2 * 3 ** len(matrices))
    log_probability = (log(success.numerator) - log(success.denominator)) / log(10)
    overlap = ((b > 0) - (b < 0)) / sqrt(1 + float(a * a / (3 * b * b)))
    return {
        "source_events_per_arm": len(matrices),
        "relative_sign_exact": (visibility > 0) - (visibility < 0),
        "source_overlap": overlap,
        "interference_visibility": float(visibility),
        "conditional_minus_probability": float(minus),
        "log10_original_success_probability": log_probability,
        "exact_visibility": str(visibility),
    }


def main():
    assert len(WEIGHTS) == 96 and WEIGHTS[-1] == (32, 16, 16)
    matrices = [native_transfer(w) for w in WEIGHTS]
    reference = native_transfer((32, 16, 16))
    one = compare(matrices, reference)
    two = compare(matrices + matrices, reference)
    fixed = compare([reference] * 96, reference)
    assert one["relative_sign_exact"] == -1
    assert two["relative_sign_exact"] == 1
    assert fixed["interference_visibility"] == 1
    print(json.dumps({
        "one_turn": one, "same_program_two_turns": two,
        "fixed_scheduler_control": fixed,
        "endpoint_columns_recomputed": 97 * 3 * 3 * 2,
        "supplied": ["exact optical preparation", "positive coherent memory",
                     "uniform coherent eraser", "path reference and analyzer",
                     "polar analyzer to native-ticket coupling"],
        "scope": "Existing finite transfer arithmetic, not a new experiment",
        "full_failure_outcome_replay": False,
        "physical_rotation_derived": False,
        "physical_row_passed": False,
    }, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
