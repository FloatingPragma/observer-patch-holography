"""Resolve where each addressed pulse sends both admitted energy lines."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from scipy.sparse.linalg import eigsh

HERE = Path(__file__).resolve().parent
FIELD = HERE.parent


def pulses(path):
    cell = path.size - 1
    mask, cw = path.field.masks[0], path.field.clockwise[0]
    bits = [path.field.bits[e] for e in ((1, 2), (4, 5))]
    ring_rows = []
    head_signs = []
    for port, words in path.basis:
        word = words[cell]
        on_ring = word & mask in (cw, cw ^ mask)
        replacement = word ^ mask if on_ring else word
        ring_rows.append(path.ids[port, (*words[:cell], replacement)])
        head_signs.append((-1) ** sum(bool(word & bit) for bit in bits))
    ring_rows, head_signs = np.array(ring_rows), np.array(head_signs)
    np.testing.assert_array_equal(ring_rows[ring_rows], np.arange(len(path.basis)))
    np.testing.assert_array_equal(head_signs**2, 1)
    local = path.source(cell)
    np.testing.assert_array_equal(local[ring_rows], local @ np.diag([-1, 1]))
    head_swapped = local[:, [1, 0]]
    phase = int(np.sign(np.sum(head_swapped * local * head_signs[:, None])))
    np.testing.assert_array_equal(local * head_signs[:, None], phase * head_swapped)
    return ring_rows, head_signs


def probabilities(trial, norms, pulse_columns, error, sqrt_bounds):
    """Integer overlaps and outward probability intervals after the full pulse."""
    overlap = trial.astype(object).T @ pulse_columns.astype(object)
    records = []
    for incoming in range(2):
        intervals = []
        for outgoing in range(2):
            square = Fraction(
                int(overlap[outgoing, incoming]) ** 2, norms[outgoing] * norms[incoming]
            )
            lo, hi = sqrt_bounds(square)
            # Each test pairs two unit vectors, each no further away than error.
            intervals.append((max(0, lo - 2 * error) ** 2, min(1, hi + 2 * error) ** 2))
        leak = (
            max(0, 1 - sum(hi for _, hi in intervals)),
            min(1, 1 - sum(lo for lo, _ in intervals)),
        )
        assert leak[0] > 0
        records.append(
            {
                "incoming_line": incoming,
                "outcome_probabilities_lower_upper": [
                    [
                        max(0, float(np.nextafter(float(lo), -np.inf))),
                        min(1, float(np.nextafter(float(hi), np.inf))),
                    ]
                    for lo, hi in (*intervals, leak)
                ],
            }
        )
    return records


def run():
    sys.path.insert(0, str(FIELD.parents[4]))
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.channel_bound import (
        dyadic_vectors,
        sqrt_bounds,
    )
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.fixture import (
        FanPath,
        archived_fixture,
    )
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.measure import (
        conjugacy,
    )
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.spectral import (
        verified_inertia,
    )

    cut = Fraction(17, 4)
    rows = []
    with archived_fixture() as archive:
        for size in (3, 4):
            path = FanPath(archive, size, 1, asymmetric=True)
            neg = FanPath(archive, size, -1, asymmetric=True)
            h = path.h
            conjugacy(path, neg)
            assert path.board.cells == 320
            mask = sum(path.field.bits.values())
            opposite = np.array(
                [
                    neg.ids[port, tuple(mask ^ w for w in words)]
                    for port, words in path.basis
                ]
            )
            ring, head = pulses(path)
            neg_ring, neg_head = pulses(neg)
            np.testing.assert_array_equal(neg_ring[opposite], opposite[ring])
            np.testing.assert_array_equal(neg_head[opposite], head)

            counts = [verified_inertia(h, e) for e in (cut, Fraction(145, 32))]
            assert all(row["count_below"] == 2 for row in counts)
            rng = np.random.default_rng(923)
            energies, columns = eigsh(
                h, k=4, which="SA", v0=rng.standard_normal(h.shape[0]), tol=1e-13
            )
            trial, norms, centers, radii, error = dyadic_vectors(
                h, energies[:2], columns[:, :2], cut
            )
            # Both local actions are exact, but the root vector overlap uses the certified error.
            rows.append(
                {
                    "fans": size,
                    "physical_hubs": [fan[0] for fan in path.cells],
                    "inertia": counts,
                    "energy_centers": list(map(str, centers)),
                    "energy_radii": list(map(str, radii)),
                    "vector_distance_error_bound": str(error),
                    "native_first_ring_pi": probabilities(
                        trial, norms, trial[ring], error, sqrt_bounds
                    ),
                    "supplied_two_head_pi_phases": probabilities(
                        trial, norms, trial * head[:, None], error, sqrt_bounds
                    ),
                }
            )
    return {
        "board_cells": 320,
        "exact_full_generator_and_pulse_conjugacy_between_charge_signs": True,
        "energy_outcome_order": [
            "lower_retained_line",
            "upper_retained_line",
            "outside_retained_band",
        ],
        "cases": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as handle:
            handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        compare(
            result, json.loads((HERE / "energy_order.json").read_text()), "energy_order"
        )
        print("Verified the full-band energy order controls.")


if __name__ == "__main__":
    main()
