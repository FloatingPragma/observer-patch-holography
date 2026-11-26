"""Produce complete finite observer records; decode only the observed terminal counts."""

import argparse
import json
from itertools import product
from math import log, pi, sqrt
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_consensus.signed_effects.sequential import (
    grid_counts,
)
from .apparatus import ANGLES, K, LAW, chamber, contact, free_records, layout, writer

HERE = Path(__file__).parent
SHOTS, SEED = 131072, 961126
PAIRS = ((0, 1), (1, 2), (0, 2))
MENUS = ("native", "square", "exchange_only")


def banks():
    for menu, direction, p, (a, b), bias in product(
        MENUS, range(2), range(12), PAIRS, range(2)
    ):
        yield (menu, direction, -1, p, a, b, bias)
    for turned, p, bias in product(range(12), range(12), range(2)):
        yield ("native", 0, turned, p, 0, 1, bias)


def flights():
    return product(("native", "park_exchange"), (0, 1), (0, 1))


def acquire():
    grid = np.random.default_rng(SEED)
    physical = layout()
    swap = np.eye(9)[[LAW.SRC.index((b, a)) for a, b in LAW.SRC]]
    for k in (0, K):
        np.testing.assert_allclose(contact(k, pi), swap, rtol=0, atol=1e-13)
    for a, b in LAW.SRC:
        # Signed heads on unit checks really read n-1, both coins included.
        column = LAW.write()[:, LAW.SRC.index((a, 1))]
        for index in np.flatnonzero(column):
            assert LAW.LABELS[index // 9][1] == a - 1
    rows = []
    for key in banks():
        menu, direction, turned, p, a, b, bias = key
        arm = physical["stations"][direction]["arms"][p]
        assert len(arm["read"]) == 3
        circuit = tuple(
            writer(angle + (pi / 2 if j == turned else 0), menu)
            for j, angle in enumerate(ANGLES)
        )
        probability = chamber(a, b).receive(p, circuit, menu, bias)
        rows.append(
            [*key, grid_counts(grid, SHOTS, probability, tolerance=1.3e-11).tolist()]
        )
    flight = []
    for menu, turn, tick in flights():
        probability = free_records(turn * pi, tick / 16, menu == "native")
        flight.append(
            [
                menu,
                turn,
                tick,
                grid_counts(grid, SHOTS, probability, tolerance=1.3e-11).tolist(),
            ]
        )
    return dict(shots=SHOTS, seed=SEED, local=rows, flights=flight)


def unpack(rows, keys, length):
    output = {}
    for *key, counts in rows:
        key = tuple(key)
        if (
            key in output
            or len(counts) != length
            or any(type(n) is not int or n < 0 for n in counts)
            or sum(counts) != SHOTS
        ):
            raise ValueError("Invalid or duplicate record bank")
        output[key] = np.array(counts) / SHOTS
    if set(output) != set(keys):
        raise ValueError("Incomplete or unexpected record bank")
    return output


def decode(raw):
    if raw["shots"] != SHOTS or raw["seed"] != SEED:
        raise ValueError("Unexpected enrollment")
    data = unpack(raw["local"], banks(), 27)
    free = unpack(raw["flights"], flights(), 3)

    def z(menu, direction, turned, p, pair):
        a, b = pair
        rates = [
            data[menu, direction, turned, p, a, b, bias][9 * a + 3 * b + 1]
            for bias in range(2)
        ]
        return complex(2 * rates[0] - 1, 1 - 2 * rates[1])

    # Both quadratures use the same recorded load word.
    def freq_error(observed, predicted):
        residual = observed - predicted
        return max(abs(residual.real), abs(residual.imag)) / 2

    errors = {}
    for menu in MENUS:
        held = []
        for direction, p in product(range(2), range(12)):
            base = z(menu, direction, -1, p, (0, 1))
            held.extend(
                [
                    freq_error(z(menu, direction, -1, p, (1, 2)), base),
                    freq_error(z(menu, direction, -1, p, (0, 2)), base**2),
                ]
            )
        errors[menu] = max(held)
    independent = max(
        freq_error(
            z("native", 0, turned, p, (0, 1)),
            z("native", 0, -1, p, (0, 1)) * (1j if p == turned else 1),
        )
        for p, turned in product(range(12), repeat=2)
    )
    guard_bad = {
        menu: max(
            1 - sum(counts[1::3]) for key, counts in data.items() if key[0] == menu
        )
        for menu in MENUS
    }
    flight_differences = {
        menu: [float(free[menu, 1, t][0] - free[menu, 0, t][0]) for t in (0, 1)]
        for menu in ("native", "park_exchange")
    }
    eps = sqrt(log(2 * (len(data) * 27 + len(free) * 3) / 0.01) / (2 * SHOTS))
    held_bound = eps + 2 * sqrt(2) * eps + 4 * eps**2
    assert (
        errors["native"] < held_bound < min(errors["square"], errors["exchange_only"])
    )
    assert independent < 2 * eps and guard_bad["native"] == 0
    assert (
        abs(flight_differences["native"][1]) > 2 * eps
        and abs(flight_differences["native"][0]) < 2 * eps
    )
    assert max(map(abs, flight_differences["park_exchange"])) < 2 * eps
    return dict(
        banks=len(data) + len(free),
        starts=SHOTS * (len(data) + len(free)),
        simultaneous_frequency_radius=eps,
        maximum_held_rate_residual=errors,
        held_bound=held_bound,
        maximum_addressed_control_residual=independent,
        maximum_nonunit_guard_fraction=guard_bad,
        phase_change_in_vacancy_at_zero=flight_differences,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    raw = acquire()
    summary = decode(raw)
    path = HERE / "records.json"
    if args.write:
        path.write_text(json.dumps(raw, separators=(",", ":")) + "\n")
    elif json.loads(path.read_text()) != raw:
        raise ValueError("Observer record replay differs")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
