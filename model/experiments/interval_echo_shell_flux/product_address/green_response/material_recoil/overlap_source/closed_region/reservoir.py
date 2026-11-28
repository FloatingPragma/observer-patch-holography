"""Microcanonical contact thermometry, preserving both native load and energy."""

import argparse
from functools import cache
from fractions import Fraction
import json
from math import comb, log

import numpy as np

from .. import measure as flights
from .measure import HERE, check
from .modular import EQUAL, EVEN, ODD, H, calibrate


@cache
def counts(n, energy, load):
    """Count sharp energy eigenwords with native energies (0,8,9,0) and loads (0,1,1,2)."""
    if min(n, energy, load) < 0 or load > 2 * n or energy > 9 * n:
        return 0
    total = 0
    for odd in range(energy % 8, energy // 9 + 1, 8):
        even = (energy - 9 * odd) // 8
        singles = even + odd
        double, rem = divmod(load - singles, 2)
        if rem or double < 0 or singles + double > n:
            continue
        total += comb(n, singles) * comb(singles, odd) * comb(n - singles, double)
    return total


def log_interval(a, b):
    """Outward rational bounds on a log ratio using a short atanh expansion."""
    if a < b:
        lo, hi = log_interval(b, a)
        return -hi, -lo
    grid = 10**12
    x = Fraction(a - b, a + b)
    floor = int(x * grid)
    left, right = Fraction(floor, grid), Fraction(floor + 1, grid)
    assert 0 <= left <= x <= right < 1

    def partial(v):
        return 2 * sum(v**j / j for j in range(1, 50, 2))

    return partial(left), partial(right) + 2 * right**51 / (51 * (1 - right**2))


def run():
    held = json.loads(HERE.parent.joinpath("formed_heads.json").read_text())
    pilots = EVEN - (8 * ODD + EQUAL) / 9
    assert abs(np.sum(H * pilots)) < 3e-15
    rows = []
    for m in [2**j for j in range(4, 13)]:
        total = counts(m, 2 * m, m)
        weights = [
            counts(m - 1, 2 * m - e, m - q) for e, q in ((0, 0), (8, 1), (9, 1), (0, 2))
        ]
        assert total == sum(weights) and weights[0] == weights[-1]
        assert 8 * weights[1] + 9 * weights[2] == 2 * total
        p = [v / total for v in weights]
        row = dict(contacts=m, local_weights=p, singular=0 in weights)
        if not row["singular"]:
            # Log ratios avoid subtracting large integer logarithms. A uniform
            # one-contact microcanonical energy and a faithful density are different requirements.
            beta8 = log(weights[0] / weights[1]) / 8
            beta1 = log(weights[1] / weights[2])
            low8, high8 = (v / 8 for v in log_interval(weights[0], weights[1]))
            low1, high1 = log_interval(weights[1], weights[2])
            assert low8 > 0
            ratios = [a / b for a in (low1, high1) for b in (low8, high8)]
            bounds = [(x - 1) / 17 * 10**9 for x in (min(ratios), max(ratios))]
            ref = p[0] * 2 * EQUAL + p[1] * EVEN + p[2] * ODD
            controls = calibrate(ref, pilots)
            np.testing.assert_allclose(
                [beta8, beta1],
                [controls["beta_equal_to_even"], controls["beta_even_to_odd"]],
                atol=3e-13,
            )
            row.update(
                fractional_error_billionths=[
                    int(bounds[0] // 1),
                    int(-((-bounds[1]) // 1)),
                ],
                beta_equal_to_even=beta8,
                beta_even_to_odd=beta1,
                matched_energy_modular=float(controls["held_modular"]),
                source=[],
            )
            for s in held["source"]:
                if s["orientation"] != "matching":
                    continue
                # After the common fair flip, sharp unlike loads put half their
                # weight on each parity. Only the total mismatch changes.
                diff = s["independent_unlike"] - s["formed_unlike"]
                result = calibrate(ref, diff * ((EVEN + ODD) / 2 - EQUAL))
                predicted = result["held_modular"] / beta8
                row["source"].append(
                    dict(
                        source_sign=s["source_sign"],
                        native_erasure_energy=result["held_energy"],
                        modular_predicted_energy=predicted,
                        fractional_error=predicted / result["held_energy"] - 1,
                    )
                )
        rows.append(row)
    outputs = []
    # For a single fresh preparation, a null ticket gives the same downstream
    # bank in both arms. A click gives one identical independently routed pulse.
    visits = flights.receiver.visits_before_drain()
    distal = []
    for probe in (1, 2):
        setting = flights.receiver.Setting("matched", (6, 5, 5, 0), probe, (2, 2, 2))
        probs = flights.receiver.table(setting)
        reward = sum(
            (p[words] - 0.5) / 2
            for p, words in zip(
                probs, flights.receiver.parent.addresses(setting.face), strict=True
            )
        )
        per_click = float((visits * reward).sum())
        np.testing.assert_allclose(
            per_click, -1.6 * (visits[5, 4, 4] - visits[6, 4, 4]) / 8
        )
        distal.append(
            dict(
                probe=probe,
                signed_endpoint_per_prep_each_bank=per_click * 8 / 512,
                predicted_difference_from_masked_log=per_click / 576,
                physical_difference=0,
            )
        )
    for s in held["source"]:
        if s["orientation"] != "matching":
            continue
        p = np.array(s["formed_sharp_pairs"])
        ref = np.diag((p + p[::-1]) / 2)
        v = calibrate(ref, pilots)
        w = s["formed_unlike"]
        np.testing.assert_allclose(v["held_modular"], log((1 - w) / w) / 9)
        implied = v["held_modular"] / v["beta_equal_to_even"]
        np.testing.assert_allclose(implied, 8 / 9)
        outputs.append(
            dict(
                source_sign=s["source_sign"],
                energy_of_either=8,
                ticket_probability_of_either=8 / 512,
                modular_difference=v["held_modular"],
                energy_difference_inferred_from_masked_log=implied,
            )
        )
    return dict(
        held_equal_energy_challenge=outputs,
        material_endpoints=distal,
        contact_bank=rows,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    path = HERE / "reservoir.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), data)
