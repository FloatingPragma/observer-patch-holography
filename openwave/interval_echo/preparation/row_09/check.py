"""Recompute the row 9 contrasts from the archived observer records."""

import json
from itertools import combinations, product
from math import fsum, isclose, isfinite
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
EVIDENCE = ROOT / "model/candidate_laws/interval_echo/openwave/row_09_neutral"


def read(name):
    return json.loads((EVIDENCE / name).read_text(encoding="utf-8"))


def passive_decay():
    data = read("passive_decay/records.json")
    shots = data["original_starts"]
    outcomes = {}
    for arm in data["arms"]:
        for appointment in arm["settings"]:
            matrix = appointment["joint"]
            bins = [count for row in matrix for count in row]
            if (
                len(matrix) != 4
                or any(len(row) != 4 for row in matrix)
                or any(type(n) is not int or n < 0 for n in bins)
            ):
                raise ValueError("Invalid core and receiver outcome bank")
            if sum(bins) != shots:
                raise ValueError("Incomplete original-start count bank")
            if appointment["wait"] == 16:
                outcomes[arm["source"], arm["joined"]] = matrix[2][1]
    native = outcomes["parent", True]
    controls = outcomes["parent", False], outcomes["ground", True]
    ratio = data["reservoir_levels"][1] / data["source_levels"][3]
    if not (native / shots > 0.29 and controls == (0, 0) and 0 < ratio < 0.1):
        raise ValueError("Passive decay reproduction differs")
    print(
        f"Joined material: {native:,}/{shots:,}, controls {controls}, "
        f"receiver/parent formation energy {ratio:.6f}"
    )


def captured_soft(arm):
    retain = arm["any_excited"]["core_by_outlet"]
    captured = [
        arm[name]["core_by_outlet"] for name in ("return_1", "return_2", "return_3")
    ]
    for matrix in [retain, *captured]:
        if len(matrix) != 4 or any(
            len(row) != 5 or not all(0 <= n <= 1 for n in row) for row in matrix
        ):
            raise ValueError("Invalid capture / core / outlet outcome bank")
    if (
        abs(fsum(map(fsum, retain)) - fsum(fsum(row) for m in captured for row in m))
        > 3e-11
    ):
        raise ValueError(
            "Coarse tag and extracted branches disagree on their joint mass"
        )
    return [retain[i][2] - fsum(matrix[i][2] for matrix in captured) for i in range(4)]


def causal_capture():
    all_pages = read("extended_decay/wing_capture.json")["settings"]
    ticks = {"162", "164"}
    if set(all_pages["native"]) != ticks:
        raise ValueError("Incomplete capture comparison")
    for tick in sorted(ticks):
        for case, stimulus in (("known", "D_phase"), ("decay", "core_pulse")):
            contrasts = {}
            for setting, pages in all_pages.items():
                panel = pages[tick]
                arms = panel["cases"][case]
                diff = [
                    a - b
                    for a, b in zip(
                        captured_soft(arms[stimulus]),
                        captured_soft(arms["idle"]),
                        strict=True,
                    )
                ]
                # Four subnormalized probabilities enter the capture contrast.
                # The omitted positive stock and 12-place receipt rounding each have a bound.
                discarded = panel["preparation"]["discarded_weight_upper"]
                if not isfinite(discarded) or discarded < 0:
                    raise ValueError("Invalid omitted probability bound")
                allowance = 4 * discarded + 4e-11
                contrasts[setting] = fsum(diff), fsum(diff[1:3]), allowance
            total, joint, error = contrasts["native"]
            cut, _, cut_error = contrasts["exit_cut"]
            if abs(cut) > cut_error:
                raise ValueError("Exit-cut capture contrast exceeds its omission bound")
            if case == "known" and total - error <= 0.01:
                raise ValueError("Known input misses the registered capture criterion")
            if case == "decay" and (total + error >= 0.01 or joint + error >= 0.0):
                raise ValueError("Heavy-source null disposition differs")
            print(
                f"Capture {tick}, {case}: {100 * total:+.5f} points, "
                f"with heavy remnant {100 * joint:+.5f} points"
            )


def separated_energy():
    data = read("extended_decay/larger_transition.json")
    for board, panel in data["boards"].items():
        differences = [
            data["strand_levels"][i] - data["strand_levels"][j] - line
            for (j, i), line in product(
                combinations(range(len(data["strand_levels"])), 2),
                panel["exterior_gaps"],
            )
        ]
        if not differences:
            if panel["closest"] is not None:
                raise ValueError("Energy match without an available exterior line")
            print(f"{board} cells: no eligible separated exterior line")
            continue
        closest = min(differences, key=abs)
        if not isclose(closest, panel["closest"][-1], rel_tol=0, abs_tol=1e-12):
            raise ValueError("Separated spectral mismatch does not reproduce")
        if abs(closest) <= 0.05:
            raise ValueError("Closest line lies within the quoted separation")
        print(
            f"{board} cells: {len(differences)} separated transitions, "
            f"closest detuning {abs(closest):.9f}"
        )


if __name__ == "__main__":
    passive_decay()
    causal_capture()
    separated_energy()
