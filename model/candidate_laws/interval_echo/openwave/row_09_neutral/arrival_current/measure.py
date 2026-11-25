"""Continue the archived downstream signed arrivals through a local contact receiver."""

import json
from functools import cache
from pathlib import Path
import numpy as np
from .apparatus import SENSORS, certify, install, standards, transp
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.measure import (
    gate,
    WORDS,
    ladder,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release import (
    measure as release,
)

HERE = Path(__file__).resolve().parent
OPS = ("live", "pulse_off", "blocked_donor", "half_contrast")


@cache
def record(sensor, guard, mode, word):
    # The setting is fixed before the input. The trigger sees only its own
    # signed contrast. Pulse_off, overloaded donor and changed law keep the input.
    u, v = SENSORS[sensor]
    donor = guard if mode == "blocked_donor" else 1
    x, a, g, d = install_once(sensor, guard, donor)["workstation"]
    flag = install_once(sensor, guard, donor)["persistent_flag"]
    diff = certify()[word]
    if abs(diff) != 2:
        return []
    # Sorted-address sign of the written head. Copy this classical read by
    # an addressed calibrated blank exchange of the orthogonal pointer rails.
    fixture = install_once(sensor, guard, donor)
    board = dict.fromkeys(range(ladder.federation(0).n), 1)
    board.update(dict(fixture["mounted_nonunit"]))
    board[u], board[v] = word
    if diff < 0:
        board[flag[0]], board[flag[1]] = board[flag[1]], board[flag[0]]
    # At this point the input echo has been undone. Transfer exactly these
    # endpoints through the opened hatches, without supplying fresh defects.
    for src, dst in [(u, x), (v, a)]:
        board[src], board[dst] = board[dst], board[src]
    assert board[u] == board[v] == 1
    angle = 0.0 if mode == "pulse_off" else np.pi
    k = 4 if mode == "half_contrast" else 8
    mat = gate(angle, guard, k)
    inp = tuple(board[t] for t in (x, a, d))
    dist = abs(mat[:, WORDS.index(inp)]) ** 2
    outcomes = []
    for idx in np.flatnonzero(dist > 1e-12):
        y = board.copy()
        y[x], y[a], y[d] = WORDS[idx]
        expected = y.copy()
        low, ref = fixture["charge_reader"]
        observed = []
        for path in fixture["read_routes"]:
            for i, j in transp(path):
                y[i], y[j] = y[j], y[i]
            assert y[ref] == 0 and path[-1] == low < ref
            observed.append(calibrated_loads()[y[low]])
            for i, j in transp(path):
                y[i], y[j] = y[j], y[i]
            assert y == expected
        assert observed == [y[z] for z in [x, a, d, *flag]]
        for src, dst in [(u, x), (v, a)]:
            y[src], y[dst] = y[dst], y[src]
        assert y[x] == y[a] == 1 and sum(y.values()) == sum(board.values())
        assert tuple(y[f] for f in flag) == ((2, 0) if diff > 0 else (0, 2))
        q = y[u] + y[v] - 2
        assert q + y[d] - donor == 0
        outcomes.append(
            dict(
                outgoing=[y[u], y[v], y[d]],
                regional_charge=q,
                donor_change=y[d] - donor,
                retained=[y[f] for f in flag],
                probability=float(dist[idx]),
            )
        )
    assert abs(sum(row["probability"] for row in outcomes) - 1) < 5e-13
    return outcomes


install_once = cache(install)
calibrated_loads = cache(ladder.read_loads)


def run():
    # Build the complete response and sign standards before opening the free-bank file.
    certify()
    assert ladder.read_loads() == [0, 1, 2]
    sites = {
        s: dict(**install_once(s, 2, 1), independent_standards=standards(s))
        for s in SENSORS
    }
    cal = []
    for sensor in SENSORS:
        for guard in (0, 2):
            for mode in OPS:
                for word in ((0, 2), (2, 0)):
                    cal.append(
                        dict(
                            sensor=sensor,
                            guard=guard,
                            mode=mode,
                            word=word,
                            outputs=record(sensor, guard, mode, word),
                        )
                    )
    original = json.loads((release.HERE / "results.json").read_text())
    assert original["count_order"] == [
        "00",
        "01",
        "02",
        "10",
        "11",
        "12",
        "20",
        "21",
        "22",
        "no far-hub charge",
        "original source refused",
    ]
    rows = []
    summary = {}
    for held in original["held"]:
        if held["starts"] != original["high_statistics_starts"]:
            continue
        for label in ("free", "two_gates_closed"):
            raw = held["data"][label]["counts"]
            assert sum(raw) == held["starts"]
            echo = [0] * 5
            for i, count in enumerate(raw[:9]):
                echo[certify()[divmod(i, 3)] + 2] += count
            rb = dict(
                bench=held["bench"],
                sensor=held["sensor"],
                sign=held["sign"],
                winding=held["k"],
                time=held["t"],
                release_control=label,
                input_counts=raw,
                contrast=echo,
                matched={},
            )
            for guard in (0, 2):
                for mode in OPS[:3]:
                    signature = []
                    joint = []
                    for word in ((0, 2), (2, 0)):
                        (out,) = record(held["sensor"], guard, mode, word)
                        assert out["probability"] > 1 - 5e-13
                        joint.append(
                            [
                                *out["retained"],
                                *out["outgoing"],
                                raw[3 * word[0] + word[1]],
                            ]
                        )
                        signature.append(
                            [
                                out["regional_charge"],
                                out["donor_change"],
                                raw[3 * word[0] + word[1]],
                            ]
                        )
                    q = [
                        sum(count for charge, _, count in signature if charge == s)
                        for s in (-1, 0, 1)
                    ]
                    assert sum(q) == echo[0] + echo[4]
                    # The full input table retains every no-trigger, miss and refusal.
                    rb["matched"][f"{guard}_{mode}"] = dict(
                        charges=q, flag_and_outgoing=joint
                    )
                    key = f"{label}/{guard}/{mode}"
                    acc = summary.setdefault(key, [0] * 6)
                    for i, v in enumerate([*q, sum(raw[:9]) - sum(q), *raw[9:]]):
                        acc[i] += v
            rows.append(rb)
    return dict(
        input_order=original["count_order"],
        source="neutral_ladder/release/results.json, independent held high-statistics panel",
        joint_columns=[
            "flag rail 1",
            "flag rail 2",
            "first receiving port",
            "second receiving port",
            "donor",
            "count",
        ],
        summary_order=[
            "negative daughter",
            "unchanged accepted pair",
            "positive daughter",
            "other downstream words",
            "no source-hub arrival",
            "original refusal",
        ],
        summary=summary,
        sites=sites,
        calibration=cal,
        held=rows,
    )


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    result = run()
    if args.write:
        (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    else:
        print(json.dumps(result["summary"], indent=2))
