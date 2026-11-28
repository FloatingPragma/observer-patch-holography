"""Offline checks of the four discriminating controls in the reviewer drafts."""

from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from math import pi, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
STUDIES = ROOT / "model/candidate_laws/interval_echo/openwave"
SOURCE = ROOT / "openwave/interval_echo/row_19/selected_collar"


def load(path):
    return json.loads(path.read_text())


def same_charge():
    energy = load(SOURCE / "parked_energy_charge.json")
    pair = load(SOURCE / "matched_charge_transport.json")
    flight = load(SOURCE / "parked_pair_flight.json")
    # Read both the reference costs and original shared-source weights from
    # the published joint meter tickets. The change vectors use calibration order.
    order = energy["meter_order"]
    tickets = energy["counts"]["shared"]
    ref = energy["counts"]["reference"]
    total = pair["denominator"]
    delta_by_price = Counter()
    flight_shift = 0.0
    return_prob = flight["unit_time_joint_return"]["native"]["by_source_word_pair"]
    assert sorted(order) == list(range(9)) and sum(tickets) > 0
    u, v = pair["u"], pair["v"]
    assert sum(u) == sum(v) == 0
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            d = u[a] * v[b] + v[a] * u[b]
            assert tickets[9 * i + j] * total % sum(tickets) == 0
            assert tickets[9 * i + j] * total // sum(tickets) + d > 0
            delta_by_price[ref[9 * i + j]] += d
            flight_shift += d * return_prob[9 * a + b] / total
    assert not any(delta_by_price.values())
    lo, hi = pair["unit_time_free_pair_return_change"]["native"]
    assert lo <= flight_shift <= hi < 0
    return flight_shift


def current_loop():
    data = load(STUDIES / "row_26_newton/current_loop.json")
    mass, radius, gain, tick = (
        Fraction(data[k]) for k in ("shell_mass", "shell_radius", "G", "tick")
    )
    spin = Fraction(2, 3) * mass * radius**2 * Fraction(data["rotation_rates"][2])
    assert spin == Fraction(1, 384)
    assert [r["reference_radius"] for r in data["rings"]] == [4, 8, 16]
    readings = []
    for ring in data["rings"]:
        r = ring["reference_radius"]
        tickets = ring["returns_in_receiver_ticks"]
        for setting in ("-1", "0", "1"):
            cw, ccw = tickets["scalar_only"][setting]
            assert cw == ccw
        pos = tickets["gem"]["1"]
        neg = tickets["gem"]["-1"]
        assert pos == neg[::-1] and tickets["gem"]["0"][0] == tickets["gem"]["0"][1]
        odd = (pos[0] - pos[1]) * tick
        expected = -8 * pi * float(gain * spin) / r / sqrt(1 - 2 * gain * mass / r)
        assert abs(float(odd) - expected) <= float(tick)
        readings.append(float(odd))
    return readings


def common_delay():
    data = load(STUDIES / "row_27_metric/parked_signals/common_delay.json")
    assert set(data["tests"]) == {"plus", "minus", "fourth_order", "shared"}
    for row in data["tests"].values():
        assert row["first_changed_time_power"] == 14
        slopes = list(map(Fraction, row["blank_slopes_in_squared_time"]))
        change = list(map(Fraction, row["first_changes"]))
        assert slopes == [Fraction(-5, 4), Fraction(1, 4), Fraction(-5, 4)]
        required = [c / s for c, s in zip(change, slopes, strict=True)]
        assert len(set(required)) > 1
    return len(data["tests"])


def writable_messages():
    data = load(STUDIES / "row_28_acceleration/cut_messages.json")
    seams = data["cut_seams"]
    assert len(set(seams)) == len(seams) == 25
    masks = [
        sum(((seed % 2) ^ (((seed // 2) & i).bit_count() % 2)) << i for i in range(25))
        for seed in range(64)
    ]
    refs = [seams.index(site) for site in data["affine_reference_seams"]]

    def pattern(word, addresses):
        return tuple((word >> i) & 1 for i in addresses)

    assert len(set(pattern(p, refs) for p in masks)) == len(masks) == 64
    # These masks are a six-dimensional subspace. Each point of the
    # quotient is a message; every coset fills each triple identically.
    for triple in combinations(range(25), 3):
        census = Counter(pattern(p, triple) for p in masks)
        assert len(census) == 8 and set(census.values()) == {8}
    messages = 2 ** len(seams) // len(masks)
    assert messages == data["addressable_messages"] == 524288
    return messages


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without -O so that the scientific checks execute.")
    print(
        json.dumps(
            dict(
                equal_charge_return_shift=same_charge(),
                rotating_loop_odd=current_loop(),
                incompatible_common_delays=common_delay(),
                distinct_addressable_messages=writable_messages(),
            ),
            indent=2,
        )
    )
