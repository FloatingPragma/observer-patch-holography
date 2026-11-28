"""Rebuild the integrated release without the producer or its matrix generator."""

from fractions import Fraction

import numpy as np
from scipy.linalg import expm

from . import raw_verify as raw
from .test_supplements import rising_bank


def test_original_words_close_the_launch_premise():
    receipt = raw.read_json(raw.HERE / "donor_release.json")
    clocks = raw.read_json(raw.HERE / "clock.json")["runs"]
    archive = raw.read_json(raw.DONOR / "analyzers.json")
    tree = raw.component(0)
    chosen = dict(raw.spectroscopy())
    finished = float((1 - Fraction(1, 2**96)) ** (len(raw.board(0)[0]) + 1))
    observed_ports = receipt["physical_port_order"]
    exposures = receipt["terminal_exposures"]
    maximum = 0.0
    for row, clock in zip(receipt["launches"], clocks, strict=True):
        strength = clock["strength_in_simulator"]
        assert row["strength_in_simulator"] == strength
        reading = rising_bank(clock["counts"])
        if strength == 0:
            assert reading is None and row["disposition"] == "inconclusive_clock"
            continue
        ratio = strength / reading[0]
        pulse = expm(-4j * ratio * archive["step"] * tree.h2)
        migrations = {t: raw.propagation(0, strength, t) for t in exposures}
        assert {run["column"] for run in row["launches"]} == chosen.keys()
        for run in row["launches"]:
            column = run["column"]
            a, b = map(tree.words.index, chosen[column])
            vector = np.eye(len(tree.h2), dtype=complex)[0]
            for angles in archive["pulses"][column]["angles"]:
                vector = pulse @ vector
                for port, angle in zip(archive["ports"], angles, strict=True):
                    vector[tree.sites == port] *= np.exp(-1j * ratio * angle)
            ray, _ = raw.joint_attempt(-1, (1, 1))
            read = np.zeros(len(vector))
            read[[a, b]] = [float(ray[1]), float(ray[2])]
            np.testing.assert_array_equal(tree.h2 @ read, read)
            click = finished * abs(read @ vector) ** 2
            assert abs(click - run["launch_witness_probabilities"][-2]) < 2e-12
            assert abs(8 * click - run["computed_all_time_floor"]) < 2e-10
            assert sum(run["simulated_witness_counts"]) == receipt["shots_per_witness"]
            p = abs(vector) ** 2
            for record, t in zip(run["terminal"], exposures, strict=True):
                assert record["exchange_exposure"] == t
                flight, swaps = migrations[t]
                for name, probability in (
                    ("coherent", abs(flight @ vector) ** 2),
                    ("dephased", abs(flight) ** 2 @ p),
                    ("rate_matched_swaps", swaps @ p),
                ):
                    expected = [
                        finished * probability[tree.sites == site].sum()
                        for site in observed_ports
                    ]
                    gap = max(
                        abs(np.subtract(record["port_probabilities"][name], expected))
                    )
                    maximum = max(maximum, gap)
    print(
        f"32 preparations, 5 free waits and 3 arms, max port-bin discrepancy {maximum:.3g}"
    )
    assert maximum < 2e-10
