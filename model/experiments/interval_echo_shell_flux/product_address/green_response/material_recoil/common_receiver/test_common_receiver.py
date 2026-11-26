import json
from itertools import product

import numpy as np

from . import measure, receiver


def test_all_native_swaps_and_held_echoes():
    for stiffness in (4, 8):
        assert receiver.swap_error(stiffness) < 1e-12
        for probe in (1, 2):
            energies, flight = receiver.native_terms(probe, stiffness)
            np.testing.assert_array_equal(np.diff(flight.diagonal()), [-0.5])
            for sense, ref in product(range(3), (1, 2)):
                phases = [
                    receiver.echo(probe, stiffness, sense, ref, t, False)
                    for t in range(3)
                ]
                observed = sum(v[1] * v[0].conj() for v in phases) / 3
                expected = float(stiffness == 4 or sense == ref)
                np.testing.assert_allclose(observed, expected, atol=2e-13, rtol=0)
            # The two pulse windows cancel every carrier-dependent onsite shift.
            diffs = energies[:, 0] - energies[:, 1]
            for s in (1, 2):
                assert 2 * (diffs[0] - diffs[s]) == 4 * stiffness * probe * s + 1 + int(
                    probe == s
                )


def test_unit_gain_all_source_words_and_both_receiver_signs():
    for source in product((1, 2), repeat=3):
        for probe, blind in product((1, 2), (False, True)):
            for word in receiver.WORDS:
                for sign in (-1, 1):
                    gain = sign * 0.4 * all(word) / (8 if blind else 1)
                    p = receiver.right_probability(
                        probe,
                        8,
                        word,
                        source,
                        blind,
                        sign,
                        receiver.RELEASE,
                        False,
                        False,
                    )
                    assert abs(p - 0.5 - gain) < 3e-13


def test_archived_records_against_the_independent_green_flux():
    data = json.loads((measure.HERE / "results.json").read_text())
    u = measure.visits_before_drain()
    setups = list(measure.settings())
    assert len(setups) == len(data["panels"]) == 32
    observed = []
    for setting, panel in zip(setups, data["panels"], strict=True):
        assert panel["setting"] == json.loads(json.dumps(measure.asdict(setting)))
        arms = panel["arms"]
        total = sum(sum(a["starts"]) for a in arms)
        assert total == data["appointments_per_panel"]
        assert all(
            0 <= n <= m
            for a in arms
            for n, m in zip(a["right"], a["starts"], strict=True)
        )
        shift = (sum(sum(a["right"]) for a in arms) - total / 2) / data[
            "source_releases"
        ]
        exact = measure.signed_readout(setting)
        error = measure.bound(exact, u, len(setups))
        assert panel["mean_per_release"] == shift
        assert abs(float(np.sum(u * exact)) - panel["predicted"]) < 1e-12
        assert abs(error - panel["simultaneous99_half_width"]) < 1e-12
        assert abs(shift - panel["predicted"]) < error
        if setting.name == "matched":
            observed.append((setting.face[0], shift, error))
    assert all(s - e > 0 if x == 4 else s + e < 0 for x, s, e in observed)
    assert 0 < data["last_drain"] <= measure.PLAN["maximum_appointments"]


def test_ordinary_response_and_independent_dimer_clock():
    record = measure.calibrations()
    assert (
        record == json.loads((measure.HERE / "results.json").read_text())["calibration"]
    )
    counts, shots = record["terminal_counts"], record["starts_per_bank"]
    error = np.sqrt(np.log(2 * len(counts) / 0.01) / (2 * shots))
    for key, p in [
        ("ordinary full match", 0.9),
        ("ordinary single miss", 0.5),
        ("trit subpulse", 0.25),
        ("full release", np.sin(receiver.RELEASE / 2) ** 2),
        ("half release", np.sin(receiver.RELEASE / 4) ** 2),
    ]:
        assert abs(counts[key] / shots - p) < error
