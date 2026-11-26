import json

import numpy as np

from .measure import ANGLES, FORMS, HERE, SHOTS, decode, run


def test_the_complete_replay_and_held_clock_records():
    data = json.loads((HERE / "records.json").read_text())
    assert data == json.loads(json.dumps(run()))
    for reader in data["readers"]:
        m, d = reader["meter"]
        reserved = {(m, d), tuple(sorted((d, reader["terminal_empty_port"])))}
        assert all(tuple(event["chord"]) not in reserved for event in reader["events"])
        assert sum(event["accepted"] for event in reader["events"]) == 6
        assert reader["events"][0]["before"] == reader["events"][0]["after"]
        scores = []
        for i, case in enumerate(reader["cases"]):
            for form in FORMS:
                banks = case["counts"][form]
                assert np.shape(banks) == (len(ANGLES), 12, 3)
                assert (np.sum(banks, axis=-1) == SHOTS).all()
            target = 24 - 2 * min(6, max(0, i - 1))
            reading = case["decoded"]["native"]
            assert abs(reading["quadratic"] - target) < data["quadratic_radius"]
            assert (
                abs(reading["held_forecast"] - reading["held_frequency"])
                < data["held_difference_radius"]
            )
            scores.append(reading["quadratic"])
        for before, after in zip(scores[1:7], scores[2:8], strict=True):
            assert before - after > 2 * data["quadratic_radius"]
        first = reader["cases"][0]["decoded"]
        assert abs(first["linear_private"]["quadratic"] - 12) < data["quadratic_radius"]
        assert first["exchange_only"]["quadratic"] == 0


def test_decoder_uses_only_bright_events_from_the_receiver():
    data = json.loads((HERE / "records.json").read_text())
    for reader in data["readers"]:
        for case in reader["cases"]:
            bins = np.array(case["counts"]["native"])
            bins[:, :, [1, 2]] = bins[:, :, [2, 1]]
            assert decode(bins, reader["single_B"]) == decode(
                case["counts"]["native"], reader["single_B"]
            )
