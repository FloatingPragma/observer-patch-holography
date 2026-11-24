"""Cross-check translation reduction and decode held contraction records independently."""

import json
from itertools import combinations

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from . import measure


def test_two_body_bloch_step_matches_full_pair_evolution():
    # Enumerate every position pair. Build both layers directly in the
    # distinguishable port basis, projecting to single occupancy at each port.
    n = 5
    pairs = list(combinations(range(2 * n), 2))
    index = {pair: i for i, pair in enumerate(pairs)}
    law = candidate.native_law
    a, b = law.local_term_parts()
    exact = np.eye(len(pairs), dtype=complex)
    for parity, t in enumerate(measure.THETA):
        gate = expm(-1j * t * (a + 8 * b))[:9, :9]
        for x in range(parity, 2 * n, 2):
            y = (x + 1) % (2 * n)
            step = np.zeros_like(exact)
            for j, pair in enumerate(pairs):
                cols = (2 * int(x in pair), 2 * int(y in pair))
                for target, amp in zip(law.SRC, gate[:, law.SRC.index(cols)]):
                    if abs(amp) < 1e-13:
                        continue
                    assert all(v in (0, 2) for v in target)
                    new = set(pair) - {x, y}
                    new.update(z for z, v in zip((x, y), target) if v)
                    step[index[tuple(sorted(new))], j] += amp
            exact = step @ exact
    entries = measure.transfers(n, 8)
    for m in range(n):
        lift = np.zeros((len(pairs), 2 * n - 1), complex)
        for col in range(2 * n - 1):
            x = col % 2
            y = x + 1 + col // 2
            for h in range(n):
                xy = tuple(sorted(((x + 2 * h) % (2 * n), (y + 2 * h) % (2 * n))))
                lift[index[xy], col] = np.exp(2j * np.pi * h * m / n) / np.sqrt(n)
        np.testing.assert_allclose(
            exact @ lift, lift @ measure.matrix(n, entries, m), atol=7e-13
        )


def test_length_and_current_from_raw_counts():
    data = json.loads((measure.HERE / "results.json").read_text())
    for study in data["studies"]:
        length = {}
        for panel in study["panels"]:
            counts = np.array(panel["separation"]["0"]["counts"], dtype=np.int64)
            assert counts.sum() == measure.SHOTS
            d = np.arange(1, len(counts))
            length[panel["winding"]] = sum(d * counts[:-1]) / counts[:-1].sum()
            means = []
            for row in panel["odd_occupancy_counts"]:
                assert sum(row) == measure.SHOTS
                means.append((row[1] + 2 * row[2]) / sum(row[:3]))
            speed = (2 * means[1] - means[0] - means[2]) / 4
            assert abs(speed - panel["velocity_cells_per_period"]) < 1e-12
        ranges = []
        for study_row in study["contraction"]:
            m = study_row["winding"]
            ratio = length[m] / length[0]
            assert abs(ratio - study_row["length_ratio"]) < 1e-13
            assert (
                study_row["length_ratio_interval"][1]
                < study_row["prediction_interval"][0]
            )
            if m > 0:
                ranges.append(study_row["compatible_limiting_speed"])
        assert max(b[0] for b in ranges) > min(b[1] for b in ranges)


def test_frame_records_and_independent_signal():
    data = json.loads((measure.HERE / "results.json").read_text())
    for study in data["studies"]:
        for row in study["front_reads"]:
            front, beyond, inside, failed = row["front_beyond_inside_failed"]
            assert front + inside == measure.SHOTS and not (failed or beyond)
            p = np.cos(measure.DETUNING) ** (2 * row["period"])
            assert abs(front / measure.SHOTS - p) < 0.0001
        for obs in study["second_orientation_center"].values():
            assert obs["selected_line_fidelity"] > 0.99998
            assert (
                abs(obs["speed"] - obs["ideal_phase_secant_speed"])
                < obs["speed_error_including_filter"]
            )
            assert min(o["visibility"] for o in obs["appointments"]) > 0.4
        for obs in study["first_orientation_center"].values():
            assert (
                min(sum(a["center_counts"][:-1]) for a in obs["appointments"])
                / measure.SHOTS
                < 0.00003
            )
            assert obs["speed_error_including_filter"] is None
