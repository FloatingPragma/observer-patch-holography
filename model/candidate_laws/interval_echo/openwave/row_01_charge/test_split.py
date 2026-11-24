"""Check deposited counts and the two compact sources independently of the producer."""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.carrier import federation

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "results.json").read_text())


def test_exact_sector_energy_and_enclosing_records():
    for level, row in enumerate(DATA["boards"]):
        board = federation(level)
        for station in row["stations"]:
            lo = station["owner"] * 12
            vector = np.zeros(board.n, dtype=int)
            vector[lo : lo + 12] = station["integer_amplitudes_in_cell"]
            for p in range(board.n):
                assert (
                    board.degree[p] * vector[p]
                    - sum(vector[q] for q in board.neighbours[p])
                    == 6 * vector[p]
                )
                if board.degree[p] != 5:
                    assert vector[p] == 0
            for name, background, amount in (
                ("vacancy_in_host", 1, 0),
                ("single_on_empty", 0, 1),
                ("double_on_empty", 0, 2),
            ):
                # One defect sees five neighbour contrasts and one degree reference.
                native = (
                    8 * (5 * (background - amount) ** 2 + amount**2 - background**2) + 3
                )
                assert (
                    abs(row["free_eigenstates"][name]["energy_above_host"] - native)
                    < 1e-12
                )
            for ring in station["cuts"]:
                edges = ring["cut_oriented_from_inside"]
                endpoints = np.array(edges).ravel()
                assert len(set(endpoints)) == len(endpoints)
                assert not np.any(vector[endpoints])
                assert ring["host_first_second_shared_signed_word"] == [0] * len(edges)
                assert ring["vacancy_at_first_inner_endpoint_signed_word"] == [-1] + [
                    0
                ] * (len(edges) - 1)


def test_daughter_gate_and_changed_law():
    # Every free local term preserves the two-port multiplicities, including dirty memory.
    exchange, contrast = law.local_term_parts()
    for amount in (0, 1, 2):
        counts = np.tile([int(a == amount) + int(b == amount) for a, b in law.SRC], 11)
        for h in (exchange, contrast):
            assert np.max(abs(h * (counts[:, None] - counts[None, :]))) < 1e-13
    first, second = DATA["isolated_symmetric_total_two_challenge"]
    assert first["signed_word_counts"][2] == 0
    # The symmetric pair has energy 32; the singly occupied seam has zero.
    h = np.array([[32.0, 8.0], [8.0, 0.0]])
    expected = 0.5 * abs(expm(-1j * h / 16)[1, 0]) ** 2
    assert abs(expected - second["split_probability"]) < 1e-14
    assert (
        abs(second["signed_word_counts"][2] / DATA["shots_per_setting"] - expected)
        < 0.001
    )


def test_exclusive_vs_independent_raw_counts():
    sign_score = np.array([0, -1, 0, 1, 0])
    radius = (2 * np.log(12 / 0.01) / DATA["precision_shots"]) ** 0.5
    initial_passes = []
    for row in DATA["boards"]:
        setting = row["precision_positive_probes"]
        initial = row["joint_probe_settings"][0]["histories"][-1]
        assert initial["probe_signs"] == setting["probe_signs"] == [1, 1]
        initial_passes.append(
            initial["absolute_score_difference"]
            > 2 * DATA["one_mean_simultaneous_99_percent_error"]
        )
        means = []
        for label in ("shared", "independent_half_mean"):
            counts = np.array(setting[label]["counts"])
            assert counts.sum() == DATA["precision_shots"]
            means.append(sign_score @ counts @ sign_score / counts.sum())
        assert means[1] - means[0] > 2 * radius
        assert abs((means[1] - means[0]) - setting["absolute_score_difference"]) < 1e-15
        pilot, removed = row["joint_probe_settings"]
        assert (
            removed["histories"][-1]["joint_total_variation"]
            < pilot["histories"][-1]["joint_total_variation"] / 100
        )
    assert initial_passes == [True, False, False]
