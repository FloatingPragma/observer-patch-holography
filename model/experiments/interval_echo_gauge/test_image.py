import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import currents
from model.executable_model import cell
from .__main__ import ordered_word, portable, run


def test_kernel_measurement_and_controls():
    result = portable(run())
    assert result == json.loads(Path(__file__).with_name("result.json").read_text())
    coarse, fine = result["convergence"].values()
    assert len(fine["resolved_kernel"]) == 6
    assert [0, 1, 0] in fine["resolved_kernel"] and [1, 1, 1] not in fine["resolved_kernel"]
    assert fine["maximum_grid_operator_error"] < coarse["maximum_grid_operator_error"] / 15
    for name, row in result["terminal_measurements"].items():
        assert sum(row["counts"]) == result["shots_per_arm"]
        assert abs(sum(row["probabilities"]) - 1) < 2e-9
        assert row["resolved_increase"] == (name == "diagonal_comparison")


def test_full_public_port_pulses_on_all_pair_states():
    # Compare a whole symmetric slice with independent public pulses in the original 66-state basis.
    receipt = json.loads(Path(__file__).with_name("result.json").read_text())
    _, _, _, _, bands, _, _, _, addresses, _ = currents._service(cell(), True, 1.)
    basis = np.column_stack(bands)
    coefficients = receipt["axis_port_coefficients"][0]
    small_angle = 2 * np.pi / 512
    expected = ordered_word(addresses, coefficients, small_angle, steps=1)
    whole_expected = np.eye(66) + basis @ (expected - np.eye(6)) @ basis.T
    half = [currents.pulse(p, -0.5 * small_angle * x) for p, x in enumerate(coefficients)]
    whole = np.eye(66)
    for unitary in half + half[::-1]:
        whole = unitary @ whole
    assert np.linalg.norm(whole - whole_expected, 2) < 1e-12
