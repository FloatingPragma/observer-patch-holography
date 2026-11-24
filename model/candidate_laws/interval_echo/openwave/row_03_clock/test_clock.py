"""Independent sector, parity-receiver and count checks for the vacancy clock."""

import json
from pathlib import Path

import numpy as np

from model.executable_model import cell, federation

from .measure import SHOTS, spectral_bounds

RECORD = Path(__file__).with_name("results.json")


def test_whole_board_energy_enclosure_from_incidence():
    for level in range(3):
        board = federation(level)
        # Construct the row-sum intervals directly from degrees and the signed defect rule.
        for sign in (-1, 1):
            potential = 48 + 16 * sign * (6 - board.degree)
            direct = [float(np.min(potential)), float(np.max(potential + board.degree))]
            assert spectral_bounds(board, sign)["energy_interval"] == direct
    # Independent diagonalization includes every gateway link of the 20-cell board.
    board = federation(0)
    laplace = np.diag(board.degree) - board.adjacency().toarray()
    levels = np.linalg.eigvalsh(np.diag(48 - 16 * (6 - board.degree)) + laplace / 2)
    assert levels[-1] - levels[0] <= 22
    assert np.min(abs(levels - 35)) < 1e-10


def test_sharp_source_invariant_subspace_and_cosine_weight():
    a = cell().adjacency().toarray().astype(int)
    u = np.array([0, 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0])
    v = a @ u
    np.testing.assert_array_equal(a @ v, 5 * u)
    assert [u @ u, u @ v, v @ v] == [4, 4, 20]
    for t in (0, 0.2, 0.7):
        # For onsite return each projected energy has a positive real weight.
        weights = [
            float(x[4] ** 2 / (x @ x)) for x in (np.sqrt(5) * u - v, np.sqrt(5) * u + v)
        ]
        p = abs(weights[0] + weights[1] * np.exp(1j * t * np.sqrt(5))) ** 2
        np.testing.assert_allclose(
            p, np.sum(np.square(weights)) + 0.025 * np.cos(t * np.sqrt(5))
        )
    for panel in json.loads(RECORD.read_text())["unmonitored_sharp_returns"]:
        np.testing.assert_allclose(panel["return_cosine_amplitude_lower_bound"], 1 / 40)
        np.testing.assert_array_equal(
            np.array(panel["return_counts"]) + panel["other_counts"], SHOTS
        )


def test_monitored_clock_and_energy_controls_from_count_record():
    panels = json.loads(RECORD.read_text())["panels"]
    for name, panel in panels.items():
        rows = np.array(panel["measurements"])
        np.testing.assert_array_equal(rows[:, 1:].sum(axis=1), SHOTS)
        assert np.all((0 <= rows) & (rows <= SHOTS))
        times = 16 * np.arcsin(np.sqrt(rows[:, 0] / SHOTS))
        prob = rows[:, 1:3].sum(axis=1) / SHOTS
        q = panel["energy_read"][
            "counts_antisymmetric_symmetric_elsewhere_source_failed"
        ]
        assert sum(q) == SHOTS
        kinetic = 30 * panel["exchange_rate"] * q[0] / sum(q[:3])
        np.testing.assert_allclose(kinetic, panel["energy_read"]["exchange_estimate"])
        f = panel["frequency_fit"]
        if name in ("no_antipodal_gateways", "exchange_off_at_release"):
            assert f["frequency_per_clock_gap_if_resolved"] is None
            assert np.ptp(prob) < 0.008
        else:
            # Predict held data at the analytically fixed frequency regardless of the fitted peak.
            design = np.column_stack(
                (
                    np.ones(len(times)),
                    np.cos(np.sqrt(5) * times),
                    np.sin(np.sqrt(5) * times),
                )
            )
            coeff = np.linalg.lstsq(design[:65], prob[:65], rcond=None)[0]
            assert np.sqrt(np.mean((design[65:] @ coeff - prob[65:]) ** 2)) < 0.006
            assert 2 * np.linalg.norm(coeff[1:]) > 0.03
            assert abs(
                f["frequency_per_clock_gap_if_resolved"] - np.sqrt(5)
            ) < 0.01 * np.sqrt(5)
    energy_difference = (
        panels["antipodal_80"]["energy_read"]["estimated_host_referenced_mean_energy"]
        - panels["half_contrast"]["energy_read"][
            "estimated_host_referenced_mean_energy"
        ]
    )
    assert abs(energy_difference - 16) < 0.12
