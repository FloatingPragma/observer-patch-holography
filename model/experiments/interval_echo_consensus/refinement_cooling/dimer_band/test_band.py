import json
from math import exp, log, sqrt

import numpy as np
from scipy.linalg import expm

from model.executable_model import CELL_SEAMS, federation
from .measure import HERE, STEPS, FLIGHT_TIMES, gates, layouts, outcomes


def test_native_levels_decoder_and_all_initial_strata():
    h, d, f = gates("native")
    modes = np.array([[1, 1], [1, -1], [0, 0]]) / sqrt(2)
    np.testing.assert_allclose(h @ modes, modes * [32, 33], atol=2e-13, rtol=0)
    np.testing.assert_allclose(abs(d @ modes) ** 2, np.eye(3)[:, :2], atol=2e-13)
    np.testing.assert_allclose(
        abs(np.array([expm(-1j * t * h)[1, 0] for t in FLIGHT_TIMES])) ** 2,
        [0, 0.25, 0.5, 1, 0.5],
        atol=2e-13,
    )
    for i in (0, 1):
        for b in STEPS:
            np.testing.assert_allclose(
                outcomes(d, f, i, b),
                [2 ** (-j - 2) for j in range(b)] + [0.5, 0.5 * 2**-b],
                atol=2e-13,
            )
    # Absolute phases can depend on the appointment. The entire conditional
    # density for the stratified launch, including off-diagonals, is Gibbs.
    for b in STEPS:
        no_click = f[:2, :2]
        v = d[:2, :2].conj().T @ np.linalg.matrix_power(no_click, b) @ d[:2, :2]
        density = v @ v.conj().T
        target = modes[:2] @ np.diag([1, 2**-b]) @ modes[:2].T
        np.testing.assert_allclose(
            density / density.trace(), target / target.trace(), atol=2e-13
        )


def test_active_geometry_and_four_to_one_physical_refinement():
    paths = layouts()
    edges = {frozenset(pair) for pair in CELL_SEAMS}
    assert len(paths) == 220
    for blocked, path in paths.items():
        assert len(set(path)) == 4 and set(blocked).isdisjoint(path)
        assert all(frozenset(pair) in edges for pair in zip(path, path[1:]))

    def contacts(board, block):
        return {
            tuple(sorted((p // (12 * block), q // (12 * block))))
            for p, q in board.external
            if p // (12 * block) != q // (12 * block)
        }

    for level in (1, 2):
        fine, coarse = federation(level), federation(level - 1)
        assert fine.cells == 4 * coarse.cells
        assert contacts(fine, 4) == contacts(coarse, 1)
        for beta in (0, log(2), 8):
            p = 1 / (1 + exp(beta))
            for sign in (p, 1 - p):
                assert 4 * sign / fine.cells == sign / coarse.cells


def test_joint_counts_and_separate_one_appointment_calibration():
    data = json.loads((HERE / "records.json").read_text())

    # At most 400 signed comparisons of bounded individual shot indicators.
    # Use two-sided Hoeffding intervals for all three hosts/control words.
    def radius(n):
        return sqrt(log(2 * 400 / 0.01) / (2 * n))

    for board in data["boards"]:
        n = 2 * data["shots_per_rail"] * board["owners"]
        for panel in board["panels"]:
            banks = []
            for row in panel["banks"]:
                by_owner = np.array(row["counts"])
                assert by_owner.shape == (board["owners"], 2, row["steps"] + 2)
                assert (by_owner.sum(axis=2) == data["shots_per_rail"]).all()
                banks.append(by_owner.sum(axis=(0, 1)) / n)
            if panel["mode"] != "native":
                assert banks[-1][-1] - radius(n) > 0.45
                continue
            e = radius(n)
            calibration = 2 * banks[1][-1]
            for b, bank in zip(STEPS, banks, strict=True):
                forecast = [2 ** (-j - 2) for j in range(b)] + [0.5, 0.5 * 2**-b]
                assert max(abs(bank - forecast)) < e
                assert 0.5 * max(0, calibration - 2 * e) ** b <= bank[-1] + e
                assert bank[-1] - e <= 0.5 * min(1, calibration + 2 * e) ** b
            flight = np.asarray(panel["flight_transfers"]) / data["flight_shots"]
            assert max(abs(flight - [0, 0.25, 0.5, 1, 0.5])) < radius(
                data["flight_shots"]
            )
