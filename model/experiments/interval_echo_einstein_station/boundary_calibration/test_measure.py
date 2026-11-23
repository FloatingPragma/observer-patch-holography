"""Verify the pair projection, coherence convention, and held scores from complete counts."""

import json
from types import SimpleNamespace

import numpy as np

from .measure import (
    HELD,
    HERE,
    SHOTS,
    PairSector,
    candidate,
    channel_capacity,
    federation,
    inputs,
    operators,
    primary,
    product_ball,
    spectral_inverse,
)
from .phase_receiver import receiver_blocks


def test_configuration_restriction_keeps_native_surface_diagonal():
    board = federation(0)
    _, launches = inputs(board)
    first, second = product_ball(board, launches)
    kinetic, contrast, keys = operators(board, first, second)
    full = PairSector(board)
    rows = np.searchsorted(full.aa * board.n + full.bb, keys)
    difference = (kinetic + contrast) - full.hamiltonian(candidate.KAPPA)[rows][:, rows]
    assert difference.nnz == 0


def test_partial_trace_preserves_positive_circular_projector_convention():
    # One material port stays outside. Inside, the two source branches occupy distinct ports.
    blocks = receiver_blocks(
        np.eye(2, dtype=complex),
        np.array([0, 0]),
        np.array([12, 13]),
        receiver=1,
        board=SimpleNamespace(n=24),
    )
    mean, imaginary = blocks[1]
    reader = np.eye(12)[0] / np.sqrt(2) + 1j * np.eye(12)[1] / np.sqrt(2)
    np.testing.assert_allclose(reader.conj() @ mean @ reader, 0.5, atol=1e-14)
    np.testing.assert_allclose(reader.conj() @ imaginary @ reader, 0.5, atol=1e-14)
    np.testing.assert_allclose(np.diag(imaginary), 0, atol=1e-14)


def test_archived_original_records_and_scores():
    receipt = json.loads((HERE / "results.json").read_text())
    assert receipt["shots_per_ensemble"] == SHOTS
    for board in receipt["boards"]:
        regional = board["regional"]
        for field in ("regional", "phase_regional"):
            np.testing.assert_equal(
                np.array(board[field]["counts"]).sum(axis=-1), SHOTS
            )
        for rows in board["studies"].values():
            for label in ("reference_counts", "calibration_counts", "held_counts"):
                np.testing.assert_equal(np.array(rows[label]).sum(axis=-1), SHOTS)
            calibration, held = (
                np.array(rows[k]) for k in ("calibration_counts", "held_counts")
            )
            np.testing.assert_equal(held.shape, (len(HELD), 2, 3, 79))
            assert rows["escape"] == primary(held, "escape", regional)
            assert rows["terminal_entropy"] == primary(held, "entropy", regional)
            np.testing.assert_allclose(
                rows["inverse"]["recovered_xyz"],
                spectral_inverse(calibration, held)["recovered_xyz"],
                atol=1e-11,
            )
            np.testing.assert_allclose(
                rows["channel"]["channel_capacity_nats"],
                channel_capacity(calibration)["channel_capacity_nats"],
                atol=1e-13,
            )
            np.testing.assert_equal(
                np.array(rows["phase_challenge"]["held_counts"]).sum(axis=-1), SHOTS
            )
        assert abs(board["phase_regional"]["slopes"][0] + 0.554035237512932) < 0.0001
        assert (
            board["studies"]["registered"]["phase_challenge"]["maximum_full_tv_bound"]
            < 2.2e-13
        )


def test_receiver_selection_and_held_contrast_derive_from_separate_counts():
    receipt = json.loads((HERE / "phase_receiver.json").read_text())
    for board in receipt["boards"]:
        for rows in board["studies"].values():
            calibration = np.array(rows["calibration_counts"])
            np.testing.assert_equal(calibration.sum(axis=-1), SHOTS)
            yes = calibration[:, :, 0] / SHOTS
            uncertainty = np.sqrt((yes * (1 - yes)).sum(axis=0) / SHOTS)
            best = np.argmax(
                np.abs(yes[1] - yes[0]) / np.maximum(uncertainty, 1 / SHOTS)
            )
            assert best == rows["selected_setting"]
            held = np.array(rows["held_counts"])
            np.testing.assert_equal(held.sum(axis=-1), SHOTS)
            assert (held[1, 0] - held[0, 0]) / SHOTS == rows["held_difference"]
            for setting in rows["settings"]:
                assert setting["occupation_tv"] <= setting["trace_distance"] + 1e-15


def test_joint_phase_replication_recomputes_selection_and_contrasts():
    from .joint_replication import REPEATS, contrast

    receipt = json.loads((HERE / "joint_replication.json").read_text())
    assert receipt["shots_per_ensemble"] == REPEATS
    for board in receipt["boards"]:
        for name, study in board["studies"].items():
            cal = np.array(study["calibration_counts"])
            held = np.array(study["held_occupations"])
            routed = np.array(study["held_routed_reader_counts"])
            for tally in (cal, held, routed):
                np.testing.assert_array_equal(tally.sum(axis=-1), REPEATS)
            p = cal / REPEATS
            criterion = 3 * np.sqrt(np.sum(p * (1 - p), axis=0) / REPEATS)
            select = p[1] - p[0] > criterion
            np.testing.assert_array_equal(
                np.flatnonzero(select), study["selected_occupations"]
            )
            assert study["occupation_contrast"] == contrast(held, select)
            assert study["routed_contrast"] == contrast(routed, [0])
            routed_signal = study["routed_contrast"]
            assert routed_signal["difference"] > 5 * routed_signal["standard_error"]
            if name == "coherent":
                signal = study["occupation_contrast"]
                assert signal["difference"] > 35 * signal["standard_error"]
