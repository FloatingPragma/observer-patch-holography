import json

import numpy as np

from model.experiments.interval_echo_gauge.arrival_ladder import measure as ladder
from model.experiments.interval_echo_angular.charged_orientation.flux_doublet import (
    measure as spin,
)
from . import measure as joint


def test_full_receiver_calibration_and_reference_transfer():
    owner, ports = joint.stations()
    columns, projector, (r, s) = spin.charged_source(owner)
    gate = ladder.interaction(np.pi, ladder.K)
    local = []
    for load in (0, 1):
        pair = [ladder.WORDS.index((load, a, b)) for a, b in ((1, 2), (2, 1))]
        local.append(gate[np.ix_(pair, pair)])
    # Repeated one-address gates have a common phase from the eleven
    # filled controls. Retain it when verifying the full amplitude.
    expected = []
    for vacancy in range(len(ports)):
        rail_gate = np.eye(2, dtype=complex)
        for addressed in range(len(ports)):
            rail_gate = local[int(vacancy != addressed)] @ rail_gate
        expected.append(rail_gate[1, 0])
    aperture = ladder.aperture_map(ports)
    k = aperture["charged_to_neutral"]
    np.testing.assert_allclose(k, np.diag(expected), atol=3e-13)
    np.testing.assert_allclose(
        columns.conj().T @ projector @ columns, np.eye(2), atol=1e-12
    )
    np.testing.assert_allclose(np.diag(projector), np.ones(12) / 6, atol=1e-13)
    # Independent test with a retained reference. Perform the neutral turn
    # after conversion and trace its stored address in the which-port control.
    bell = np.reshape(columns.T, (24,)) / np.sqrt(2)
    rotated = np.kron(np.eye(2), r @ k) @ bell
    desired = np.kron(np.eye(2), r) @ bell
    np.testing.assert_allclose(abs(np.vdot(desired, rotated)) ** 2, 1, atol=1e-12)
    weights = np.sum(abs(bell.reshape(2, 12)) ** 2, axis=0)
    np.testing.assert_allclose(np.sum(weights**2), 1 / 12)
    np.testing.assert_allclose((r @ s + s @ r) @ k @ columns, 0, atol=1e-12)


def test_tomography_and_dark_settings_against_archive():
    actual = joint.measure()
    saved = json.loads((joint.HERE / "results.json").read_text())
    bound = actual["simultaneous_99_process_entry_radius"]
    for name, setting in actual["settings"].items():
        for acquisition in setting["records"]:
            for read in acquisition["reads"]:
                assert sum(read["terminal_counts"]) == actual["shots_per_read"]
                assert len(read["terminal_counts"]) == 24
        assert setting["records"] == saved["settings"][name]["records"]
        expected = np.array(setting["predicted_pauli_transfer"])
        estimated = np.array(setting["sampled_pauli_transfer"])
        np.testing.assert_allclose(estimated, expected, atol=bound)
    changed = actual["zero_contrast_receiver"]
    np.testing.assert_allclose(changed["predicted_pauli_transfer"], 0, atol=1e-12)
    for launch in changed["records"]:
        for bank in launch["reads"]:
            assert sum(bank["terminal_counts"][:12]) == 0
            assert sum(bank["terminal_counts"][12:]) == actual["shots_per_read"]
    settings = actual["settings"]
    np.testing.assert_allclose(
        settings["coherent"]["predicted_pauli_transfer"], np.eye(4), atol=1e-12
    )
    np.testing.assert_allclose(
        settings["address_recorded"]["predicted_pauli_transfer"],
        np.diag([1 / 6, 1 / 18, 1 / 18, 1 / 18]),
        atol=1e-12,
    )
    np.testing.assert_allclose(
        settings["no_drive"]["predicted_pauli_transfer"], 0, atol=1e-12
    )
    np.testing.assert_allclose(
        settings["neutral_half_turn"]["predicted_pauli_transfer"],
        actual["calibrated_turn_matrix"],
        atol=1e-11,
    )
    np.testing.assert_allclose(actual["family_neutral_success"], [1, 1, 1], atol=1e-12)


def test_energy_herald_and_both_charge_banks():
    from model.experiments.interval_echo_lepton_register.pole_meter.measure import setup

    owner, ports = joint.stations()
    frame, _, turns = spin.charged_source(owner)
    _, numerator, denominator, *_ = setup()
    projection = numerator[np.ix_(ports, ports)] / denominator
    weights, rotation = np.linalg.eigh(frame.conj().T @ projection @ frame)
    source = projection @ frame @ rotation
    saved = json.loads((joint.HERE / "results.json").read_text())
    data = saved["energy_selected_round_trip"]
    np.testing.assert_allclose(data["initial_compact_admission"], weights)
    assert 0 < min(weights) < max(weights) < 0.05
    for bank in data["banks"]:
        i = bank["initial_principal_mode"]
        op = dict(park=np.eye(12), first=turns[0], second=turns[1])[
            bank["neutral_turn"]
        ]
        back = np.linalg.norm(projection @ op @ source[:, i]) ** 2
        expected = [1 - weights[i], 0, 0, back, weights[i] - back]
        np.testing.assert_allclose(bank["original_probabilities"], expected, atol=2e-12)
        assert sum(bank["counts"]) == saved["shots_per_read"]
        np.testing.assert_allclose(
            np.array(bank["counts"]) / saved["shots_per_read"], expected, atol=0.001
        )
