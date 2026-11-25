import json
from math import pi

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from . import measure as run
from .. import measure as pole


def exchanged(state, route, tau):
    swaps = list(zip(route[:-2], route[1:-1]))
    for a, b in swaps:
        state[[a, b]] = state[[b, a]]
    state[route[-2:]] = expm(-1j * tau * run.rails.seam_block(1)) @ state[route[-2:]]
    for a, b in reversed(swaps):
        state[[a, b]] = state[[b, a]]


def test_owner_specific_full_board_ramsey_paths_and_free_projected_qubit():
    owner, ports, board, num, den, c, _, inputs, _, u, read, ref, route = run.encoding()
    h, e, _, _, _ = pole.spectrum(1)
    assert owner == 6 and pole.storage_record()["home_cells"][1] == 13
    np.testing.assert_allclose(read[:2].conj().T @ read[:2], num / den, atol=2e-12)
    code = np.zeros((board.n, 2), complex)
    code[ports] = u
    np.testing.assert_allclose(h @ code, 35 * code, atol=6e-12)
    # A coherent pair in owner 6 is prepared separately from the reference in 13.
    ref_cell = np.flatnonzero(
        np.linalg.norm(ref.reshape(-1, 12), axis=1) > 1e-10
    ).item()
    remote = slice(12 * ref_cell, 12 * ref_cell + 12)
    ref_word = pole.compact_preparation(ref[remote].real, route[-1] % 12)
    for k in (0, 6, 7):
        pulse = run.transfer.preparation(inputs[:, k], owner)
        state = np.zeros(board.n, complex)
        state[route[0]] = 1
        exchanged(state, route, pi / 2)
        state[ports] = pulse @ state[ports]
        state[remote] = ref_word @ state[remote]
        for time in (0, 0.375):
            trial = expm_multiply(-1j * time * h, state)
            trial[ports] = pulse.conj().T @ trial[ports]
            trial[remote] = ref_word.conj().T @ trial[remote]
            exchanged(trial, route, 3 * pi / 2)
            probabilities = abs(trial) ** 2
            source = np.zeros(board.n, complex)
            source[ports] = inputs[:, k]
            z = np.vdot(source, expm_multiply(-1j * time * h, source)) * np.exp(
                1j * time * e
            )
            np.testing.assert_allclose(
                probabilities[[route[0], route[-1]]], pole.outcomes(z)[:2], atol=2e-12
            )


def test_counts_principal_weight_certificate_and_all_turns():
    record = json.loads((run.HERE / "results.json").read_text())
    fresh = run.measure()
    owner, _, _, num, den, c, _, _, _, _, _, _, _ = run.encoding()
    p = num / den
    for name in record["ramsey_banks"]:
        r = record["ramsey_banks"][name]
        v = fresh["ramsey_banks"][name]
        assert sum(r["counts"]) == run.SHOTS
        assert (
            r["simultaneous_interval"][0]
            <= v["exact_pole_weight"]
            <= r["simultaneous_interval"][1]
        )
        np.testing.assert_allclose(
            r["observed_quadrature"], v["observed_quadrature"], atol=5e-5
        )
    for weight, span in zip(
        np.linalg.eigvalsh(c.conj().T @ p @ c),
        fresh["interval_certificate"]["exact_weight_intervals"],
        strict=True,
    ):
        assert 0 < span[0] < weight < span[1] < 1
    assert (
        np.max(
            abs(
                np.array(fresh["measured_effect"]["real"])
                - np.array(fresh["exact_effect"]["real"])
            )
        )
        < 0.001
    )
    links = run.mag.links(1)
    b = np.linalg.eigh(p)[1][:, -2:]
    for i, g in enumerate(run.rotations()):
        material = run.mag.turn_word(g, run.mag.compensated_turn(g, links), owner=owner)
        remain = np.linalg.eigvalsh(b.conj().T @ material.conj().T @ p @ material @ b)
        np.testing.assert_allclose(
            remain, record["all_rotation_conditional_retention"][i], atol=3e-12
        )
    for row in record["turn_banks"].values():
        for bank in row["banks"]:
            assert (
                len(bank["complete_counts"]) == 21
                and sum(bank["complete_counts"]) == run.SHOTS
            )
            assert (
                abs(bank["complete_counts"][-1] / run.SHOTS - bank["both_pass"])
                < 0.0004
            )
    board = run.pole.setup()[0]
    checks = fresh["compact_complement_checks"]
    for port, route in enumerate(checks["checked_routes"], 2):
        assert route[0] == 12 * owner + port
        assert checks["double_reference_port"] in board.neighbours[route[-1]]
        assert all(b in board.neighbours[a] for a, b in zip(route, route[1:]))
    equal = record["unbiased_transfer"]
    assert equal["conditional_reference_fidelity"] > 1 - 1e-12
    for bank in equal["complete_banks"]:
        assert len(bank["outcomes"]) == 12 and sum(bank["outcomes"]) == run.SHOTS
        assert abs(bank["outcomes"][-1] / run.SHOTS - equal["common_success"]) < 0.00006
    assert equal["common_donor_ports"] == [72, 73]
    assert equal["common_donor_residue"] < 2e-12
    for bank in equal["common_donor_banks"]:
        assert len(bank["outcomes"]) == 3 and sum(bank["outcomes"]) == run.SHOTS
        assert abs(bank["outcomes"][-1] / run.SHOTS - equal["common_success"]) < 0.00006
