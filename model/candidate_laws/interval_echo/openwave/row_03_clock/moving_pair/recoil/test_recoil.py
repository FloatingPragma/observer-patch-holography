"""Microscopic apparatus checks, an independent finite-board flow identity, and raw record tests."""

import json

import numpy as np
from scipy import sparse
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_03_clock.moving_pair.recoil import (
    measure,
)
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import federation
from model.experiments.material_inertia.neutral_defect import generators, ordered_index


def test_local_apparatus_and_owner_phase():
    h = law.seam_term(candidate.KAPPA)
    permutation = np.array([law.SRC.index((b, a)) for a, b in law.SRC])
    np.testing.assert_allclose(
        expm(-1j * np.pi * h), np.eye(9)[permutation], atol=8e-14
    )
    unitary = expm(-1j * np.pi * h / 2)
    for n in (1, 2):
        positions = [law.SRC.index(p) for p in [(n, 0), (0, n)]]
        block = unitary[np.ix_(positions, positions)]
        ideal = (np.eye(2) + 1j * np.fliplr(np.eye(2))) / np.sqrt(2)
        np.testing.assert_allclose(block, ideal * block[0, 0] / ideal[0, 0], atol=8e-14)
    board = federation(0)
    for start, _ in measure.fixtures(board):
        owner = start[0] // 12
        lookup = measure.deposits(board, owner)
        for a, ref in lookup:
            assert board.degree[ref] == 5 and ref // 12 == owner
            assert a == ref or ref in board.neighbours[a]
        for theta in (measure.ANGLE, 2 * measure.ANGLE):
            np.testing.assert_allclose(
                np.exp(-1j * theta * np.array([0, 1, 4])),
                np.exp(-1j * theta * np.array([0, 1, 1])),
                atol=2e-14,
            )


def test_native_generator_independently_predicts_meter_and_energy():
    board = federation(0)
    x, y, adjacent, kinetic, _, _ = generators(board)
    h = kinetic - sparse.diags(4 * candidate.KAPPA * adjacent)
    outside = (x // 12 != 0).astype(int) + (y // 12 != 0)
    cuts = measure.cut_indices(board, 0)
    rng = np.random.default_rng(7421)
    # Include all two-load configurations; deliberate gateway coherence makes
    # the signal well above roundoff and unrelated to the archived source states.
    w = rng.normal(size=len(x)) + 1j * rng.normal(size=len(x))
    for left, right in cuts:
        w[right] += (0.5 + 0.7j) * w[left]
    w /= np.linalg.norm(w)
    energies, currents = [], []
    for phi in (-measure.ANGLE, 0, measure.ANGLE):
        u = np.exp(1j * phi * outside) * w
        du = -1j * (h @ u)
        j = 2 * np.vdot(outside * u, du).real
        p = measure.boundary_probabilities(u, cuts, "current")
        np.testing.assert_allclose(1.5 * (p[0] - p[1]), j, atol=3e-13)
        energies.append(float(np.vdot(u, h @ u).real))
        currents.append(j)
    q = measure.boundary_probabilities(w, cuts, "parity")
    kinetic_coherence = 1.5 * (q[0] - q[1])
    np.testing.assert_allclose(
        (currents[2] - currents[0]) / 2,
        np.sin(measure.ANGLE) * kinetic_coherence,
        atol=3e-13,
    )
    np.testing.assert_allclose(
        (energies[2] + energies[0]) / 2 - energies[1],
        (1 - np.cos(measure.ANGLE)) * kinetic_coherence,
        atol=3e-13,
    )
    for start, _ in measure.fixtures(board):
        col = ordered_index(*start, board.n)
        assert kinetic[col, col] + 208 == 240 + h[col, col] == 213


def test_archived_observer_counts_and_same_energy_comparison():
    archive = json.loads((measure.HERE / "records.json").read_text())
    for row in archive["settings"].values():
        found = measure.evaluate(row)
        pulse = measure.calibrated_sine(archive["phase_dimer_transfer_stay"])
        prediction = measure.interval_product(pulse, found["kinetic_coherence"])
        residual = measure.minus(np.array(found["opposite_kick_current"]), prediction)
        assert residual[0] <= 0 <= residual[1]
        for label in ("forecast_residual", "cosine_residual"):
            lo, hi = found[label]
            assert lo <= 0 <= hi
        for sign in ("plus", "minus", "sham"):
            assert all(sum(b) == measure.SHOTS for b in row["flight"][sign])
    native_lower = min(
        row["adjacent_lower_after_either_kick"]
        for key, row in archive["settings"].items()
        if key.startswith("native/")
    )
    assert 0.84319 < native_lower < 0.84320
    for i in range(3):
        for tag in ("sharp", "quarter", "reverse_quarter"):

            def read(law, age):
                return measure.evaluate(archive["settings"][f"{law}/{i}/{age}/{tag}"])

            unexpanded = read("native", 0)
            assert (
                unexpanded["opposite_kick_current"][0]
                < 0
                < unexpanded["opposite_kick_current"][1]
            )
            for age in (1, 4):
                active, control = read("native", age), read("exchange_only", age)
                assert 0 < active["opposite_kick_current"][0]
                assert (
                    active["opposite_kick_current"][1]
                    < control["opposite_kick_current"][0]
                )
            assert read("native", 1)["departure"][1]["centre_difference"][0] > 0
            assert read("native", 4)["departure"][1]["centre_difference"][1] < 0
            assert read("exchange_only", 4)["departure"][1]["centre_difference"][0] > 0
        sensitive = [
            measure.interval(
                archive["settings"][f"native/{i}/4/{tag}"]["species_selective_current"],
                [1.5, -1.5, 0],
            )
            for tag in ("quarter", "reverse_quarter")
        ]
        assert sensitive[0][0] > sensitive[1][1]


def test_held_intact_arrivals_by_independent_exponential():
    from scipy.sparse.linalg import expm_multiply

    board = federation(0)
    x, y, adjacent, kinetic, _, _ = generators(board)
    h = kinetic - sparse.diags(4 * candidate.KAPPA * adjacent)
    source, receiver = measure.fixtures(board)[0]
    a, b = source
    state = np.zeros(len(x), dtype=complex)
    state[ordered_index(a, b, board.n)] = (1 - 1j) / 2
    state[ordered_index(b, a, board.n)] = (1 + 1j) / 2
    expanded = expm_multiply(-1j * h, state)
    # Apply the individual material-square phase at the twelve emitter sites.
    inside_square = (x // 12 == a // 12).astype(int) + 4 * (y // 12 == a // 12)
    nearby = (x // 12 == receiver[0] // 12) & (y // 12 == receiver[0] // 12) & adjacent
    archive = json.loads((measure.HERE / "records.json").read_text())
    record = archive["settings"]["native/0/1/quarter"]
    for sign, theta in [("plus", measure.ANGLE), ("minus", -measure.ANGLE)]:
        final = expm_multiply(-0.5j * h, expanded * np.exp(-1j * theta * inside_square))
        expected = np.sum(abs(final[nearby]) ** 2)
        band = measure.interval(record["flight"][sign][1], [0, 0, 0, 0, 1])
        assert band[0] <= expected <= band[1]
