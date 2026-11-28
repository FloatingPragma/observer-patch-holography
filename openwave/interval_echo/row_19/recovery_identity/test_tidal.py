"""Independent checks of release blindness and a counted finite-time comparison."""

import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from openwave.interval_echo.row_19.recovery_identity import tidal


def test_disjoint_contacts_and_a_joint_potential():
    # The two mixtures have equal one-head marginals and distinct ZZ.
    delta = np.diag([1, -1, -1, 1]) / 4
    z = np.diag([1, -1])
    identity = np.eye(2)
    edge = np.diag(np.ones(3), 1)
    hopping = edge + edge.T
    left = np.zeros((4, 4))
    left[:2, :2] = [[1, 2], [2, -1]]
    right = left[::-1, ::-1].copy()
    h = (
        np.kron(hopping, np.eye(4))
        + np.kron(left, np.kron(z, identity))
        + np.kron(right, np.kron(identity, z))
    )
    joint_contact = np.kron(np.diag([1, 0, 0, 0]), np.kron(z, z))
    psi = np.array([1, 1, 0, 0]) / np.sqrt(2)
    source_difference = np.kron(np.outer(psi, psi), delta)
    x = np.kron(np.diag([0, 1, 2, 3]), np.eye(4))

    def read(generator):
        v = generator @ x - x @ generator
        return np.trace(source_difference @ (v @ generator - generator @ v))

    np.testing.assert_allclose(read(h), 0, atol=1e-14)
    assert abs(read(h + joint_contact)) > 0.9


def test_forward_pure_flights_recover_held_outer_minus_inner():
    # Diagonalize each actual stock, propagate the pure components separately
    # under a second assembly of the floating native seam Hamiltonian.
    stock, steps, _, _ = tidal.transfer.inputs()
    ranging = tidal.front.optics.ranging
    board = ranging.formation.federation(1)
    route = ranging.formation.region(board)
    radii, paths, _ = tidal.chart(board, route)
    observed = np.zeros((2, 2, 3))
    for parity in (0, 1):
        h = ranging.generator(board, route, -1, 8, bool(parity))
        for s, key in enumerate(("held_minus", "held_plus")):
            weights, vectors = np.linalg.eigh(stock[key][parity])
            for ray in (0, 1):
                wave = np.zeros((board.n, 4, 4), complex)
                wave[paths[0][ray : ray + 2]] = vectors / np.sqrt(2)
                evolved = np.stack(
                    [
                        expm_multiply(-1j * time * h, wave.reshape(-1, 4), traceA=0)
                        for time in (1.75, 2.0, 2.25)
                    ]
                )
                np.testing.assert_allclose(
                    np.sum(abs(evolved) ** 2, axis=1), 1, atol=2e-12, rtol=0
                )
                probs = (abs(evolved.reshape(3, board.n, 4, 4)) ** 2).sum(
                    axis=2
                ) @ weights
                observed[s, ray] += probs @ radii
    slope = np.diff(observed, axis=0).squeeze()
    relative = slope[1] - slope[0]
    accel = (relative[0] - 2 * relative[1] + relative[2]) / (
        0.25**2 * 2 * steps["held"]
    )
    saved = json.loads(tidal.RECEIPT.read_text())["banks"]["-1"]
    np.testing.assert_allclose(accel, saved["population"]["2.0"][0]["held"], atol=2e-9)
    reading = saved["sampled"]["2.0"][0]
    assert abs(accel - reading["held"]) <= reading["held_radius"]
