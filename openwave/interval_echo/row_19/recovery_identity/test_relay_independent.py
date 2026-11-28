"""Reconstruct the receiver operators from the original 99-state flight restriction."""

import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_einstein_station.mobile_collar_source import (
    measure as native,
)

HERE = Path(__file__).parent


def single_written_edge(board, receiver, q):
    # The original two-contact solver includes a spectator. Project it on coin
    # zero, then replace precisely its four spatial blocks by those from a second
    # flight having the spectator elsewhere. Both restrictions start with the
    # raw 99-state operator. Remove the common scalar at an untouched address.
    unavailable = set(receiver)
    spectators = []
    for a, b in board.seams:
        if a // 12 == b // 12 and not unavailable.intersection((a, b)):
            spectators.append((a, b))
            unavailable.update((a, b))
            if len(spectators) == 2:
                break
    anchor = 2 * next(p for p in range(board.n) if p not in unavailable)
    candidates = []
    for edge in spectators:
        full = native.hamiltonian(board, (*receiver, *edge), q, 8)
        one = full[::2, ::2].tocsr()
        one -= one[anchor, anchor] * sparse.eye(one.shape[0], format="csr")
        candidates.append(one)
    a, b = candidates
    sites = [2 * p + coin for p in spectators[0] for coin in (0, 1)]
    mask = np.zeros(a.shape[0])
    mask[sites] = 1
    selector = sparse.diags(mask)
    doubled = 2 * (a + selector @ (b - a) @ selector)
    np.testing.assert_allclose(doubled.data, np.rint(doubled.data), atol=2e-12, rtol=0)
    doubled.data = np.rint(doubled.data).astype(np.int64)
    return doubled


def launch(board, path):
    rails = np.zeros((board.n, 2), dtype=int)
    rails[path[:2], 0] = 1
    rails[path[1:], 1] = 1
    return np.kron(rails, np.eye(2, dtype=int))


def test_receipt_with_original_contact_restriction():
    receipt = json.loads((HERE / "relay.json").read_text())
    precision = json.loads((HERE / "relay_precision.json").read_text())
    board = native.federation(1)
    route = native.region(board)
    radii = (
        shortest_path(
            board.adjacency(), directed=False, unweighted=True, indices=route[1:3]
        )
        .min(axis=0)
        .astype(int)
    )
    x = sparse.diags(np.repeat(radii, 2), dtype=np.int64)
    operators = []
    for station in receipt["stations"]:
        path = station["path"]
        np.testing.assert_array_equal(radii[path], station["radii"])
        launch_columns = launch(board, path)
        reports = station["appointments"]["8"]
        by_charge = []
        for q in (-1, 1):
            k = single_written_edge(board, path[:2], q)
            commit = k @ x - x @ k
            second = commit @ k - k @ commit
            projected = launch_columns.T @ (second @ launch_columns)
            # K=2H, and the raw two-site rails have squared norm two.
            eight_a = projected[2:, 2:] - projected[:2, :2]
            np.testing.assert_allclose(
                eight_a / 8,
                reports["probes"][str(q)]["acceleration_operator"],
                atol=3e-12,
                rtol=0,
            )
            by_charge.append(eight_a)
        operators.append(by_charge)
    np.testing.assert_array_equal(operators, precision["operators_times_8"])
    # Exact matrix entries for (A_- + A_+)/2 and for its traceless part.
    near, far = (sum(pair) for pair in operators)
    np.testing.assert_array_equal(near, [[-236, -16], [-16, -236]])
    np.testing.assert_array_equal(far, -232 * np.eye(2, dtype=int))


def test_near_eight_clock_contrast_from_receiver_probabilities():
    receipt = json.loads((HERE / "relay.json").read_text())
    station = receipt["stations"][0]
    path = station["path"]
    report = station["appointments"]["8"]
    bloch = np.asarray(report["bloch"])
    differences = bloch[1::2] - bloch[::2]
    steps = receipt["half_steps"]
    divisors = [2 * steps["ordinary"], 2 * steps["held"], 1, 1]
    release = -differences[:, 0] / divisors
    np.testing.assert_allclose(release[1] / release[0], -0.01766866947, rtol=1e-8)

    board = native.federation(1)
    route = native.region(board)
    radii = shortest_path(
        board.adjacency(), directed=False, unweighted=True, indices=route[1:3]
    ).min(axis=0)
    x = np.repeat(radii - radii[path[0]], 2)[:, None]
    columns = launch(board, path) / np.sqrt(2)
    per_charge = []
    for q in (-1, 1):
        k = single_written_edge(board, path[:2], q)
        relative = []
        for t in (0, 0.25, 0.5):
            wave = expm_multiply(
                (-0.5j * t) * k, columns, traceA=(-0.5j * t) * k.diagonal().sum()
            )
            np.testing.assert_allclose(
                np.diag(wave.conj().T @ wave), 1, atol=1e-12, rtol=0
            )
            read = wave.conj().T @ (x * wave)
            relative.append(read[2:, 2:] - read[:2, :2])
        stencil = 16 * (relative[0] - 2 * relative[1] + relative[2])
        pauli = np.array([np.trace(stencil @ sigma).real / 2 for sigma in native.PAULI])
        slopes = (differences @ pauli) / divisors
        old = report["probes"][str(q)]["finite_prediction"]
        np.testing.assert_allclose(
            slopes[:2], [old["ordinary"], old["held"]], atol=7e-10, rtol=0
        )
        assert slopes[1] / slopes[0] < 0
        per_charge.append(slopes)
    mean = np.mean(per_charge, axis=0)
    np.testing.assert_allclose(mean[1] / mean[0], -0.022052, atol=1e-6, rtol=0)
    assert mean[0] > 0 and np.all(mean[1:] < 0)
    source = receipt["source"]
    error = source["held_contrast_nats"]["halfwidth"]
    estimate = source["held_contrast_nats"]["estimate"]
    source_ratio = (estimate + np.array([-error, error])) / (
        2 * steps["held"] * source["ordinary_modular_slope"]
    )
    np.testing.assert_allclose(source_ratio, [0.03502848867, 0.03650307523], rtol=1e-9)
