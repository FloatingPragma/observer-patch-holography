"""Exterior insertion, its two-station interference read, and two-mode rotations."""

from itertools import combinations

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.face_frame import (
    shared_report,
)
from model.executable_model import federation
from model.experiments.interval_echo_consensus.signed_effects import record_transfer

from . import material


def at(operator, j):
    return np.kron(np.eye(2 ** (4 - j)), np.kron(operator, np.eye(2**j)))


def endpoint_for_unit_branch(state):
    weights = abs(state) ** 2
    which = int(np.argmax(weights))
    np.testing.assert_allclose(weights, np.eye(len(state))[which], rtol=0, atol=2e-12)
    return which


def insertion_trials(sharp_reads):
    flip, phase, hadamard = material.dimer_gates()
    # The vacuum response fixes the phase convention of each mode.
    flip /= flip[0, 1]
    phase /= phase[0, 0]
    eye = np.eye(32)
    occupation = np.diag([0, 1])
    add, ordinary = [], []
    prefix = eye.copy()
    for j in range(5):
        vacancy = eye - at(occupation, j)
        ordinary.append(at(flip, j) @ vacancy)
        add.append(ordinary[j] @ prefix)
        prefix = at(phase, j) @ prefix
    for i in range(5):
        for j in range(5):
            np.testing.assert_allclose(add[i] @ add[j] + add[j] @ add[i], 0, atol=2e-13)
            remove = add[j].conj().T
            np.testing.assert_allclose(
                add[i] @ remove + remove @ add[i], int(i == j) * eye, atol=2e-13
            )

    # All ten physical rails are permuted with no phases on the full load sector.
    # Read the code permutation from the sharp-flight calibration.
    delivery = eye[:, [material.decode(r["receiver"]) for r in sharp_reads]]
    records = []
    for i, j in combinations(range(5), 2):
        rotate = at(hadamard, i)
        for base in range(32):
            if base & ((1 << i) | (1 << j)):
                continue
            launch = rotate @ eye[:, base]
            clear = sharp_reads[base]["sender"][4 * j : 4 * j + 4]
            row = [i, j, base, clear]
            for pulse in (add[j], ordinary[j]):
                received = rotate.conj().T @ delivery @ pulse @ launch
                out = endpoint_for_unit_branch(received)
                row.append(sharp_reads[out]["receiver"])
            records.append(row)
    refused = []
    for j in range(5):
        for base in range(32):
            if base & (1 << j):
                np.testing.assert_allclose(add[j] @ eye[:, base], 0, atol=2e-13)
                code = sharp_reads[base]["sender"]
                refused.append([j, base, code[4 * j : 4 * j + 4]])
    return records, refused


def connector_routes(geometry):
    board = federation(0)
    sender = geometry["owners"]["bob"]
    graph = nx.Graph(board.seams).subgraph(range(12 * sender, 12 * (sender + 1)))
    gate = min(
        (min(a, b), max(a, b))
        for a, b in graph.edges
        if board.degree[a] == board.degree[b] == 5
    )
    routes = []
    for i, j in ((0, 1), (1, 2), (3, 4)):
        labels = list(range(board.n))
        trails = []
        for mode, target in zip((i, j), gate, strict=True):
            label = geometry["dimers"]["bob"][mode][1]
            trail = shared_report.path(graph, labels.index(label), target)
            shared_report.transpose(labels, trail)
            trails.append(trail)
        assert [labels[p] for p in gate] == [
            geometry["dimers"]["bob"][m][1] for m in (i, j)
        ]
        for trail in reversed(trails):
            shared_report.transpose(labels, trail)
        assert labels == list(range(board.n))
        routes.append(dict(modes=[i, j], connector=gate, transfers=trails))
    return routes


def rotations():
    _, _, h = material.dimer_gates()
    law = candidate.native_law
    ex, con = law.local_term_parts()
    full = ex + candidate.KAPPA * con
    # Two marked ports, logical 0 at load one and logical 1 at a vacancy.
    local = np.eye(law.DIM)[:, [law.SRC.index((a, b)) for a in (1, 0) for b in (1, 0)]]
    contact = local.T @ full @ local
    np.testing.assert_allclose(full @ local, local @ contact, atol=2e-13)
    ref = candidate.KAPPA * np.diag([1, 0])
    refs = [np.kron(ref, np.eye(2)), np.kron(np.eye(2), ref)]
    cz = np.eye(4)
    for kind, a, _, dwell in record_transfer.controlled_phase(0, 1):
        generator = contact if kind == "contact" else refs[a]
        cz = expm(-1j * np.pi * float(dwell) * generator) @ cz
    np.testing.assert_allclose(cz, -np.diag([1, 1, 1, -1]), atol=2e-13)
    h_a, h_b = np.kron(h, np.eye(2)), np.kron(np.eye(2), h)
    cx_a = h_b.conj().T @ cz @ h_b
    cx_b = h_a.conj().T @ cz @ h_a
    dimer = contact[np.ix_((1, 2), (1, 2))]
    quarter = expm(-1j * np.pi / 16 * ref)
    examples = []
    for numerator, denominator in ((1, 5), (-1, 7)):
        angle = np.pi * numerator / denominator
        # The cross-parity rotation uses a half-angle on one dimer.
        dwell = (-angle) % (2 * np.pi)
        turn = quarter @ expm(-1j * dwell * dimer) @ quarter.conj().T
        turn_a = np.kron(turn, np.eye(2))
        conditioned = cx_b @ turn_a.conj().T @ cx_b @ turn_a
        rotate = cx_a @ conditioned @ cx_a
        expected = np.eye(4)
        c, s = np.cos(angle), np.sin(angle)
        expected[1:3, 1:3] = [[c, -s], [s, c]]
        np.testing.assert_allclose(rotate, expected, atol=2e-13)
        examples.append(
            dict(
                angle_over_pi=[numerator, denominator],
                probabilities_by_input=[
                    [round(float(p), 12) for p in col] for col in (abs(rotate) ** 2).T
                ],
            )
        )
    return examples
