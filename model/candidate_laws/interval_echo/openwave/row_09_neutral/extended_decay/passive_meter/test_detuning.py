"""Replay the stored laws and check a direct native flight for both detunings."""

from hashlib import sha256
import json

import numpy as np
from scipy.linalg import expm

from . import calorimeter, detuning
from ..material import route
from ..measure import rebuild
from ...joined_decay.measure import core_pulse
from ...passive_decay.measure import hamiltonians, physical


def test_complete_replay_and_controls():
    archived = json.loads((detuning.HERE / "detuning.json").read_text())
    assert (
        archived["protocol_sha256"]
        == sha256((detuning.HERE / "DETUNING.md").read_bytes()).hexdigest()
    )
    rebuilt = detuning.acquire()
    assert archived["spectator_edges"] == [[72, 77], [73, 81]]
    assert archived["waits"] == list(calorimeter.WAITS)
    for saved, actual in zip(archived["mounts"], rebuilt["mounts"], strict=True):
        assert saved["preparation"] == actual["preparation"]
        assert saved["costs"] == actual["costs"]
        assert actual["preparation"]["total_load"] == 234
        gap = np.hypot(1, saved["costs"][0] - saved["costs"][1])
        np.testing.assert_allclose(saved["clock"]["gap"], gap, atol=1e-10, rtol=0)
        np.testing.assert_allclose(
            saved["clock"]["wait_cosine_sine"],
            actual["clock"]["wait_cosine_sine"],
            atol=1e-10,
            rtol=0,
        )
        for name, arms in saved["settings"].items():
            for arm, joint in arms.items():
                laws = np.asarray(joint)
                assert laws.shape == (6, 4, 2)
                assert laws.min() >= 0
                np.testing.assert_allclose(laws.sum(axis=(1, 2)), 1, atol=1e-10, rtol=0)
                np.testing.assert_allclose(
                    laws, actual["settings"][name][arm], atol=1e-10, rtol=0
                )
            if name != "native":
                pulse, idle = (
                    np.asarray(arms[a])[:, :, 1].sum(axis=1) for a in ("pulse", "idle")
                )
                np.testing.assert_allclose(pulse, idle, atol=1e-10, rtol=0)
            if name == "detector_cut":
                np.testing.assert_allclose(
                    np.asarray(arms["pulse"])[:, :, 1], 0, atol=1e-10
                )


def test_direct_physical_flight_at_primary_appointment():
    record = json.loads((detuning.HERE / "detuning.json").read_text())
    board, dimers, heads, _, _ = calorimeter.apparatus()
    parent = rebuild(board, dimers[:-1], route((3,), crossing=True)[2], False, True)
    _, source = np.linalg.eigh(hamiltonians(4)[0])
    _, guide = np.linalg.eigh(parent)
    stock = np.kron(core_pulse(source), np.eye(8)) @ guide[:, 0]
    for mount in record["mounts"]:
        costs = mount["costs"]
        contacts = dict(zip(map(tuple, record["spectator_edges"]), costs, strict=True))
        diagonal = physical(board, dimers, heads | contacts)
        full = np.diag((diagonal - diagonal[0]).astype(float))
        for bit in range(len(dimers)):
            labels = np.arange(len(full))
            full[labels, labels ^ (1 << bit)] -= 0.5
        meter = np.array([[costs[0], -0.5], [-0.5, costs[1]]])
        _, basis = np.linalg.eigh(meter)
        flight = expm(-32j * full) @ np.kron(stock, basis[:, 0])
        read = np.kron(source.T, np.kron(np.eye(8), basis.T))
        probabilities = abs((read @ flight).reshape(4, 8, 2)) ** 2
        target = mount["settings"]["native"]["pulse"][record["waits"].index(32)]
        np.testing.assert_allclose(
            probabilities.sum(axis=1), target, atol=1e-10, rtol=0
        )


def test_reversing_bias_preserves_isolated_absorption_strength():
    rail = np.diag([0.0, 1.0])
    for bias in (-1, 1):
        h = calorimeter.native_meter() + bias * rail
        levels, v = np.linalg.eigh(h)
        np.testing.assert_allclose(np.ptp(levels), np.sqrt(2), atol=1e-13)
        np.testing.assert_allclose(
            abs(v[:, 1] @ rail @ v[:, 0]) ** 2, 1 / 8, atol=1e-13
        )
