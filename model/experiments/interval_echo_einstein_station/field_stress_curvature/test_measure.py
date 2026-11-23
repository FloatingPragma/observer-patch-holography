import hashlib
import json

import numpy as np
from scipy.optimize import linprog
from scipy.sparse.linalg import expm_multiply

from .measure import HERE, SHOTS, THETA, components, metric, ring_laplacian, search


def test_recorded_curvature_contrast_with_independent_primal_solver():
    data = json.loads((HERE / "results.json").read_text())
    held = data["selected"]
    distance = metric(np.array(held["held_vacuum_calibration_counts"]) / SHOTS)
    i, j = held["ray"]
    kernels = np.array(held["held_source_counts"]).mean(axis=0) / SHOTS
    n = len(distance)
    origins = np.repeat(np.eye(n), n, axis=1)
    destinations = np.tile(np.eye(n), n)
    curvatures = []
    for branch in kernels[:2]:
        solution = linprog(
            distance.ravel(),
            A_eq=np.vstack((origins, destinations)),
            b_eq=np.r_[branch[i], branch[j]],
            bounds=(0, None),
            method="highs",
        )
        assert solution.success
        curvatures.append(1 - solution.fun / distance[i, j])
    np.testing.assert_allclose(
        curvatures, held["curvature_estimates"][2][:2], atol=1e-11, rtol=0
    )
    difference = curvatures[1] - curvatures[0]
    assert difference > held["conservative_simultaneous_99_percent_radius"]


def test_ring_history_and_phase_sensitive_control():
    patch = search.Patch()
    words = np.array(patch.field_component())
    field_h = ring_laplacian(patch, words)
    rings = components(patch, words)
    ground = np.ones(len(words)) / np.sqrt(len(words))
    z = np.where(words & patch.bit[4, 11], -1, 1)
    excited = expm_multiply(-0.5j * field_h, z * ground)
    waves = np.array(
        [np.cos(THETA) * ground + s * 1j * np.sin(THETA) * excited for s in (-1, 1)]
    ).T
    np.testing.assert_allclose(
        waves.conj().T @ ground, (waves.conj().T @ ground).conj()[::-1], atol=1e-13
    )
    for click1, click2 in ((False, False), (False, True), (True, False), (True, True)):
        projected = rings[10] @ waves
        branch = projected if click1 else waves - projected
        later = expm_multiply(-0.31j * field_h, branch)
        second = rings[12] @ later
        endpoints = second if click2 else later - second
        law = np.sum(abs(endpoints) ** 2, axis=0)
        np.testing.assert_allclose(law[0], law[1], atol=2e-14, rtol=0)
    assert abs(np.sum(z * (abs(waves[:, 1]) ** 2 - abs(waves[:, 0]) ** 2))) > 0.01


def test_replay_counts_and_registered_criteria():
    primary = json.loads((HERE / "results.json").read_text())
    replay = json.loads((HERE / "replay.json").read_text())
    for data, name in ((primary, "PROTOCOL.md"), (replay, "REPLAY_PROTOCOL.md")):
        assert (
            data["protocol_sha256"]
            == hashlib.sha256((HERE / name).read_bytes()).hexdigest()
        )
    trials = replay["cohort"]
    for name in (
        "calibration_counts",
        "vacuum_deadline_counts",
        "incoherent_flight_counts",
    ):
        assert np.all(np.array(trials[name]).sum(axis=-1) == SHOTS)
    assert not any(trials["each_probe_passes"])
    assert replay["continuity"]["resolved"] and replay["continuity"]["agrees"]
    assert sum(replay["continuity"]["current_counts"]) == SHOTS


def test_addressed_ring_interferometry_against_energy_change():
    result = json.loads((HERE / "interference.json").read_text())
    assert (
        result["protocol_sha256"]
        == hashlib.sha256((HERE / "INTERFERENCE_PROTOCOL.md").read_bytes()).hexdigest()
    )
    probabilities = np.array(result["pulse_and_detector_probabilities"])
    endpoint = result["slope_detector_probabilities"]
    net = ((probabilities[:-1, 0] - probabilities[:-1, 1]) / 2).sum()
    np.testing.assert_allclose(
        net, (endpoint[1] - endpoint[0]) * 16, atol=2e-14, rtol=0
    )
    assert (
        result["sign_resolved"]
        and result["disjoint_control_passes"]
        and result["continuity_agrees"]
    )
