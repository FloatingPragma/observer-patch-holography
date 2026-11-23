"""Measure local energy currents with addressed ring pulses and binary dark-port clicks."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import expm_multiply

from .measure import HERE, SHOTS, THETA, ring_laplacian, search


def flips(patch, words):
    lookup = {int(w): i for i, w in enumerate(words)}
    perms = np.tile(np.arange(len(words)), (len(patch.faces), 1))
    for i, word in enumerate(words):
        for face, other in patch.ring_moves(int(word)):
            perms[face, i] = lookup[other]
    np.testing.assert_array_equal(
        np.take_along_axis(perms, perms, axis=1),
        np.broadcast_to(np.arange(len(words)), perms.shape),
    )
    return perms


def pulse(state, perm, angle):
    """Symmetric amplitude stays fixed, antisymmetric amplitude acquires a phase."""
    return (state + state[perm]) / 2 + np.exp(-1j * angle) * (state - state[perm]) / 2


def dark_port(state, perm):
    return np.linalg.norm((state - state[perm]) / 2) ** 2


def run():
    seed = 517814
    rng = np.random.default_rng(seed)
    patch = search.Patch()
    words = np.array(patch.field_component())
    h = ring_laplacian(patch, words)
    permutations = flips(patch, words)
    head = (4, 11)
    b = patch.faces.index((2, 4, 11))
    neighbors = [
        i
        for i, face in enumerate(patch.faces)
        if len(set(face) & set(patch.faces[b])) == 2
    ]
    control = next(
        i for i, face in enumerate(patch.faces) if not set(face) & set(patch.faces[b])
    )
    assert len(neighbors) == 3
    ground = np.ones(len(words)) / np.sqrt(len(words))
    source = np.exp(1j * THETA * np.where(words & patch.bit[head], -1, 1)) * ground
    half_window = 1 / 32
    center = expm_multiply(-0.5j * h, source)
    ends = [
        expm_multiply(-1j * offset * h, center)
        for offset in (-half_window, half_window)
    ]
    before_after = [dark_port(state, permutations[b]) for state in ends]
    endpoint_clicks = rng.binomial(SHOTS, before_after)
    nodes, weights = leggauss(16)
    appointment_states = [
        expm_multiply(-1j * t * h, center) for t in half_window * nodes
    ]
    probabilities = []
    for a in neighbors + [control]:
        probabilities.append(
            [
                sum(
                    weight
                    * dark_port(pulse(state, permutations[a], angle), permutations[b])
                    for weight, state in zip(weights / 2, appointment_states)
                )
                for angle in (np.pi / 2, 3 * np.pi / 2)
            ]
        )
    clicks = rng.binomial(SHOTS, probabilities)
    currents = (clicks[:, 0] - clicks[:, 1]) / (2 * SHOTS)
    net = sum(currents[:-1])
    slope = (endpoint_clicks[1] - endpoint_clicks[0]) / (2 * SHOTS * half_window)
    n_rates = 2 + len(probabilities) * 2
    eps = np.sqrt(np.log(2 * n_rates / 0.01) / (2 * SHOTS)) + 1e-31
    net_radius, slope_radius = len(neighbors) * eps, eps / half_window
    return dict(
        seed=seed,
        shots_per_setting=SHOTS,
        head=head,
        face=b,
        neighbors=neighbors,
        disjoint_control=control,
        protocol_sha256=hashlib.sha256(
            (HERE / "INTERFERENCE_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        slope_detector_probabilities=before_after,
        slope_detector_clicks=endpoint_clicks,
        pulse_and_detector_probabilities=probabilities,
        pulse_and_detector_clicks=clicks,
        measured_face_currents=currents,
        measured_net_influx=net,
        measured_energy_slope=slope,
        simultaneous_99_percent_net_radius=net_radius,
        simultaneous_99_percent_slope_radius=slope_radius,
        simultaneous_99_percent_single_current_radius=eps,
        sign_resolved=bool(abs(net) > net_radius and abs(slope) > slope_radius),
        continuity_agrees=bool(abs(net - slope) <= net_radius + slope_radius),
        disjoint_control_passes=bool(abs(currents[-1]) <= eps),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as output:
        output.write(
            json.dumps(
                run(),
                default=lambda x: x.tolist(),
                allow_nan=False,
                separators=(",", ":"),
            )
            + "\n"
        )
