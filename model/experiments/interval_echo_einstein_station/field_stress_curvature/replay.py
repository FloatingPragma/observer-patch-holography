"""Independently count the full energy–curvature cohort and a local energy continuity witness."""

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from .measure import (
    HERE,
    SHOTS,
    THETA,
    Transport,
    components,
    metric,
    multinomial_rows,
    ring_laplacian,
    search,
    signs_and_pooled,
)


def cohort(primary, rng):
    energy = np.array(primary["energy_profiles"])
    clicks = rng.binomial(SHOTS, np.clip(energy, 0, 1))
    calibration = multinomial_rows(np.array(primary["vacuum_calibration_kernels"]), rng)
    quiet = multinomial_rows(np.array(primary["vacuum_deadline_kernels"]), rng)
    pulses = multinomial_rows(np.array(primary["phase_flights"])[:, :, 2], rng)
    reader = Transport(metric(calibration / SHOTS))
    edges = [tuple(e) for e in primary["head_order"]]
    face_edges = [
        [edges.index(tuple(sorted(edge))) for edge in search.around(face)]
        for face in primary["face_order"]
    ]
    baseline = np.array(
        [reader.curvature(k, edges) for k in signs_and_pooled(quiet) / SHOTS]
    )
    readings = np.array(
        [
            [reader.curvature(k, edges) for k in heads]
            for heads in signs_and_pooled(pulses) / SHOTS
        ]
    )
    face_response = np.take(readings - baseline[:, None], face_edges, axis=-1).mean(
        axis=-1
    )
    train = primary["local_fit"]["trained_head_indices"]
    held = primary["local_fit"]["held_head_indices"]
    observed_energy = clicks / SHOTS
    coefficient = (observed_energy[train] * face_response[2, train]).sum() / (
        observed_energy[train] ** 2
    ).sum()
    forecast = coefficient * observed_energy[held]
    errors = [
        np.linalg.norm(values[held] - forecast) / np.linalg.norm(values[held])
        for values in face_response
    ]
    return dict(
        calibration_counts=calibration,
        vacuum_deadline_counts=quiet,
        incoherent_flight_counts=pulses,
        ring_energy_clicks=clicks,
        training_coefficient=coefficient,
        face_response=face_response,
        relative_held_errors=errors,
        each_probe_passes=[bool(v <= 0.05) for v in errors],
    )


def local_continuity(primary, rng):
    patch = search.Patch()
    words = np.array(patch.field_component())
    ground = np.ones(len(words)) / np.sqrt(len(words))
    field_h = ring_laplacian(patch, words)
    operators = components(patch, words)
    head = primary["selected"]["head_index"]
    face = int(np.abs(np.array(primary["energy_currents"])[head].sum(axis=0)).argmax())
    ring = operators[face]
    rotated = np.where(words & patch.bit[patch.edges[head]], -1, 1) * ground
    theta_state = np.cos(THETA) * ground + 1j * np.sin(THETA) * rotated
    center_state = expm_multiply(-0.5j * field_h, theta_state)
    transfer = 1j * (field_h @ ring - ring @ field_h)
    eigenvalues, eigenvectors = eigh(transfer.toarray())
    levels, groups = np.unique(np.round(eigenvalues, 12), return_inverse=True)
    half_window = 1 / 32
    ends = [
        expm_multiply(-1j * offset * field_h, center_state)
        for offset in (-half_window, half_window)
    ]
    ring_prob = np.array([np.linalg.norm(ring @ state) ** 2 for state in ends])
    ring_counts = rng.binomial(SHOTS, ring_prob)
    nodes, weights = leggauss(16)
    mean_prob = np.zeros(len(levels))
    for time, weight in zip(half_window * nodes, weights / 2):
        state = expm_multiply(-1j * time * field_h, center_state)
        mean_prob += weight * np.bincount(
            groups, abs(eigenvectors.conj().T @ state) ** 2
        )
    current_counts = rng.multinomial(SHOTS, mean_prob / mean_prob.sum())
    energy_slope = float(np.diff(ring_counts)[0]) / (2 * half_window * SHOTS)
    mean_current = current_counts @ levels / SHOTS
    finite_slope = float(np.diff(ring_prob)[0]) / (2 * half_window)
    np.testing.assert_allclose(finite_slope, mean_prob @ levels, atol=1e-12, rtol=0)
    eps = np.sqrt(np.log(6 / 0.01) / (2 * SHOTS))
    radius_slope = eps / half_window
    quad_error = 2 * (40 * half_window) ** 32 / factorial(32)
    radius_current = np.ptp(levels) * eps + quad_error * abs(levels).sum() + 5.1e-13
    return dict(
        head_index=head,
        face_index=face,
        half_appointment_width=half_window,
        endpoint_ring_probabilities=ring_prob,
        click_counts=ring_counts,
        current_levels=levels,
        current_probabilities=mean_prob,
        current_counts=current_counts,
        expected_windowed_current=finite_slope,
        measured_energy_slope=energy_slope,
        measured_windowed_current=mean_current,
        simultaneous_99_percent_slope_radius=radius_slope,
        simultaneous_99_percent_current_radius=radius_current,
        uniform_appointment_probability_quadrature_bound=quad_error,
        resolved=bool(
            abs(energy_slope) > radius_slope and abs(mean_current) > radius_current
        ),
        agrees=bool(abs(energy_slope - mean_current) <= radius_slope + radius_current),
    )


def run():
    primary_path = HERE / "results.json"
    primary = json.loads(primary_path.read_text())
    seed = 517813
    rng = np.random.default_rng(seed)
    return dict(
        seed=seed,
        independent_simulated_shots_per_setting=SHOTS,
        protocol_sha256=hashlib.sha256(
            (HERE / "REPLAY_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        primary_protocol_sha256=primary["protocol_sha256"],
        cohort=cohort(primary, rng),
        continuity=local_continuity(primary, rng),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = run()
    with args.output.open("x") as output:
        output.write(
            json.dumps(
                data,
                default=lambda x: x.tolist(),
                allow_nan=False,
                separators=(",", ":"),
            )
            + "\n"
        )
