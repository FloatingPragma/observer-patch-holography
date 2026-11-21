"""Calibrated joint detector records and an interacting clock–velocity envelope."""

import math

import numpy as np
from scipy import sparse

from model.executable_model import carrier
from model.experiments.radiative_repair.source import evolve
from ..mobility import PORTS_PER_CELL
from . import clock, geometry, source


G = 1.0
WAIT = 0.75
GAP = 1 / 128
SAMPLES = 16
TOLERANCE = 3e-9


def _positions(model, level, radius):
    sender_port = model.addressed_seam[0]
    sender = sender_port // PORTS_PER_CELL
    centers = np.asarray(level["cell_centers"])
    origin = centers[sender]
    arc = radius * np.arctan2(np.linalg.norm(np.cross(origin, centers), axis=1), centers @ origin)
    midpoint = [(arc[a // PORTS_PER_CELL] + arc[b // PORTS_PER_CELL]) / 2
                for a, b in model.seam_ids]
    return np.r_[arc[list(model.cell_ids)], midpoint]


def _velocity_observable(operator, positions, cells, speed):
    position = sparse.diags(positions, format="csr")
    velocity = 1j * (operator @ position - position @ operator)
    grading = sparse.diags(np.r_[np.ones(cells), -np.ones(len(positions) - cells)], format="csr")
    anticommutator = grading @ velocity + velocity @ grading
    if np.max(abs(anticommutator.data), initial=0.0) > TOLERANCE:
        raise ValueError("the physically scaled velocity is not odd in the clock grading")
    # For a 1-Lipschitz cell read and seam midpoints, each half-link
    # entry is bounded by c/2. Cell degree is at most 3 and edge degree 2.
    off_diagonal = velocity[:cells, cells:]
    if (np.max(abs(off_diagonal.data), initial=0.0) > speed / 2 + TOLERANCE
            or np.max(np.asarray(abs(off_diagonal).sum(axis=1)), initial=0.0) > 1.5 * speed + TOLERANCE
            or np.max(np.asarray(abs(off_diagonal).sum(axis=0)), initial=0.0) > speed + TOLERANCE):
        raise ValueError("the atlas read violates its half-link Lipschitz certificate")
    speed_bound = speed * math.sqrt(1.5)
    norm = float(np.max(abs(np.linalg.eigvalsh(velocity.toarray())), initial=0.0))
    if norm > speed_bound + TOLERANCE:
        raise ValueError("the velocity exceeds the shared envelope")
    return velocity, norm, speed_bound


def _groups(model):
    sender, receiver = (p // PORTS_PER_CELL for p in model.addressed_seam)
    modes = len(model.cell_ids) + len(model.seam_ids)
    locations = np.full(modes, 3, dtype=int)
    locations[model.cell_ids.index(sender)] = 0
    locations[model.cell_ids.index(receiver)] = 1
    locations[len(model.cell_ids) + model.seam_ids.index(model.addressed_seam)] = 2
    groups = np.arange(2)[:, None, None] * 16 + locations[None, :, None] * 4 + locations[None, None, :]
    return groups.ravel()


def _mode_marginal(state, shape, axis):
    vectors = np.moveaxis(state.reshape(shape), axis, 0).reshape(shape[axis], -1)
    return np.sum(abs(vectors) ** 2, axis=1)


def _spectral_comparison(operator, mass, states, shape, axis):
    energies, vectors = np.linalg.eigh(operator.toarray())
    negatives, shell_weights = energies < 0, mass / abs(energies)
    leakage = []
    shell_rates = []
    for state in states:
        columns = np.moveaxis(state.reshape(shape), axis, 0).reshape(shape[axis], -1)
        weights = np.sum(abs(vectors.conj().T @ columns) ** 2, axis=1)
        leakage.append(float(weights[negatives].sum()))
        shell_rates.append(float(shell_weights @ weights))

    # In the free interaction picture the negative amplitude is G-Lipschitz
    # and the expectation of m/|D| is 2G-Lipschitz. The midpoint enclosure
    # bounds the between-sample error without assuming a smooth trajectory.
    step = WAIT / SAMPLES
    amplitude_upper = np.minimum(1, np.sqrt(leakage) + G * step / 2)
    leakage_bound = float(np.mean(2 * amplitude_upper ** 2 + 2 * amplitude_upper))
    midpoint_error = G * step / 2
    return {
        "negative_band_probabilities_at_midpoints": leakage,
        "mass_shell_rate_midpoint_mean": float(np.mean(shell_rates)),
        "midpoint_mean_error_bound": midpoint_error,
        "clock_to_mass_shell_error_bound": leakage_bound + midpoint_error,
    }


def _joint_replay(model, additive, level, radius, speed):
    groups = _groups(model)
    terminal = next(evolve(model.hamiltonian, model.initial, (WAIT,)))
    no_contact_terminal = next(evolve(additive.hamiltonian, model.initial, (WAIT,)))
    records = np.bincount(groups, weights=abs(terminal) ** 2, minlength=32)
    control = np.bincount(groups, weights=abs(no_contact_terminal) ** 2, minlength=32)
    np.testing.assert_allclose([records.sum(), control.sum()], 1, atol=TOLERANCE, rtol=0)
    tv = float(np.abs(records - control).sum() / 2)
    if tv < 1e-7:
        raise ValueError("the joint total-variation control did not resolve the contact")

    difference = model.hamiltonian - additive.hamiltonian
    difference.eliminate_zeros()
    if difference.nnz != 1 or np.count_nonzero(difference.diagonal()) != 1:
        raise ValueError("the control changed more than the addressed contact cost")
    contact_commutator = additive.hamiltonian @ difference - difference @ additive.hamiltonian
    if np.max(abs(contact_commutator.data), initial=0.0) <= TOLERANCE:
        raise ValueError("the contact commutes with the rest of the source")

    step = WAIT / SAMPLES
    midpoints = tuple((i + 0.5) * step for i in range(SAMPLES))
    states = list(evolve(model.hamiltonian, model.initial, midpoints))
    positions = _positions(model, level, radius)
    if np.any(model.initial.imag):
        raise ValueError("the zero-initial-current certificate needs a real preparation")
    results = {}
    for name, axis, operator, grading, rest in (
        ("charge", 1, model.charge_operator, model.charge_grading, model.charge_rest),
        ("cage", 2, model.cage_operator, model.cage_grading, model.cage_rest),
    ):
        velocity, velocity_norm, speed_bound = _velocity_observable(operator, positions, len(model.cell_ids), speed)
        initial_columns = np.moveaxis(model.initial.reshape(model.shape), axis, 0).reshape(model.shape[axis], -1)
        initial_current = float(np.vdot(initial_columns, velocity @ initial_columns).real)
        initial_grade = float(grading @ (abs(model.initial) ** 2))
        if abs(initial_current) > TOLERANCE or not 0 < initial_grade < 1 - TOLERANCE:
            raise ValueError("the initial clock read did not distinguish the centroid proper-time rule")
        initial_distance = float(positions @ _mode_marginal(model.initial, model.shape, axis))
        final_distance = float(positions @ _mode_marginal(terminal, model.shape, axis))
        drift = (final_distance - initial_distance) / WAIT
        moving_clock = clock.audit(model.hamiltonian, grading, model.initial, WAIT, gap=GAP, groups=groups)
        rest_clock = clock.audit(model.hamiltonian, grading, rest, WAIT, gap=GAP, groups=groups)
        np.testing.assert_allclose(rest_clock["phase_clock_rate"], 1, atol=TOLERANCE, rtol=0)
        additive_clock = clock.audit(additive.hamiltonian, grading, model.initial, WAIT, gap=GAP, groups=groups)
        exact_rate = moving_clock["exact_grade_rate"]
        phase_shift = moving_clock["phase_clock_rate"] - additive_clock["phase_clock_rate"]
        phase_resolution = (
            moving_clock["normalized_phase_error_bound"] + additive_clock["normalized_phase_error_bound"]
        )
        if abs(drift) <= 1e-5 or abs(phase_shift) <= phase_resolution + TOLERANCE:
            raise ValueError("the joint replay did not distinguish transport and a contact clock shift")
        disk = exact_rate ** 2 + (drift / speed_bound) ** 2
        if disk > 1 + TOLERANCE:
            raise ValueError("the co-active history violates its clock–velocity envelope")
        mass = float(operator[0, 0])
        spectral = _spectral_comparison(operator, mass, states, model.shape, axis)
        if (abs(exact_rate - spectral["mass_shell_rate_midpoint_mean"])
                > spectral["clock_to_mass_shell_error_bound"] + TOLERANCE):
            raise ValueError("the integrated clock and its instantaneous spectral comparison disagree")
        results[name] = {
            "initial_and_final_mean_radial_distance": [initial_distance, final_distance],
            "initial_radial_velocity": initial_current,
            "initial_mass_response_rate": initial_grade,
            "initial_centroid_proper_time_rate": 1.0,
            "mean_radial_drift": drift,
            "radial_velocity_operator_norm": velocity_norm,
            "shared_velocity_envelope": speed_bound,
            "clock_velocity_disk_read": disk,
            "motion_clock": moving_clock,
            "rest_clock": rest_clock,
            "rest_normalized_detector_rate": moving_clock["phase_clock_rate"] / rest_clock["phase_clock_rate"],
            "additive_clock": additive_clock,
            "contact_phase_rate_change": phase_shift,
            "contact_phase_rate_change_error_bound": phase_resolution,
            "contact_clock_rate_change": exact_rate - additive_clock["exact_grade_rate"],
            "spectral_comparison": spectral,
        }
    return {
        "seam": list(model.addressed_seam),
        "swap_seam": list(model.swap_seam),
        "native_load_pairs": model.pairs,
        "shape": list(model.shape),
        "interaction_joint_record": records.reshape(2, 4, 4).tolist(),
        "additive_joint_record": control.reshape(2, 4, 4).tolist(),
        "joint_record_total_variation": tv,
        "objects": results,
    }


def run(*, include_finest=True):
    carriers = [carrier.federation(level) for level in range(3)]
    atlas = geometry.run(carriers)
    speed, radius = atlas["calibrated_speed"], atlas["radius"]
    lengths = [{(a, b): d for a, b, d in level["seam_lengths"]} for level in atlas["levels"]]
    patches = []
    for car, links, level in zip(carriers, lengths, atlas["levels"], strict=True):
        active = source.build(car, links, speed, selected_seam_only=True, coupling=G)
        additive = source.build(car, links, speed, selected_seam_only=True, coupling=0)
        patch = _joint_replay(active, additive, level, radius, speed)
        patches.append({"cells": car.cells, "positive_band_leakage": source.leakage_witness(active), **patch})

    full_levels = len(carriers) if include_finest else 2
    whole = []
    for i in range(full_levels):
        car, links, level = carriers[i], lengths[i], atlas["levels"][i]
        active = source.build(car, links, speed, coupling=G)
        additive = source.build(car, links, speed, coupling=0)
        whole.append({"cells": car.cells, **_joint_replay(active, additive, level, radius, speed)})
    return {
        "atlas_and_calibration": atlas,
        "contact_strength": G,
        "free_wait": WAIT,
        "mass_tag_gap": GAP,
        "position_regions": ["sender", "receiver", "addressed_seam", "elsewhere"],
        "addressed_patches": patches,
        "all_links": whole,
    }
