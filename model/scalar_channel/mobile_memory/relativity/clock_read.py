"""Two successive clock reads with a definite initial internal label."""

import math

import numpy as np

from ....contrast_memory_spherical_v1 import carrier
from ...radiative_repair.source import evolve
from .incidence import (
    CAGE_MASS, CHARGE_MASS, SPEED, cell_incidence, dirac_operator,
    laplacian_modes, positive_mode,
)


SIGNS = (-1, 1)
ANALYZERS = (1, -1j)
GAP = 1 / 8
CLOCK_WAITS = (2, 4, 8)
SPECTRAL_WAIT = 1 / 4
TOLERANCE = 2e-10


def two_read_records(lower_arm, upper_arm, observed_cell):
    """Accumulate both pointer reads and the complementary position outcomes.

    Each arm includes amplitude 1/2 from the first X read of a definite label.
    The positive first outcome is represented by the given arms. The negative
    first outcome reverses the sign of the second arm.
    """
    complement = np.ones(len(lower_arm), dtype=bool)
    complement[observed_cell] = False
    accepted_one = np.vdot(lower_arm, lower_arm).real + np.vdot(upper_arm, upper_arm).real
    if accepted_one > 0.5 + TOLERANCE:
        raise ValueError("the filter probabilities exceed the first-read probability")
    rejected_first = max(0.0, float(0.5 - accepted_one))
    joint = np.empty((2, 2, 2, 2))
    for axis, phase in enumerate(ANALYZERS):
        for i, first in enumerate(SIGNS):
            for j, terminal in enumerate(SIGNS):
                occupation = abs(lower_arm + first * terminal * phase * upper_arm) ** 2 / 2
                joint[axis, i, j] = (occupation[observed_cell], occupation[complement].sum())
    np.testing.assert_allclose(joint.sum(axis=(1, 2, 3)) + 2 * rejected_first, 1, atol=TOLERANCE, rtol=0)

    parity = np.outer(SIGNS, SIGNS)
    real, imag = (joint.sum(axis=-1) * parity).sum(axis=(1, 2))
    visibility = math.hypot(real, imag)
    if visibility < TOLERANCE:
        raise ValueError("the clock record does not have a measurable phase")
    return {
        "joint_records": joint.tolist(),
        "rejected_by_first_read": [rejected_first, rejected_first],
        "visibility": visibility,
        "phase": -math.atan2(imag, real),
    }


def _phase_rows(lower, upper, waits, read_cell):
    return [
        {"wait": t, **two_read_records(left, right, read_cell)}
        for t, left, right in zip(waits, lower, upper, strict=True)
    ]


def mass_clock(B, mass, gap, eigenvalue, mode, read_cell):
    """Read the centered mass splitting after the first measured event.

    The same spatial state initially enters both arms of the first read.
    A heralded filter admits their respective positive bands in this mode.
    Rejections remain in every joint record.
    """
    if not 0 <= gap < 2 * mass:
        raise ValueError("both clock masses must be positive")
    low, high = mass - gap / 2, mass + gap / 2
    E_low, seed = positive_mode(B, low, eigenvalue, mode)
    H_low, H_high = dirac_operator(B, low), dirac_operator(B, high)
    E_high = math.hypot(high, SPEED * math.sqrt(eigenvalue))
    # On the declared incidence doublet this polynomial is the positive projector.
    lower = seed / 2
    upper = (seed + (H_high @ seed) / E_high) / 4
    np.testing.assert_allclose(H_high @ upper, E_high * upper, atol=TOLERANCE, rtol=0)
    rows = _phase_rows(evolve(H_low, lower, CLOCK_WAITS),
                       evolve(H_high, upper, CLOCK_WAITS), CLOCK_WAITS, read_cell)
    observed = [row["phase"] / row["wait"] for row in rows]
    # All enrolled waits give |phase| ≤ 1 for the declared gap.
    expected = gap * (2 * mass) / (E_high + E_low)
    np.testing.assert_allclose(observed, expected, atol=TOLERANCE, rtol=0)
    return {"rows": rows, "observed_frequency": observed[0]}


def spectral_read(B, mass, eigenvalue, mode, read_cell):
    """Measure one line against a parked reference using the same pointer counts."""
    energy, seed = positive_mode(B, mass, eigenvalue, mode)
    active = next(evolve(dirac_operator(B, mass), seed / 2, (SPECTRAL_WAIT,)))
    row = two_read_records(seed / 2, active, read_cell)
    observed = row["phase"] / SPECTRAL_WAIT
    np.testing.assert_allclose(observed, energy, atol=TOLERANCE, rtol=0)
    return {"wait": SPECTRAL_WAIT, "observed_energy": observed, **row}


def run():
    """Read one nonzero mode and its rest control for both transported objects."""
    car = carrier.federation(0)
    B, _ = cell_incidence(car)
    levels, modes = laplacian_modes(B)
    index = int(np.flatnonzero(levels > 1e-8)[-1])
    eigenvalue = float(levels[index])
    rest_mode = np.ones(car.cells) / math.sqrt(car.cells)
    moving_mode = modes[:, index]
    read_cell = int(np.argmax(abs(moving_mode)))
    records = {}
    for name, mass in (("charge", CHARGE_MASS), ("cage", CAGE_MASS)):
        rest = mass_clock(B, mass, GAP, 0.0, rest_mode, read_cell)
        moving = mass_clock(B, mass, GAP, eigenvalue, moving_mode, read_cell)
        equal_labels = mass_clock(B, mass, 0.0, eigenvalue, moving_mode, read_cell)
        rest_line = spectral_read(B, mass, 0.0, rest_mode, read_cell)
        moving_line = spectral_read(B, mass, eigenvalue, moving_mode, read_cell)
        c_squared = (moving_line["observed_energy"] ** 2 - rest_line["observed_energy"] ** 2) / eigenvalue
        rate = moving["observed_frequency"] / rest["observed_frequency"]

        energy, state = positive_mode(B, mass, eigenvalue, moving_mode)
        grading_counts = [float(np.vdot(part, part).real) for part in (state[:car.cells], state[car.cells:])]
        grade_rate = grading_counts[0] - grading_counts[1]
        velocity = SPEED ** 2 * math.sqrt(eigenvalue) / energy
        finite_gap_error = GAP ** 2 / (8 * (mass - GAP / 2) ** 2)
        np.testing.assert_allclose(c_squared, SPEED ** 2, atol=TOLERANCE, rtol=0)
        np.testing.assert_allclose(grade_rate, mass / energy, atol=TOLERANCE, rtol=0)
        np.testing.assert_allclose(grade_rate ** 2 + (velocity / SPEED) ** 2, 1, atol=TOLERANCE, rtol=0)
        if not 0 < rate < 1 or abs(rate - grade_rate) > finite_gap_error + TOLERANCE:
            raise ValueError("the clock record failed its finite-gap proper-time rate bound")
        records[name] = {
            "rest_clock": rest,
            "moving_mode_clock": moving,
            "equal_mass_control": equal_labels,
            "rest_line_read": rest_line,
            "moving_line_read": moving_line,
            "squared_speed_from_records": c_squared,
            "rest_normalized_clock_rate": rate,
            "cell_and_edge_occupation": grading_counts,
            "infinitesimal_rate_from_occupation": grade_rate,
            "group_rate_in_declared_momentum": velocity,
            "finite_gap_rate_error_bound": finite_gap_error,
        }
    return {
        "carrier_cells": car.cells,
        "read_cell": read_cell,
        "selected_graph_eigenvalue": eigenvalue,
        "declared_momentum": math.sqrt(eigenvalue),
        "mass_gap": GAP,
        "first_and_last_pointer_signs": list(SIGNS),
        "analyzers": ["X", "Y"],
        "detector_regions": ["addressed_cell", "elsewhere"],
        "objects": records,
    }
