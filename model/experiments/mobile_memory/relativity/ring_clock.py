"""Compare internal clock counts with measured circulation on the same carrier loop."""

import math

import numpy as np

from model.executable_model import carrier
from ..current import one_charge_read
from ..mobility import PORTS_PER_CELL
from .clock_read import GAP, mass_clock
from .incidence import CAGE_MASS, CHARGE_MASS, SPEED, dirac_operator, positive_mode, ring_incidence


SPACING = 1
TOLERANCE = 2e-10


def _half_link_record(source, state, left, right):
    """Read the local continuity operator, including its zero eigenspace."""
    hop = float(source[left, right])
    if hop == 0:
        raise ValueError("a ring half-link has no hopping term")
    branches = one_charge_read(state, (left, right), -1 if hop > 0 else 1)
    probabilities = [float(np.vdot(arm, arm).real) for arm in branches.values()]
    np.testing.assert_allclose(sum(probabilities), 1, atol=TOLERANCE, rtol=0)
    return probabilities, abs(hop) * (probabilities[2] - probabilities[0])


def run():
    """Measure both mass clocks and their circulation at the same Fourier label."""
    car = carrier.federation(0)
    ring, incidence, seams = ring_incidence(car, SPACING)
    phase_step = math.tau / len(ring)
    momentum = phase_step / SPACING
    eigenvalue = (2 * math.sin(phase_step / 2) / SPACING) ** 2
    traveling = np.zeros(car.cells, complex)
    traveling[list(ring)] = np.exp(1j * phase_step * np.arange(len(ring))) / math.sqrt(len(ring))
    rest = np.zeros(car.cells)
    rest[list(ring)] = 1 / math.sqrt(len(ring))
    seam_modes = {tuple(sorted((left // PORTS_PER_CELL, right // PORTS_PER_CELL))): car.cells + i
                  for i, (left, right) in enumerate(seams)}
    read_cell = ring[0]
    objects = {}
    for name, mass in (("charge", CHARGE_MASS), ("cage", CAGE_MASS)):
        operator = dirac_operator(incidence, mass)
        energy, state = positive_mode(incidence, mass, eigenvalue, traveling)
        half_links = []
        for a, b in zip(ring, ring[1:] + ring[:1], strict=True):
            seam = seam_modes[tuple(sorted((a, b)))]
            for left, right in ((a, seam), (seam, b)):
                probabilities, flow = _half_link_record(operator, state, left, right)
                half_links.append({"oriented_modes": [left, right], "probabilities": probabilities, "current": flow})
        # Each positive current here is along the declared loop, at half the cell spacing.
        measured_velocity = SPACING * sum(record["current"] for record in half_links) / 2
        expected_velocity = SPEED ** 2 * math.sin(phase_step) / (SPACING * energy)
        np.testing.assert_allclose(measured_velocity, expected_velocity, atol=TOLERANCE, rtol=0)
        if not 0 < measured_velocity <= SPEED:
            raise ValueError("the loop current failed the shared velocity envelope")

        reference_clock = mass_clock(incidence, mass, GAP, 0.0, rest, read_cell)
        moving_clock = mass_clock(incidence, mass, GAP, eigenvalue, traveling, read_cell)
        rate = moving_clock["observed_frequency"] / reference_clock["observed_frequency"]
        cell_occupation = float(np.vdot(state[:car.cells], state[:car.cells]).real)
        grading_rate = 2 * cell_occupation - 1
        measured_defect = 1 - grading_rate ** 2 - (measured_velocity / SPEED) ** 2
        theory_defect = 4 * SPEED ** 2 * math.sin(phase_step / 2) ** 4 / (SPACING * energy) ** 2
        np.testing.assert_allclose(measured_defect, theory_defect, atol=TOLERANCE, rtol=0)
        if not 0 < measured_defect <= (SPACING * momentum / 2) ** 2:
            raise ValueError("the measured finite-spacing defect is outside its bound")

        gap_bound = GAP ** 2 / (8 * (mass - GAP / 2) ** 2)
        if abs(rate - grading_rate) > gap_bound + TOLERANCE:
            raise ValueError("the measured loop clock violates its finite-gap bound")
        objects[name] = {
            "half_link_records": half_links,
            "measured_circulation_velocity": measured_velocity,
            "moving_clock": moving_clock,
            "rest_clock": reference_clock,
            "rest_normalized_rate": rate,
            "rate_from_cell_and_seam_occupation": grading_rate,
            "finite_spacing_rate_squared_defect": measured_defect,
            "finite_gap_rate_error_bound": gap_bound,
        }
    return {
        "ordered_carrier_cells": list(ring),
        "circulation_current_outcomes": [-1, 0, 1],
        "supplied_cell_spacing": SPACING,
        "fourier_momentum": momentum,
        "objects": objects,
    }
