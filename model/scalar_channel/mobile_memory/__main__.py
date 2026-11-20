"""Replay the mobility and memory-controlled records of the finite candidate."""

from fractions import Fraction
import json

import numpy as np

from ...contrast_memory_spherical_v1 import carrier
from ..radiative_repair.source import evolve
from . import coactive, mobility, refinement
from .cage_transfer import transfer
from .contrast_binding import run as contrast_cage
from .current import current_experiment, one_charge_read
from .law import charged_operator
from .memory import coherent_check
from .register_link import compile_link


CHARGE_WAIT = 24
CURRENT_OUTCOMES = (-1, 0, 1)
SIGNS = (-1, 1)


def current_controlled_cage(car, seam):
    """Feed the actual signed read of a low-band charge into the cage link."""
    accepted = mobility.lowest_band_preparation(car)
    acceptance = float(np.vdot(accepted, accepted).real)
    waited = next(evolve(charged_operator(car), accepted, (CHARGE_WAIT,)))
    enabled = transfer(seam=seam, memory=True)
    disabled = transfer(seam=seam, memory=False)
    active_fuel = np.array(enabled["receiver_record"])
    blank_fuel = np.array(disabled["receiver_record"])
    rows = []
    for charge in SIGNS:
        branches = one_charge_read(waited, seam, charge)
        current = np.array([np.vdot(branches[r], branches[r]).real for r in CURRENT_OUTCOMES])
        feedback = current[:, None] * np.stack((blank_fuel, blank_fuel, active_fuel))
        blocked = current[:, None] * blank_fuel
        enable_rate = current[-1] / acceptance
        shuffled = current[:, None] * (enable_rate * active_fuel + (1 - enable_rate) * blank_fuel)
        if any(abs(record.sum() - acceptance) > 1e-10 for record in (feedback, blocked, shuffled)):
            raise ValueError("the joint cage detector lost accepted trials")
        np.testing.assert_allclose(feedback.sum(axis=0), shuffled.sum(axis=0), atol=1e-10, rtol=0)
        np.testing.assert_allclose(feedback.sum(axis=1), shuffled.sum(axis=1), atol=1e-10, rtol=0)
        distance = float(np.abs(feedback - shuffled).sum() / 2)
        if distance < 1e-4:
            raise ValueError("the stored current did not distinguish the joint detector records")
        rows.append({
            "charge": charge,
            "feedback": feedback.tolist(),
            "blocked": blocked.tolist(),
            "shuffled": shuffled.tolist(),
            "feedback_shuffled_total_variation": distance,
        })
    if active_fuel[3] < blank_fuel[3] + 0.2:
        raise ValueError("enabling the cage link did not shift the receiver fuel read")
    return {
        "charge_wait": CHARGE_WAIT,
        "filter_rejection_probability": 1 - acceptance,
        "current_records": list(CURRENT_OUTCOMES),
        "receiver_records": ["00", "01", "10", "11"],
        "joint_probabilities": rows,
        "cage_enabled": enabled,
        "cage_blocked": disabled,
    }


def run():
    car = carrier.federation(0)
    seam = mobility.home_seam(car)
    pair = current_experiment(car, seam)
    for rows in pair["joint_probabilities"].values():
        np.testing.assert_allclose(np.sum(rows, axis=(1, 2)), (1, 1), atol=1e-10, rtol=0)
    if pair["feedback_advantage"] <= float(Fraction(pair["feedback_advantage_lower_bound"])):
        raise ValueError("the neutral receiver did not distinguish the feedback pulse")
    return {
        "mobility": mobility.run(),
        "pair_current": pair,
        "current_controlled_cage": current_controlled_cage(car, seam),
        "contrast_cage": contrast_cage(seam=seam),
        "coactive": coactive.run(),
        "refinement": refinement.run(),
        "register_link": compile_link(seam=seam),
        "coherent_memory_errors": coherent_check(),
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
