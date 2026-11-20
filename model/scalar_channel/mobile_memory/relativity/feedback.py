"""Use a measured signed half-link current to enable the incidence cage link."""

import numpy as np

from ....contrast_memory_spherical_v1 import carrier
from ...radiative_repair.source import evolve
from ..current import one_charge_read
from ..mobility import PORTS_PER_CELL, home_seam
from .incidence import (
    CAGE_MASS, CHARGE_MASS, addressed_records, cell_incidence, dirac_operator, positive_embedding,
)


SIGNS = (-1, 1)
CURRENTS = (-1, 0, 1)
PRE_READ_WAIT = 2
TOLERANCE = 2e-10


def run():
    car = carrier.federation(0)
    seam = home_seam(car)
    sender, receiver = (port // PORTS_PER_CELL for port in seam)
    incidence, seams = cell_incidence(car)
    edge_mode = car.cells + seams.index(seam)
    source = dirac_operator(incidence, CHARGE_MASS)
    hop = float(source[sender, edge_mode])
    if hop <= 0:
        raise ValueError("the declared half-link orientation changed")
    positive = positive_embedding(incidence, CHARGE_MASS)
    waited = next(evolve(source, positive[:, sender], (PRE_READ_WAIT,)))

    scalar = addressed_records(car, CAGE_MASS, seam)
    blocked, written = np.asarray(scalar["blocked"][-1]), np.asarray(scalar["written"][-1])
    wait = scalar["waits"][-1]
    receiver_law = np.stack((blocked, blocked, written))
    reports = []
    for charge in SIGNS:
        # The native difference-hop reader has a negative off-diagonal.
        # A positive incidence half-hop reverses its oriented sign.
        branches = one_charge_read(waited, (sender, edge_mode), -charge)
        branch_probabilities = np.array([np.vdot(branches[r], branches[r]).real for r in CURRENTS])
        local_current = -2 * charge * hop * np.imag(waited[sender].conj() * waited[edge_mode])
        np.testing.assert_allclose(branch_probabilities @ CURRENTS, local_current / hop,
                                   atol=TOLERANCE, rtol=0)

        admitted = np.column_stack([positive @ (positive.conj().T @ branches[r]) for r in CURRENTS])
        charge_at_read = abs(next(evolve(source, admitted, (wait,))).T) ** 2
        other_charge_sites = np.ones(charge_at_read.shape[1], dtype=bool)
        other_charge_sites[[sender, receiver]] = False
        charge_record = np.column_stack((
            charge_at_read[:, sender], charge_at_read[:, receiver],
            charge_at_read[:, other_charge_sites].sum(axis=1),
        ))
        rejections = branch_probabilities - charge_record.sum(axis=1)
        if rejections.min() < -TOLERANCE:
            raise ValueError("the charge filter gained probability")
        rejections = np.maximum(rejections, 0.0)

        feedback = charge_record[:, :, None] * receiver_law[:, None, :]
        enable_rate = charge_record[-1].sum() / charge_record.sum()
        shuffled = charge_record[:, :, None] * ((1 - enable_rate) * blocked + enable_rate * written)
        np.testing.assert_allclose(feedback.sum(axis=(0, 1)), shuffled.sum(axis=(0, 1)),
                                   atol=TOLERANCE, rtol=0)
        np.testing.assert_allclose(feedback.sum(axis=2), shuffled.sum(axis=2),
                                   atol=TOLERANCE, rtol=0)
        np.testing.assert_allclose(feedback.sum(axis=(1, 2)) + rejections, branch_probabilities,
                                   atol=TOLERANCE, rtol=0)
        distance = float(abs(feedback - shuffled).sum() / 2)
        if distance < 1e-4:
            raise ValueError("the incidence cage read did not distinguish feedback from shuffling")
        reports.append({
            "charge": charge,
            "measured_current_probabilities": branch_probabilities.tolist(),
            "positive_band_rejections_by_current": rejections.tolist(),
            "feedback": feedback.tolist(),
            "shuffled": shuffled.tolist(),
            "total_variation": distance,
        })
    return {
        "seam": list(seam),
        "oriented_half_link_modes": [sender, edge_mode],
        "pre_read_wait": PRE_READ_WAIT,
        "post_read_wait": wait,
        "current_outcomes": list(CURRENTS),
        "charge_regions": ["sender", "receiver", "elsewhere"],
        "cage_mode_order": {"cells": list(range(car.cells)), "seams": list(seams)},
        "written_cage_receiver_gain": scalar["written_receiver_gain"],
        "records": reports,
    }
