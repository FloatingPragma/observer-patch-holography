"""Joint charge, load-pair and cage records with both seam controls active."""

from fractions import Fraction
import json

import numpy as np

from model.executable_model import carrier
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.experiments.radiative_repair.source import evolve
from .cage_transfer import PORTS_PER_CELL
from .contrast_binding import (
    _activation, _addressed_pair, _all_links_probabilities, _all_links_record, _external_laplacian,
    _native_controls,
)
from .current import one_charge_read
from .law import charged_operator
from .mobility import home_seam, lowest_band_preparation


CHARGE_SIGNS = (-1, 1)
CURRENT_OUTCOMES = (-1, 0, 1)
CHARGE_WAIT = 24
POST_WRITE_WAIT = 0.5
KAPPA = 128
COUPLING = 1
CERTIFICATE_WAIT = Fraction(1, 128)


def _joint_arms(charge_record, load_cage_records):
    """Combine the exact conditional factors for one prepared charge sign."""
    pointer_on = np.array(CURRENT_OUTCOMES) > 0
    feedback = charge_record[:, :, None, None] * load_cage_records[pointer_on.astype(int), None]
    p_on = charge_record[pointer_on].sum() / charge_record.sum()
    mixed_load_cage = (1 - p_on) * load_cage_records[0] + p_on * load_cage_records[1]
    shuffled = charge_record[:, :, None, None] * mixed_load_cage
    return feedback, shuffled


def _short_wait_certificate(sector, load_source, contacts, initial_load, seam):
    """Certify a pointer effect in the joint original-pair and receiver outcome."""
    pair = int(np.flatnonzero(initial_load).item())
    sender, receiver = (port // PORTS_PER_CELL for port in seam)
    a, b = sector.aa[pair], sector.bb[pair]
    swap_bound = int(sector.car.degree[a] + sector.car.degree[b] - 2)
    column = load_source.getcol(pair).tocoo()
    moves = column.data[column.row != pair]
    if not (len(moves) <= swap_bound <= 10 and np.all(abs(moves) == 0.5)):
        raise ValueError("the initial pair's outgoing swaps exceed the local degree bound")
    off_diagonal_squared = Fraction(len(moves), 4)
    source_idx = pair * sector.car.cells + sender
    target_idx = pair * sector.car.cells + receiver
    slopes = [-int(contact[target_idx, source_idx]) * COUPLING for contact in contacts]
    if slopes != [COUPLING, 2 * COUPLING]:
        raise ValueError("the initial pair and positive current do not give the required link weights")

    # On the original-load subspace the cage block, after the common load
    # phase, has norm at most 6g. Two integrations bound excursions off and
    # back to that load subspace by ||B||² t² / 2. Both contacts are load diagonal.
    t = CERTIFICATE_WAIT
    remainder = (18 * COUPLING ** 2 + off_diagonal_squared / 2) * t ** 2
    on_lower = max(Fraction(0), slopes[1] * t - remainder) ** 2
    off_upper = (slopes[0] * t + remainder) ** 2
    margin = on_lower - off_upper
    if margin <= 0:
        raise ValueError("the short wait does not certify a positive pointer effect")
    probabilities = [
        float(_all_links_probabilities(
            sector, load_source, contact, initial_load, sender, COUPLING, float(t)
        )[pair, receiver])
        for contact in contacts
    ]
    if probabilities[1] - probabilities[0] < float(margin) - 1e-12:
        raise ValueError("the computed original-pair receiver record violates the certificate")
    return {
        "wait": str(t),
        "receiver_cell": receiver,
        "initial_pair_port_degrees": [int(sector.car.degree[a]), int(sector.car.degree[b])],
        "outgoing_native_swaps": len(moves),
        "outgoing_native_swap_bound": swap_bound,
        "off_diagonal_norm_squared": str(off_diagonal_squared),
        "addressed_amplitude_slopes_without_and_with_pointer": slopes,
        "amplitude_remainder_bound": str(remainder),
        "probabilities_without_and_with_pointer": probabilities,
        "probability_increase_lower_bound": str(margin),
    }


def run():
    """Keep every outcome of a signed read followed by the common free generator."""
    car = carrier.federation(0)
    hypotheses = carrier.federation_checks(car)
    if not hypotheses["every_cell_complete"] or not hypotheses["external_seams_form_a_matching"]:
        raise ValueError("the binding bound requires complete cells and a port matching")
    addressed = home_seam(car)
    sender = addressed[0] // PORTS_PER_CELL
    sector = PairSector(car)
    blank_weights, _ = _native_controls(KAPPA)
    initial_pair, initial_load, _ = _addressed_pair(sector, addressed, blank_weights)
    contrasts = _external_laplacian(
        sector.P, car.cells,
        [(seam, _activation(sector, seam, blank_weights)) for seam in car.external],
    )
    pointer = _external_laplacian(
        sector.P, car.cells, [(addressed, np.ones(sector.P, dtype=int))],
    )
    contact_on = contrasts + pointer
    if float(abs(contact_on).sum(axis=1).max()) > 6:
        raise ValueError("the co-active contact exceeds the binding proof's norm bound")
    load_source = sector.hamiltonian(KAPPA)
    contacts = (contrasts, contact_on)
    load_cage_records = np.array([
        _all_links_record(sector, load_source, contact, initial_load, sender, COUPLING, POST_WRITE_WAIT)
        for contact in contacts
    ])
    separation_bound = 40 / (2 * KAPPA - 6 - 6 * COUPLING) ** 2
    if np.any(load_cage_records[:, 0].sum(axis=1) > separation_bound + 1e-10):
        raise ValueError("the co-active source exceeds the all-time separation bound")

    # The actual signed-current read writes a stationary label. For each
    # label the subsequent charge propagator commutes with load and cage.
    # The unnormalized branches keep the pre-read spectral filter's failures.
    charge_source = charged_operator(car)
    accepted = lowest_band_preparation(car)
    acceptance = float(np.vdot(accepted, accepted).real)
    waited = next(evolve(charge_source, accepted, (CHARGE_WAIT,)))
    projected = [one_charge_read(waited, addressed, sign) for sign in CHARGE_SIGNS]
    columns = np.column_stack([branch[r] for branch in projected for r in CURRENT_OUTCOMES])
    terminal = next(evolve(charge_source, columns, (POST_WRITE_WAIT,)))
    charge_records = (abs(terminal.T) ** 2).reshape(2, 3, car.cells, PORTS_PER_CELL).sum(axis=-1)
    arms = [_joint_arms(record, load_cage_records) for record in charge_records]
    feedback, shuffled = (np.array(records) for records in zip(*arms, strict=True))

    for table in (feedback, shuffled):
        np.testing.assert_allclose(table.sum(axis=(1, 2, 3, 4)), acceptance, atol=1e-9, rtol=0)
    np.testing.assert_allclose(feedback.sum(axis=(3, 4)), shuffled.sum(axis=(3, 4)), atol=1e-9, rtol=0)
    np.testing.assert_allclose(feedback.sum(axis=(1, 2)), shuffled.sum(axis=(1, 2)), atol=1e-9, rtol=0)
    distances = abs(feedback - shuffled).sum(axis=(1, 2, 3, 4)) / 2
    if np.any(distances < 1e-4):
        raise ValueError("the co-active joint records did not distinguish the actual write from shuffling")

    return {
        "seam": list(addressed),
        "initial_pair": initial_pair,
        "charge_wait": CHARGE_WAIT,
        "post_write_wait": POST_WRITE_WAIT,
        "filter_rejection_probability": 1 - acceptance,
        "record_axes": ["charge_sign", "actual_current", "charge_cell", "pair_presence", "cage_cell"],
        "charge_signs": list(CHARGE_SIGNS),
        "actual_currents": list(CURRENT_OUTCOMES),
        "pair_presence": ["no_adjacent_pair", "adjacent_pair"],
        "cells": list(range(car.cells)),
        "feedback": feedback.tolist(),
        "shuffled": shuffled.tolist(),
        "joint_total_variation_by_sign": distances.tolist(),
        "all_time_no_pair_bound": separation_bound,
        "conditional_load_cage_without_and_with_pointer": load_cage_records.tolist(),
        "short_wait_certificate": _short_wait_certificate(sector, load_source, contacts, initial_load, addressed),
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
