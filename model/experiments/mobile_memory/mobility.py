"""One-charge detector records from a preparation in the lowest spectral band."""

import numpy as np

from model.executable_model import carrier
from model.candidate_laws.scalar_channel.carrier import require
from .law import EXTERNAL_WEIGHT, MASS, charged_operator


WAITS = (0, 8, 24, 48, 96)
PORTS_PER_CELL = 12
HOME_CELL = 0
CELL_GAP = 5 - np.sqrt(5)
TOLERANCE = 1e-10


def _detector_records(energies, vectors, initial, reject):
    spectral = vectors.conj().T @ initial
    # A global phase has no effect on the detector readout.
    amplitudes = vectors @ (
        spectral[:, None] * np.exp(-1j * (energies[:, None] - MASS) * WAITS)
    )
    cell_probabilities = (np.abs(amplitudes) ** 2).reshape(-1, PORTS_PER_CELL, len(WAITS))
    records = cell_probabilities.sum(axis=1).T
    require(np.max(np.abs(records.sum(axis=1) + reject - 1)) < TOLERANCE,
            "detector records lost probability")
    return records


def _filtered_home(car, energies, vectors):
    require(MASS + 2 * EXTERNAL_WEIGHT < MASS + CELL_GAP, "the band intervals overlap")
    require(abs(energies[0] - MASS) < TOLERANCE
            and energies[car.cells - 1] <= MASS + 2 * EXTERNAL_WEIGHT + TOLERANCE
            and energies[car.cells] >= MASS + CELL_GAP - TOLERANCE,
            "the computed spectrum violates the matching bound")
    home_orbital = (np.arange(car.n) // PORTS_PER_CELL == HOME_CELL) / np.sqrt(PORTS_PER_CELL)
    low_vectors = vectors[:, :car.cells]
    return low_vectors @ (low_vectors.conj().T @ home_orbital)


def lowest_band_preparation(car):
    """Project the uniform home-cell orbital, keeping its filter probability in the norm."""
    energies, vectors = np.linalg.eigh(charged_operator(car, EXTERNAL_WEIGHT).toarray())
    return _filtered_home(car, energies, vectors)


def home_seam(car):
    """Choose the oriented seam from cell zero to its first neighboring cell."""
    return min(
        ((a, b) for a, b in car.external if a // PORTS_PER_CELL == HOME_CELL),
        key=lambda seam: seam[1],
    )


def run():
    """Return joint filter and cell-read probabilities for every registered wait."""
    car = carrier.federation(0)
    geometry = carrier.federation_checks(car)
    require(geometry["every_cell_complete"] and geometry["external_seams_form_a_matching"],
            "spectral bound requires complete cells and an external matching")
    source = HOME_CELL
    receiver = home_seam(car)[1] // PORTS_PER_CELL
    active = charged_operator(car, EXTERNAL_WEIGHT).toarray()
    parked = charged_operator(car, EXTERNAL_WEIGHT, parked=True).toarray()
    cells = np.arange(car.n) // PORTS_PER_CELL
    within_cells = cells[:, None] == cells
    require(np.array_equal(active[within_cells], parked[within_cells]),
            "parked law changed the within-cell blocks")
    require(not np.any(parked[~within_cells]), "parked law transfers between cells")

    energies, vectors = np.linalg.eigh(active)
    band = energies[:car.cells]
    floor = MASS
    low_ceiling = MASS + 2 * EXTERNAL_WEIGHT
    upper_floor = MASS + CELL_GAP
    require(band[1] > band[0] + TOLERANCE, "the lowest eigenvalue is not simple")

    accepted = _filtered_home(car, energies, vectors)
    acceptance = float(np.vdot(accepted, accepted).real)
    minimum_acceptance = 1 - (2 * EXTERNAL_WEIGHT / CELL_GAP) ** 2
    require(minimum_acceptance - TOLERANCE <= acceptance <= 1,
            "the filter acceptance violates the high-band norm bound")
    reject = 1 - acceptance

    active_records = _detector_records(energies, vectors, accepted, reject)
    parked_energies, parked_vectors = np.linalg.eigh(parked)
    parked_records = _detector_records(parked_energies, parked_vectors, accepted, reject)
    require(np.max(np.abs(active_records[0] - parked_records[0])) < TOLERANCE,
            "the two laws have different initial detector records")
    parked_change = float(np.max(np.abs(parked_records - parked_records[0])))
    active_gain = float(active_records[:, receiver].max() - active_records[0, receiver])
    require(parked_change < TOLERANCE and active_gain > 0.02,
            "the receiver did not distinguish active and parked transport")

    return {
        "carrier_cells": car.cells,
        "source_cell": source,
        "receiver_cell": receiver,
        "external_weight": EXTERNAL_WEIGHT,
        "band": {
            "dimension": car.cells,
            "allowed_interval": [floor, low_ceiling],
            "next_level_lower_bound": float(upper_floor),
            "lowest": float(band[0]),
            "second": float(band[1]),
            "highest": float(band[-1]),
            "next_level": float(energies[car.cells]),
            "filter_acceptance_lower_bound": float(minimum_acceptance),
        },
        "filter_acceptance_probability": acceptance,
        "filter_rejection_probability": reject,
        "waits": list(WAITS),
        "active_cell_probabilities": active_records.tolist(),
        "parked_cell_probabilities": parked_records.tolist(),
        "receiver_gain": active_gain,
        "maximum_parked_change": parked_change,
    }
