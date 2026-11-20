"""Finite graph mass-shell and receiver audits for the unchanged source."""

from fractions import Fraction

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from ....contrast_memory_spherical_v1 import carrier
from ...cage import construction as cage
from ...carrier import require
from ..cage_transfer import PORTS_PER_CELL
from ..contrast_binding import _covered_cells
from ..law import EXTERNAL_WEIGHT, MASS, charged_operator
from ..refinement import _cell_seams, _laplacian, _neighbors


WAIT = Fraction(1, 128)
CAGE_REST = 2 * cage.ENERGY
CAGE_REST_SQUARED = 24
SPECTRAL_AUDIT_SPEED_SQUARED = 2 * CAGE_REST
TOLERANCE = 1e-12


def _check_uniform_compression(federation, laplacian):
    ports = np.arange(federation.n)
    cells = ports // PORTS_PER_CELL
    lift = sparse.csr_matrix(
        (np.ones(federation.n, dtype=int), (ports, cells)),
        shape=(federation.n, federation.cells),
    )
    charge_kinetic = charged_operator(federation) - MASS * sparse.eye(
        federation.n, format="csr"
    )
    # The unnormalized lift removes the square roots of twelve. All
    # coefficients are dyadic for the declared charged operator.
    defect = lift.T @ charge_kinetic @ lift - EXTERNAL_WEIGHT * laplacian
    require(not np.any(defect.data), "the exact uniform charged compression changed")


def _receiver_bounds(two_step_paths):
    # Both generators have norm at most six after subtracting the rest phase.
    # Their cubic Taylor remainders are therefore at most 6^3 t^3 / 6.
    # Subordination bounds the comparator's off-diagonal entry below -a/36.
    remainder = 36 * WAIT ** 3
    native_amplitude_upper = two_step_paths * WAIT ** 2 / 2 + remainder
    audit_amplitude_lower = two_step_paths * WAIT / 36 - remainder
    require(audit_amplitude_lower > native_amplitude_upper,
            "the short-wait receiver certificate lost its positive margin")
    return native_amplitude_upper ** 2, audit_amplitude_lower ** 2


def _cell_probabilities(energies, modes, sender):
    # A common rest phase drops out of the complete occupation read.
    amplitudes = modes @ (np.exp(-1j * float(WAIT) * energies) * modes[sender].conj())
    probabilities = abs(amplitudes) ** 2
    require(abs(float(probabilities.sum()) - 1) < TOLERANCE,
            "the cell detector probabilities lost normalization")
    return probabilities


def _federation_record(level):
    federation = carrier.federation(level)
    geometry = carrier.federation_checks(federation)
    require(geometry["every_cell_complete"] and geometry["external_seams_form_a_matching"],
            "the spectral comparison requires complete cells and a port matching")
    _covered_cells(federation, level)
    neighbors = _neighbors(federation.cells, _cell_seams(federation))
    laplacian = _laplacian(neighbors)
    require(connected_components(laplacian, directed=False, return_labels=False) == 1,
            "the cell graph must be connected")
    _check_uniform_compression(federation, laplacian)

    sender = 0
    receiver = min(
        cell
        for neighbor in neighbors[sender]
        for cell in neighbors[neighbor] - neighbors[sender] - {sender}
    )
    laplacian_entry = int(laplacian[receiver, sender])
    two_step_paths = int((laplacian @ laplacian)[receiver, sender])
    require(laplacian_entry == 0 and two_step_paths == len(neighbors[receiver] & neighbors[sender])
            and two_step_paths > 0,
            "the chosen detector is not a graph distance-two receiver")
    native_upper, audit_lower = _receiver_bounds(two_step_paths)

    eigenvalues, modes = np.linalg.eigh(laplacian.toarray())
    audit_energies = np.sqrt(CAGE_REST_SQUARED + SPECTRAL_AUDIT_SPEED_SQUARED * eigenvalues)
    native = _cell_probabilities(eigenvalues, modes, sender)
    audit = _cell_probabilities(audit_energies - CAGE_REST, modes, sender)
    require(native[receiver] <= float(native_upper) + TOLERANCE
            and audit[receiver] >= float(audit_lower) - TOLERANCE
            and audit[receiver] > native[receiver],
            "the receiver read disagrees with the certified arrival separation")
    parked = np.zeros(federation.cells)
    parked[sender] = 1
    return {
        "level": level,
        "cells": federation.cells,
        "sender_cell": sender,
        "receiver_cell": receiver,
        "two_step_paths": two_step_paths,
        "receiver_laplacian_element": laplacian_entry,
        "receiver_native_squared_shell_defect": two_step_paths,
        "native_cell_probabilities": native.tolist(),
        "slope_matched_spectral_audit_cell_probabilities": audit.tolist(),
        "parked_cell_probabilities": parked.tolist(),
        "receiver_native_probability_upper_bound": str(native_upper),
        "receiver_spectral_audit_probability_lower_bound": str(audit_lower),
        "receiver_spectral_audit_excess_lower_bound": str(audit_lower - native_upper),
    }


def run():
    """Read the finite obstruction using the shared cell-Laplacian dictionary."""
    require(abs(CAGE_REST ** 2 - CAGE_REST_SQUARED) < TOLERANCE,
            "the declared cage frequency changed")
    charged_slope = Fraction(EXTERNAL_WEIGHT) / PORTS_PER_CELL
    # The jth ordered charged level lies below MASS + charged_slope * lambda_j
    # by min-max on the uniform-cell subspace. The cubic Laplacian norm is six.
    charged_upper = 2 * Fraction(MASS) * charged_slope + 6 * charged_slope ** 2
    require(0 <= charged_upper and charged_upper ** 2 < 4 * CAGE_REST_SQUARED,
            "the certified squared-speed intervals are no longer disjoint")
    return {
        "momentum_dictionary": "ordered eigenvalues of the common cell dual Laplacian",
        "source_preparation": "a localized cage with all scalar links written",
        "cage_link_weight": 1,
        "cage_rest_gap": CAGE_REST,
        "charged_rest_gap": MASS,
        "charged_kinetic_compression_factor": str(charged_slope),
        "charged_squared_speed_upper_bound": str(charged_upper),
        "cage_squared_speed_lower_bound": 2 * CAGE_REST,
        "wait": str(WAIT),
        "spectral_audit_comparator": "positive square root with shared cage rest gap and low-momentum slope",
        "federations": [_federation_record(level) for level in range(3)],
    }
