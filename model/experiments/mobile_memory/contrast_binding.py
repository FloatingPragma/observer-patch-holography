"""Join native contrast binding and whole-cage transfer in one interaction."""

from math import isfinite
from numbers import Real

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.executable_model import carrier as load_carrier
from model.candidate_laws.contrast_memory_spherical_v1 import law as contrast
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.experiments.cage import construction as cage
from model.candidate_laws.scalar_channel.carrier import require
from model.experiments.contact.law import collar_states
from .cage_transfer import PORTS_PER_CELL, TOLERANCE, transfer


def _native_controls(kappa):
    """Read the link weight from the full echoed contrast-square operator."""
    exchange, square = contrast.local_term_parts()
    blank = np.eye(contrast.DIM)[:, :contrast.N_SRC]
    expected = np.array([(a - b) ** 2 for a, b in contrast.SRC])
    blank_weights = np.rint(np.diag(blank.T @ square @ blank)).astype(int)
    require(
        np.array_equal(blank_weights, expected)
        and np.linalg.norm(square @ blank - blank * blank_weights) < TOLERANCE,
        "the contact lost the native blank-memory restriction",
    )
    arms = collar_states(1)
    collar_weights = np.rint(np.diag(arms.T @ square @ arms)).astype(int)
    require(
        np.array_equal(collar_weights, [0, 1])
        and np.linalg.norm(square @ arms - arms * collar_weights) < TOLERANCE
        and np.linalg.norm((exchange + kappa * square) @ arms - kappa * arms * collar_weights)
        < TOLERANCE,
        "the empty-load collar arms do not give stationary link weights",
    )
    return blank_weights, collar_weights


def _activation(sector, seam, blank_weights):
    a, b = seam
    a_present = (sector.aa == a) | (sector.bb == a)
    b_present = (sector.aa == b) | (sector.bb == b)
    local_index = 3 * a_present.astype(int) + b_present
    return blank_weights[local_index]


def _covered_cells(carrier, level):
    scalar_carrier = cage.federation(level)
    require(scalar_carrier["edges"] == carrier.seams, "load and scalar carriers differ")
    cells = set()
    for ring in cage.exact_rings(scalar_carrier):
        cell = ring[0] // PORTS_PER_CELL
        require(
            all(port // PORTS_PER_CELL == cell for port in ring),
            "a degree-five ring crosses an external seam",
        )
        cells.add(cell)
    require(cells == set(range(carrier.cells)), "some cells have no exact cage")
    return len(cells)


def _evolve(hamiltonian, start, wait):
    action = -1j * wait
    return expm_multiply(action * hamiltonian, start, traceA=action * hamiltonian.diagonal().sum())


def _joint_record(sector, symmetric, antisymmetric):
    """Ideal terminal presence and cage records with one retained external term.

    The two eigenvalues of the single cage link are 0 and 2. Each block
    propagates all binary load configurations under the existing native law.
    """
    sender = abs((symmetric + antisymmetric) / 2) ** 2
    receiver = abs((symmetric - antisymmetric) / 2) ** 2
    near = sector.adjacent
    record = np.array([
        [sender[~near].sum(), receiver[~near].sum()],
        [sender[near].sum(), receiver[near].sum()],
    ])
    require(abs(record.sum() - 1.0) < 1e-9, "the joint projective records lost outcomes")
    return record


def _external_laplacian(configurations, cells, activations):
    """Sparse weighted cell Laplacian on each binary-load configuration."""
    rows, columns, values = [], [], []
    for (a, b), weights in activations:
        occupied = np.flatnonzero(weights)
        left = occupied * cells + a // PORTS_PER_CELL
        right = occupied * cells + b // PORTS_PER_CELL
        weight = weights[occupied]
        rows.extend((left, left, right, right))
        columns.extend((left, right, left, right))
        values.extend((weight, -weight, -weight, weight))
    dimension = configurations * cells
    contact = sparse.coo_matrix(
        (np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
        shape=(dimension, dimension),
    ).tocsr()
    require(
        not np.any(np.asarray(contact.sum(axis=1)))
        and float(abs(contact).sum(axis=1).max()) <= 4,
        "the full link action is outside the weighted Laplacian bound",
    )
    return contact


def _all_links_probabilities(sector, source, contact, start, sender, coupling, wait):
    cells = sector.car.cells
    ket = np.zeros((sector.P, cells))
    ket[:, sender] = start
    hamiltonian = sparse.kron(source, sparse.eye(cells), format="csr") + coupling * contact
    amplitude = _evolve(hamiltonian, ket.ravel(), wait)
    probability = abs(amplitude.reshape(ket.shape)) ** 2
    require(abs(probability.sum() - 1.0) < 1e-9, "the all-cell projective records lost outcomes")
    return probability


def _all_links_record(sector, source, contact, start, sender, coupling, wait):
    probability = _all_links_probabilities(sector, source, contact, start, sender, coupling, wait)
    near = sector.adjacent
    return np.stack([probability[~near].sum(axis=0), probability[near].sum(axis=0)])


def _addressed_pair(sector, seam, blank_weights):
    """Prepare one adjacent load pair that opens the addressed external term."""
    carrier = sector.car
    a = seam[0]
    sender = a // PORTS_PER_CELL
    neighbor = next(
        port for port in carrier.neighbours[a] if port // PORTS_PER_CELL == sender
    )
    pair_state = sector.index(a, neighbor)
    activation = _activation(sector, seam, blank_weights)
    require(
        sector.adjacent[pair_state] and activation[pair_state] == 1,
        "the supplied adjacent preparation does not activate the selected link",
    )
    start = np.zeros(sector.P)
    start[pair_state] = 1.0
    return sorted((a, neighbor)), start, activation


def run(kappa=128.0, coupling=1.0, wait=0.5, *, seam=None):
    """Certify coexistence and replay its 20-cell joint instruments.

    The common source adds g times the sum of the external contrasts A_s
    multiplied by their positive cage-exchange terms. The ordered-block bound
    applies when all terms are active or a specified subset is enabled.
    The instruments share an ideal adjacent load pair and an ideal sharp cage.
    They compare one enabled term with all terms, and read pair presence
    jointly with cage location in every available terminal cell.
    """
    require(
        isinstance(kappa, Real) and isfinite(kappa) and kappa >= 0,
        "kappa must be nonnegative and finite",
    )
    require(
        isinstance(coupling, Real) and isfinite(coupling) and coupling >= 0,
        "coupling must be nonnegative and finite",
    )
    require(isinstance(wait, Real) and isfinite(wait) and wait >= 0, "invalid wait")
    gap = 2 * kappa - 6 - 4 * coupling
    require(isfinite(gap) and gap > 0, "the ordered-block separation bound needs a positive gap")
    bound = min(1.0, 40 / gap ** 2)

    level = 0
    carrier = load_carrier.federation(level)
    sector = PairSector(carrier)
    covered = _covered_cells(carrier, level)
    endpoints = [port for seam in carrier.external for port in seam]
    require(
        len(endpoints) == len(set(endpoints)),
        "external seams do not form the required port matching",
    )
    blank_weights, collar_weights = _native_controls(kappa)
    activations = [(link, _activation(sector, link, blank_weights)) for link in carrier.external]
    active_count = sum(weights for _, weights in activations)
    require(
        sector.P == len(active_count) and 0 <= active_count.min()
        and active_count.max() <= 2,
        "a binary pair activates more than two external cage links",
    )

    seam = carrier.external[0] if seam is None else seam
    require(
        isinstance(seam, (tuple, list)) and len(seam) == 2
        and all(type(port) is int for port in seam)
        and tuple(sorted(seam)) in carrier.external,
        "the selected seam must be an oriented external seam on the 20-cell carrier",
    )
    sender = seam[0] // PORTS_PER_CELL
    initial_pair, start, activation = _addressed_pair(sector, seam, blank_weights)
    source = sector.hamiltonian(kappa)
    symmetric = _evolve(source, start, wait)
    drive = source + sparse.diags(2 * coupling * activation, format="csr")
    antisymmetric = _evolve(drive, start, wait)
    joint = _joint_record(sector, symmetric, antisymmetric)
    additive = _joint_record(sector, symmetric, symmetric)
    require(
        float(joint[0].sum()) <= bound + TOLERANCE
        and float(additive[:, 1].sum()) < TOLERANCE,
        "joint records violate the pair bound or additive cage control",
    )
    contact = _external_laplacian(sector.P, carrier.cells, activations)
    all_additive = np.zeros((2, carrier.cells))
    all_additive[:, sender] = additive[:, 0]
    all_joint = (
        _all_links_record(sector, source, contact, start, sender, coupling, wait)
        if coupling else all_additive
    )
    require(
        float(all_joint[0].sum()) <= bound + TOLERANCE,
        "all-link records violate the pair bound",
    )

    empty_load_control = {}
    relay_wait = np.pi / (2 * coupling) if coupling else np.pi / 2
    for name, eigenvalue in zip(("blank", "written"), collar_weights):
        relay = transfer(
            level=level, seam=seam, memory=bool(eigenvalue),
            coupling=coupling, wait=relay_wait,
        )
        empty_load_control[name] = {
            "native_contrast_eigenvalue": int(eigenvalue),
            "receiver_record": relay["receiver_record"],
            "reference_record": relay["reference_record"],
            "exciton_record": relay["exciton_record"],
        }
    return {
        "cells_with_exact_cages": covered,
        "checked_binary_configurations": sector.P,
        "active_link_count_histogram": np.bincount(active_count, minlength=3).tolist(),
        "contact_norm_bound": float(4 * coupling),
        "ordered_block_gap_lower_bound": float(gap),
        "all_time_no_pair_bound": float(bound),
        "seam": list(seam),
        "initial_load_pair": initial_pair,
        "wait": float(wait),
        "rows": ["no_adjacent_pair", "adjacent_pair"],
        "columns": ["sender_cage", "receiver_cage"],
        "joint_record": joint.tolist(),
        "additive_record": additive.tolist(),
        "all_link_cells": list(range(carrier.cells)),
        "all_link_basis_dimension": contact.shape[0],
        "all_link_joint_record": all_joint.tolist(),
        "all_link_additive_record": all_additive.tolist(),
        "empty_load_control": empty_load_control,
    }
