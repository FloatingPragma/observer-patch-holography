"""Exact one-cage refinement certificate for the stored geodesic tower."""

from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply

from model.executable_model import carrier
from model.candidate_laws.scalar_channel.carrier import require
from .cage_transfer import PORTS_PER_CELL
from .contrast_binding import _covered_cells


CHILDREN = 4
CENTER = 3
WAIT = Fraction(1, 32)
TOLERANCE = 1e-12


def _cell_seams(federation):
    return [
        tuple(sorted((a // PORTS_PER_CELL, b // PORTS_PER_CELL)))
        for a, b in federation.external
    ]


def _neighbors(size, seams):
    edges = set(seams)
    require(len(edges) == len(seams), "the cell dual has a repeated seam")
    neighbors = [set() for _ in range(size)]
    for a, b in edges:
        require(0 <= a < b < size, "the cell dual has an invalid seam")
        neighbors[a].add(b)
        neighbors[b].add(a)
    require(all(len(adjacent) == 3 for adjacent in neighbors), "the cell dual is not cubic")
    return neighbors


def _laplacian(neighbors):
    count = len(neighbors)
    rows = [i for i, adjacent in enumerate(neighbors) for _ in adjacent]
    columns = [j for adjacent in neighbors for j in sorted(adjacent)]
    adjacency = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.int64), (rows, columns)), shape=(count, count)
    )
    return 3 * sparse.eye(count, dtype=np.int64, format="csr") - adjacency


def _zero(matrix):
    return not np.any(matrix.data)


def certify_step(coarse_cells, coarse_seams, fine_cells, fine_seams):
    """Check the parent map and both exact integer matrix identities."""
    require(coarse_cells > 0 and fine_cells == CHILDREN * coarse_cells,
            "one refinement must have four children for each parent")
    coarse_neighbors = _neighbors(coarse_cells, coarse_seams)
    fine_neighbors = _neighbors(fine_cells, fine_seams)
    crossings = defaultdict(list)
    for a, b in fine_seams:
        left, right = a // CHILDREN, b // CHILDREN
        if left != right:
            crossings[left, right].append((a, b))
    require(
        {edge: len(children) for edge, children in crossings.items()}
        == {edge: 2 for edge in coarse_seams},
        "a coarse seam must have exactly two actual fine crossings",
    )

    for parent, neighbors in enumerate(coarse_neighbors):
        first = CHILDREN * parent
        corners = set(range(first, first + CENTER))
        middle = first + CENTER
        require(fine_neighbors[middle] == corners, "the central child lost its internal star")
        missing = []
        for child in corners:
            require(middle in fine_neighbors[child], "the corner lacks its central seam")
            other_parents = {j // CHILDREN for j in fine_neighbors[child] - {middle}}
            require(len(other_parents) == 2 and other_parents < neighbors,
                    "each corner must cross to two different adjacent parents")
            missing.extend(neighbors - other_parents)
        require(Counter(missing) == Counter(neighbors),
                "the three coarse seams must omit three different crossing corners")

    l_coarse = _laplacian(coarse_neighbors)
    l_fine = _laplacian(fine_neighbors)
    require(
        connected_components(l_coarse, directed=False, return_labels=False) == 1
        and connected_components(l_fine, directed=False, return_labels=False) == 1,
        "the cell dual must be connected at both levels",
    )
    assignment = sparse.csr_matrix(
        (np.ones(fine_cells, dtype=np.int64),
         (np.arange(fine_cells), np.arange(fine_cells) // CHILDREN)),
        shape=(fine_cells, coarse_cells),
    )
    # P = assignment/2 is the isometric lift. Numerator equals 4 times
    # (1 - P P^T) L_fine P, with integer coefficients.
    require(_zero(assignment.T @ l_fine @ assignment - 2 * l_coarse),
            "the cell lift failed exact half-Laplacian compression")
    leakage = 2 * l_fine @ assignment - assignment @ l_coarse
    require(_zero(assignment.T @ leakage)
            and _zero(leakage.T @ leakage - 8 * l_coarse),
            "the cell lift failed the exact leakage Gram identity")
    return l_coarse, l_fine, crossings


def _record(laplacian, state):
    step = -1j * float(WAIT) * laplacian
    result = expm_multiply(step, state, traceA=step.diagonal().sum())
    probabilities = abs(result) ** 2
    require(abs(probabilities.sum() - 1) < TOLERANCE, "the cage record lost probability")
    return probabilities


def _detector_gap_bound():
    """Rational lower bound from imaginary and real Taylor remainders."""
    t = WAIT
    # Two fine crossing children each receive imaginary amplitude t/2
    # at first order. Use 7/10 < 1/sqrt(2) for its receiver-vector norm.
    fine_lower = max(Fraction(0), Fraction(7, 10) * t - (6 * t) ** 3 / 6) ** 2
    coarse_upper = (t / 2 + (3 * t) ** 3 / 6) ** 2 + ((3 * t) ** 2 / 2) ** 2
    return fine_lower - coarse_upper


def _step_record(level, coarse, fine):
    coarse_seams = _cell_seams(coarse)
    fine_seams = _cell_seams(fine)
    l_coarse, l_fine, crossings = certify_step(coarse.cells, coarse_seams, fine.cells, fine_seams)

    sender = 0
    receiver = min(b for a, b in coarse_seams if a == sender)
    coarse_prep = np.zeros(coarse.cells)
    coarse_prep[sender] = 1
    fine_prep = np.repeat(coarse_prep / 2, CHILDREN)
    # In the one-cage sector all link records written give this Laplacian,
    # up to a scalar offset. The control sets the link coefficient to zero.
    fine_counts = _record(l_fine, fine_prep)
    coarse_counts = _record(l_coarse / 2, coarse_prep)
    parked_counts = abs(fine_prep) ** 2
    fine_parents = fine_counts.reshape(coarse.cells, CHILDREN).sum(axis=1)
    bound = _detector_gap_bound()
    difference = float(fine_parents[receiver] - coarse_counts[receiver])
    require(bound > 0 and difference >= float(bound) - TOLERANCE,
            "the neighbor detector failed the exact separation bound")
    return {
        "coarse_level": level,
        "coarse_cells": coarse.cells,
        "fine_cells": fine.cells,
        "coarse_seams": len(coarse_seams),
        "fine_crossings": 2 * len(coarse_seams),
        "children_per_parent": CHILDREN,
        "parent_of_fine_cell": [cell // CHILDREN for cell in range(fine.cells)],
        "seam_fibers": [
            {"coarse": list(edge), "fine": [list(child) for child in sorted(children)]}
            for edge, children in sorted(crossings.items())
        ],
        "within_parent_edges": CENTER * coarse.cells,
        "compression_factor": "1/2",
        "leakage_gram_factor": "1/2",
        "no_proper_parent_local_invariant_subspace": True,
        "wait": str(WAIT),
        "sender_parent": sender,
        "receiver_parent": receiver,
        "fine_cell_records": fine_counts.tolist(),
        "fine_parent_records": fine_parents.tolist(),
        "lifted_coarse_cell_records": np.repeat(coarse_counts / CHILDREN, CHILDREN).tolist(),
        "energy_matched_coarse_records": coarse_counts.tolist(),
        "zero_link_cell_records": parked_counts.tolist(),
        "zero_link_parent_records": parked_counts.reshape(coarse.cells, CHILDREN).sum(axis=1).tolist(),
        "receiver_difference": difference,
        "receiver_difference_lower_bound": str(bound),
    }


def run():
    """Return exact graph certificates and complete finite occupation records."""
    federations = [carrier.federation(level) for level in range(3)]
    coverage = [_covered_cells(federation, level) for level, federation in enumerate(federations)]
    return {
        "cells_with_exact_cages_by_level": coverage,
        "steps": [
            _step_record(level, federations[level], federations[level + 1])
            for level in (0, 1)
        ],
    }
