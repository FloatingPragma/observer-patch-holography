"""Use native face incidence and compare two-current couplings on its three ports."""

import numpy as np
from numpy.testing import assert_allclose

from model.executable_model import cell, rotations

from ..phase_feedback import circuit


def analyzer_rail_census(adjacency):
    order_three = [
        p for p in rotations() if all(p[i] != i and p[p[p[i]]] == i for i in range(12))
    ]
    observed = []
    for seam in circuit.marked.CELL_SEAMS:
        a, b, plus, minus = circuit.marked.addresses(seam)
        for tile in ((a, b, plus, minus), (b, a, plus, minus)):
            for turn in order_three:
                copy = [turn[i] for i in tile]
                last = [turn[i] for i in copy]
                if len(set((*tile, *copy, *last))) == 12:
                    center = [tile[0], copy[0], last[0]]
                    observed.append(int(adjacency[np.ix_(center, center)].sum()) // 2)
    assert observed and not any(observed)
    return dict(
        order_three_turns=len(order_three),
        oriented_edge_tiles=60,
        disjoint_tilings=len(observed),
        largest_number_of_central_seams=max(observed),
    )


def reencode(address, turn, face, mass, modes, current, gap):
    """Move the same three characters to an actual triangular face."""
    face_basis = np.eye(12)[:, face]
    transport = face_basis @ address.T
    assert_allclose(turn @ transport, transport @ turn, atol=0)
    edge = cell().adjacency().toarray().astype(bool)
    p_face = face_basis @ face_basis.T
    direct = transport @ mass @ transport.T + gap * (np.eye(12) - p_face)
    support = np.triu(abs(direct) > 1e-12, 1)
    assert support.sum() == 3 and not np.any(support & ~edge)
    assert_allclose(
        transport.T @ direct @ transport,
        address @ (address.T @ mass @ address) @ address.T,
        atol=2e-14,
    )
    triangle = direct[np.ix_(face, face)]
    edge_indices = (np.array([1, 2, 0]), np.arange(3))
    loop = np.prod(triangle[edge_indices])
    gauge = np.diag(np.exp(1j * np.array([0.1, 0.3, -0.2])))
    rephased = gauge @ triangle @ gauge.conj().T
    assert_allclose(np.prod(rephased[edge_indices]), loop, atol=2e-16)
    zero_flux = np.linalg.eigvalsh(abs(triangle))
    assert_allclose(zero_flux[0], zero_flux[1], atol=2e-14)
    records = dict(
        ports=face,
        face_mass_seams=np.argwhere(support).tolist(),
        analyzer_census=analyzer_rail_census(edge.astype(int)),
        mass_hopping_loop_phase=float(np.angle(loop)),
        zero_loop_phase_energies=zero_flux.tolist(),
    )
    return direct, p_face, transport @ modes, current @ transport.T, records


def two_current_records(face, modes, neutral, current):
    """One-body family selection does not determine the two-current operator."""
    read = neutral.conj().T @ current[:, face]
    currents = np.einsum("ap,pi->pai", read, modes[list(face)])
    joined = currents.sum(axis=0)
    assert_allclose(joined, np.diag(joined.diagonal()), atol=2e-14)
    assert_allclose(abs(joined.diagonal()), 1, atol=2e-14)

    def contamination(kernel):
        # Indices are incoming charge, outgoing charge, neutral, antineutral.
        amplitude = np.einsum("pq,pbd,qai->idab", kernel, currents.conj(), currents)
        powers = abs(amplitude) ** 2
        intended = np.einsum("idid->id", powers)
        total = powers.sum(axis=(-2, -1))
        return 1 - intended / total

    coherent = contamination(np.ones((3, 3)))
    independent = contamination(np.eye(3))
    assert_allclose(coherent, 0, atol=2e-14)
    assert_allclose(independent, 2 / 3, atol=2e-14)

    laplacian = 3 * np.eye(3) - np.ones((3, 3))
    comparisons = []
    for stiffness in (1, 8):
        # Unit squared mass is a demonstration setting of the mediator.
        kernel = np.linalg.inv(np.eye(3) + stiffness * laplacian)
        leakage = contamination(kernel)
        relative = 1 / (1 + 3 * stiffness)
        prediction = 2 * relative**2 / (1 + 2 * relative**2)
        assert_allclose(leakage, prediction, atol=2e-14)
        comparisons.append(dict(stiffness=stiffness, wrong_neutral_fraction=prediction))
    return dict(
        coherent_wrong_fraction=float(abs(coherent).max()),
        independent_wrong_fraction=float(independent[0, 0]),
        massive_mediator_examples=comparisons,
    )
