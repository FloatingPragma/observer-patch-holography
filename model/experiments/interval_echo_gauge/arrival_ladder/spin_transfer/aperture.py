"""A common-donor Lüders measurement of a real rank-two vacancy subspace."""

from math import atan2, pi, sin, cos

import numpy as np

from model.executable_model import federation
from model.experiments.interval_echo_lepton_register import measure as rails
from model.experiments.interval_echo_gauge.arrival_ladder import measure as ladder


def _real_basis(projector):
    p = np.asarray(projector)
    if (
        p.shape != (12, 12)
        or not np.isfinite(p).all()
        or np.max(np.abs(p.imag)) > 1e-12
    ):
        raise ValueError("Expected a finite real 12-port projector")
    p = p.real
    if (
        not np.allclose(p, p.T, rtol=0, atol=1e-12)
        or not np.allclose(p @ p, p, rtol=0, atol=1e-12)
        or abs(p.trace() - 2) > 1e-12
    ):
        raise ValueError("Expected an orthogonal projector of rank two")
    # Choose two well-conditioned actual projected port reads. Eigensolvers
    # may choose unrelated bases of the exactly degenerate accepted space.
    pivot = int(np.argmax(np.diag(p)))
    u = p[:, pivot] / np.linalg.norm(p[:, pivot])
    other = p - np.outer(u, u)
    pivot = int(np.argmax(np.diag(other)))
    v = other[:, pivot] / np.linalg.norm(other[:, pivot])
    return np.column_stack((u, v))


def _compile_basis(projector, owner):
    coordinates = _real_basis(projector)
    operations = []
    # The order of the two surviving rails is arbitrary. The common meter
    # sees their union, so no eigenvector label is ever needed on acceptance.
    for column in range(2):
        for row in range(column + 1, len(coordinates)):
            a, b = coordinates[column, column], coordinates[row, column]
            if abs(b) < 1e-14:
                continue
            angle = atan2(-b, a)
            rotation = np.array([[cos(angle), -sin(angle)], [sin(angle), cos(angle)]])
            coordinates[[column, row]] = rotation @ coordinates[[column, row]]
            operations.append((column, row, angle))
    forward = np.eye(12, dtype=complex)
    backward = np.eye(12, dtype=complex)
    for i, j, angle in operations:
        rails.rotate(forward, i, j, angle, 1, owner)
    for i, j, angle in reversed(operations):
        rails.rotate(backward, i, j, -angle, 1, owner)
    accepted = np.diag([1, 1] + [0] * 10)
    np.testing.assert_allclose(
        forward @ projector @ forward.conj().T, accepted, atol=3e-12, rtol=0
    )
    np.testing.assert_allclose(backward @ forward, np.eye(12), atol=3e-12, rtol=0)
    return forward, backward, operations


def compact_instrument(projector, owner, *, resolve_address=False):
    """Measure the accepted two-rail subspace, then return the common donor.

    The accepted and refused matrices act in the original twelve-port basis.
    They are the actual recovered-donor Kraus amplitudes. If resolve_address
    is set, the acceptance operation has two Kraus matrices, supplied in
    accepted_kraus, and accepted is None. The refusal has one Kraus matrix.
    Numerical residuals on the other final donor rail are bounded separately.
    """
    board = federation(0)
    if owner < 0 or owner >= board.cells:
        raise ValueError("Owner outside the twenty-cell federation")
    ports = range(12 * owner, 12 * owner + 12)
    if set(ports) & set(ladder.layout()["offline"]):
        raise ValueError("The synthesized cell overlaps the donor workstation")
    rotate, restore, commands = _compile_basis(projector, owner)
    flag_ports = list(ports[:2])
    hit = ladder.aperture_map(flag_ports)
    if ladder.read_loads() != [0, 1, 2]:
        raise ArithmeticError("The donor's complete load pointer failed calibration")

    # When the vacancy is outside both rails, both visits see an ordinary
    # unit. The successive local words have the common splitter and its
    # inverse adjoining each other, so their product realizes the same sweep.
    ix = [ladder.WORDS.index(w) for w in ((1, 1, 2), (1, 2, 1))]
    filled = ladder.interaction(pi, ladder.K)[np.ix_(ix, ix)]
    miss = filled @ filled
    blocks = {}
    for key, (i, j) in {"cc": (0, 0), "cn": (0, 1), "nc": (1, 0), "nn": (1, 1)}.items():
        names = ("charged", "neutral")
        label = names[j] + "_to_" + names[i]
        diagonal = np.full(12, miss[i, j], complex)
        diagonal[:2] = np.diag(hit[label])
        blocks[key] = np.diag(diagonal)

    # The intermediate pointer reads only the common donor's sharp load.
    # Applying the same physical two-visit word resets it on both outcomes.
    accept_in_rails = blocks["cn"] @ blocks["nc"]
    refuse_in_rails = blocks["cc"] @ blocks["cc"]
    accepted = restore @ accept_in_rails @ rotate
    refused = restore @ refuse_in_rails @ rotate
    residues = [
        restore @ blocks["nn"] @ blocks["nc"] @ rotate,
        restore @ blocks["nc"] @ blocks["cc"] @ rotate,
    ]
    residual = max(np.linalg.norm(z, 2) for z in residues)
    completeness = sum(z.conj().T @ z for z in (accepted, refused, *residues))
    np.testing.assert_allclose(completeness, np.eye(12), rtol=0, atol=8e-12)
    np.testing.assert_allclose(accepted, projector, rtol=0, atol=5e-12)
    np.testing.assert_allclose(refused, np.eye(12) - projector, rtol=0, atol=5e-12)
    branches = [accepted]
    if resolve_address:
        # Deliberately write the address within the accepted pair and forget it.
        # The physical refusal still carries no record of any of the ten rails.
        branches = []
        for r in range(2):
            single = np.zeros((12, 12), complex)
            single[r, r] = accept_in_rails[r, r]
            branches.append(restore @ single @ rotate)
    return dict(
        accepted=None if resolve_address else accepted,
        refused=refused,
        accepted_kraus=branches,
        refused_kraus=[refused],
        max_donor_residue=float(residual),
        flag_ports=flag_ports,
        rotation_count=len(commands),
        routes=hit["routes"],
    )
