"""Exact registered formation generators for the geometric release and a written face."""

from fractions import Fraction

import numpy as np
from scipy.sparse import block_diag, csr_matrix, diags, eye, kron

from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.face_frame import (
    measure as face_frame,
)
from model.executable_model import cell

from ..phase_feedback import circuit

KAPPA = int(circuit.marked.KAPPA)


def geometry():
    tiles, rotation, _, _, _, fourier, _, _ = circuit.build()
    face = next(
        f for f in face_frame.faces() if list(face_frame.face_turn(f)) == rotation
    )
    tags = np.eye(12)[:, tiles[:, 0]] @ fourier
    bands, _ = face_frame.theorem_projectors()
    frames = face_frame.character_frames(face, bands[1][0])
    # These are the same three rotation eigenvalues as the original Fourier tags.
    chosen = frames[[0, 2, 1], :, face[0]].T
    band_tags = chosen / np.linalg.norm(chosen, axis=0)
    turn = face_frame.permutation(rotation)
    orientation = np.diag(tags.conj().T @ turn @ tags)
    np.testing.assert_allclose(turn @ band_tags, band_tags * orientation, atol=2e-14)
    return face, tags, band_tags


def blank_generator(board, sign, private):
    """Eight times the one-defect operator relative to the matched filled host."""
    reference = (6 - board.degree) if private else np.zeros(board.n, dtype=int)
    diagonal = 8 * KAPPA * (board.degree + reference * ((1 + sign) ** 2 - 1))
    laplacian = board.laplacian().astype(np.int64)
    return (4 * laplacian + diags(diagonal, dtype=diagonal.dtype)).tocsr()


def written_generator(board, face, sign, private):
    """Attach the checked written face to an otherwise all-open blank federation."""
    local, _ = face_frame.native(face, sign, kappa=KAPPA)
    blank_cell = 4 * cell().laplacian().toarray() + 40 * KAPPA * np.eye(12, dtype=int)
    change = local - np.kron(blank_cell, np.eye(8, dtype=int))
    np.testing.assert_array_equal(change, change.astype(np.int64))
    outside = 8 * (board.n - 12)
    attachment = block_diag(
        (
            csr_matrix(change, dtype=np.int64),
            csr_matrix((outside, outside), dtype=np.int64),
        )
    )
    blank = kron(blank_generator(board, sign, private), eye(8, dtype=np.int64))
    result = (blank + attachment).tocsr()
    result.eliminate_zeros()
    assert (result - result.T).nnz == 0
    return result


def energy_certificate(matrix, soft_ports=(), head_count=1):
    """Bound every possible formation frequency, including frequencies lost to filtering."""
    diagonal = matrix.diagonal()
    off_diagonal = matrix - diags(diagonal, dtype=diagonal.dtype)
    assert off_diagonal.data.max(initial=0) == 0
    row_radius = np.asarray(abs(off_diagonal).sum(axis=1)).ravel()
    upper = Fraction(int(max(diagonal + row_radius)), 8)
    row_lower = Fraction(int(min(diagonal - row_radius)), 8)
    weight = np.ones(matrix.shape[0], dtype=np.int64)
    weight.reshape(-1, head_count)[list(soft_ports)] = 10
    applied = matrix @ weight
    scaled_lower = min(
        Fraction(int(min(applied[weight == value])), 8 * int(value))
        for value in np.unique(weight)
    )
    lower = max(row_lower, scaled_lower)
    assert 0 < lower <= upper
    ratio = upper / lower
    # (4 + 3 sqrt(2))^2 is the smallest possible span for exact positive Koide balance.
    assert ratio <= 64
    r = float(ratio)
    return dict(
        row_lower_exact=str(row_lower),
        weighted_lower_exact=str(scaled_lower),
        formation_support_exact=[str(lower), str(upper)],
        formation_support=[float(lower), float(upper)],
        heaviest_to_lightest_ceiling=float(ratio),
        koide_ceiling=(r + 2) / (np.sqrt(r) + 2) ** 2,
    )


def preparations(board, tags, head_angle=None):
    if head_angle is None:
        result = np.zeros((board.n, 3), complex)
        result[:12] = tags
    else:
        head = np.exp(
            1j * head_angle * np.array([value.bit_count() for value in range(8)])
        )
        result = np.zeros((8 * board.n, 3), complex)
        result[:96] = np.kron(tags, head[:, None] / np.sqrt(8))
    np.testing.assert_allclose(result.conj().T @ result, np.eye(3), atol=2e-14)
    return result
