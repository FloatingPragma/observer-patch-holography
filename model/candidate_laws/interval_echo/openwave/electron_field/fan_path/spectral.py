"""Calibrate separated ring channels and certify the released level count."""

from fractions import Fraction
from itertools import product
from math import floor

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from ..binding_audit.measure import count_below


def charged_factor(path, cell):
    field = path.field
    basis = [(x, w) for x, words in enumerate(path.words_at) for w in words]
    ids = dict(zip(basis, range(len(basis)), strict=True))
    h = np.zeros((len(basis), len(basis)))
    for i, (x, w) in enumerate(basis):
        h[i, i] = 8 * (6 - int(path.board.degree[path.cells[cell][x]]))
        for mask, cw in zip(field.masks, field.clockwise, strict=True):
            if w & mask in (cw, cw ^ mask):
                h[i, i] += 0.5
                h[i, ids[x, w ^ mask]] -= 0.5
        for (a, b), bit in field.bits.items():
            if x in (a, b) and ((x == a) == bool(w & bit)) == (path.sign == 1):
                other = b if a == x else a
                h[i, i] += 4
                h[i, ids[other, w ^ bit]] -= 4
    np.testing.assert_array_equal(h, h.T)
    return basis, h


def separated_reference(path):
    """Check each factorization directly, including the altered transit words."""
    sectors = []
    thresholds = []
    brackets = []
    for c in range(path.size):
        charged_basis, charged = charged_factor(path, c)
        minimum = float(np.linalg.eigvalsh(charged)[0])
        lo = Fraction(floor(minimum * 1024), 1024)
        hi = lo + Fraction(1, 1024)
        assert count_below(charged, lo) == 0 and count_below(charged, hi) >= 1
        neutral_cells = [i for i in range(path.size) if i != c]
        neutral_words = [path.words_at[path.references(c, 0)[i]] for i in neutral_cells]
        separate = sparse.csr_matrix(charged)
        gaps = []
        for words in neutral_words:
            ring = path.field.ring_laplacian(words)
            # The checked spectator sectors have two or four sharp words.
            # Integer polynomial and exact trace certify their spectra.
            allowed = (0, 1) if len(words) == 2 else (0, 1, 2)
            polynomial = np.eye(len(words))
            for e in allowed:
                polynomial = polynomial @ (ring - e * np.eye(len(words)))
            np.testing.assert_array_equal(polynomial, 0)
            assert (len(words), np.trace(ring), np.trace(ring @ ring)) in (
                (2, 1, 1),
                (4, 4, 6),
            )
            gaps.append(1)
            separate = sparse.kron(
                separate, sparse.eye(len(words)), format="csr"
            ) + sparse.kron(
                sparse.eye(separate.shape[0]), sparse.csr_matrix(ring), format="csr"
            )
        order = []
        for x, w in charged_basis:
            for nv in product(*neutral_words):
                combined = nv[:c] + (w,) + nv[c:]
                order.append(path.ids[path.cells[c][x], combined])
        assert len(set(order)) == sum(path.owner == c)
        assert (path.parked[order][:, order] != separate).nnz == 0
        assert path.parked[order][:, path.owner != c].nnz == 0
        gap = min(gaps)
        sectors.append(
            {
                "occupied_hub": path.cells[c][0],
                "charged_ground": minimum,
                "charged_ground_dyadic_bracket": [str(lo), str(hi)],
                "remote_hubs": [path.cells[i][0] for i in neutral_cells],
                "remote_word_counts": list(map(len, neutral_words)),
                "lowest_remote_positive_ring_energy": gap,
            }
        )
        thresholds.append(minimum + gap)
        brackets.append((lo + gap, hi + gap))
    return (
        min(thresholds),
        (min(a for a, _ in brackets), min(b for _, b in brackets)),
        sectors,
    )


def inf_norm(matrix):
    return float(np.max(abs(matrix).sum(axis=1)))


def verified_inertia(h, cut):
    """Certify inertia via a symmetric congruence and a verified approximate inverse."""
    n = h.shape[0]
    assert n <= 8192 and cut.denominator & (cut.denominator - 1) == 0
    assert np.all(np.isfinite(h.data)) and np.max(abs(h.data)) < 2**20
    np.testing.assert_array_equal(h.data * 2, np.rint(h.data * 2))
    scale = max(2, cut.denominator)
    assert scale * (np.max(abs(h.data)) + abs(float(cut))) < 2**52
    matrix = h - float(cut) * sparse.eye(n, format="csr")
    scaled = h.multiply(scale).astype(np.int64) - int(cut * scale) * sparse.eye(
        n, format="csr", dtype=np.int64
    )
    assert (matrix.multiply(scale) != scaled).nnz == 0
    assert (matrix != matrix.T).nnz == 0
    lu = splu(
        matrix.tocsc(),
        permc_spec="MMD_AT_PLUS_A",
        diag_pivot_thresh=0,
        options={"SymmetricMode": True, "Equil": False},
    )
    np.testing.assert_array_equal(lu.perm_r, lu.perm_c)
    order = np.argsort(lu.perm_r)
    lower = lu.L.tocsr()
    np.testing.assert_array_equal(lower.diagonal(), 1)
    assert sparse.triu(lower, 1).nnz == 0
    d = lu.U.diagonal()
    assert np.min(abs(d)) > 0 and np.max(abs(d)) < 2**30
    assert np.max(abs(lower.data)) < 2**30
    factor = (lower.multiply(d) @ lower.T).tocsr()
    measured_error = inf_norm(factor - matrix[order][:, order])
    abs_factor = np.asarray(
        abs(lower) @ (abs(d) * np.asarray(abs(lower).sum(axis=0)).ravel())
    )
    norm = inf_norm(matrix)
    # For n<=8192 this dominates eight ordinary binary64 dot or sum
    # error factors, including sparse products and sums accumulated in batches.
    guard = 2.0**-32
    pad = 2.0**-800
    factor_error = (measured_error + guard * (norm + np.max(abs_factor)) + pad) * (
        1 + guard
    )
    inverse_sums = np.zeros(n)
    residual_sums = np.zeros(n)
    for start in range(0, n, 64):
        cols = np.arange(start, min(n, start + 64))
        unit = np.zeros((n, len(cols)), order="F")
        unit[cols, np.arange(len(cols))] = 1
        inverse = lu.solve(unit)
        assert np.max(abs(inverse)) < 2**30
        inverse_sums += np.sum(abs(inverse), axis=1)
        residual_sums += np.sum(abs(unit - matrix @ inverse), axis=1)
    xnorm = (1 + guard) * (np.max(inverse_sums) + pad)
    residual = (1 + guard) * (np.max(residual_sums) + guard * (1 + norm * xnorm) + pad)
    assert residual < 1
    margin = (1 + guard) * xnorm * factor_error / (1 - residual) + pad
    assert margin < 1
    return {
        "dyadic_cut": str(cut),
        "count_below": int(np.count_nonzero(d < 0)),
        "factor_error_infinity_upper": float(factor_error),
        "inverse_residual_infinity_upper": float(residual),
        "inverse_times_factor_error_upper": float(margin),
    }
