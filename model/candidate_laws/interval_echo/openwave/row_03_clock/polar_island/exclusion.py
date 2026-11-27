"""A posteriori eigenvalue enclosures for all orientation-space terminal reads."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from model.candidate_laws.interval_echo.openwave.row_22_confinement.mobile_string import (
    material,
)
from model.experiments.material_inertia.pair_spectrum.neutral_guide.polar_ring import (
    measure as ring,
)


def max_abs_rows(matrix):
    return float(np.max(np.asarray(abs(matrix).sum(axis=1))))


def certify(n, kappa):
    _, blocks = ring.spectrum(n, kappa)
    dim = 1 << n
    words = material.bits(n)
    walls = (words != np.roll(words, 1, axis=1)).sum(axis=1)
    diagonal = kappa * (2 * walls + words.sum(axis=1))
    h = material.generate(diagonal.astype(float), n, material.H, True)
    assert abs(material.H) == 0.5
    # F is a permuted direct sum of small orbit Fourier matrices. Its
    # near-orthogonality is checked on the stored numbers rather than assumed
    # from the accuracy of the trigonometric library.
    f = sparse.hstack([b[0] for b in blocks], format="csr")
    assert f.shape == (dim, dim) and max(np.diff(f.indptr)) <= n
    assert float(abs(f.data).max()) < 1.01
    unit = 2.0**-53
    fourier_error = (
        2 * max_abs_rows(f @ f.conj().T - sparse.eye(dim)) + 64 * n**3 * unit
    )
    vector_error = 0
    squared_residual = 0
    for t, e, v, _, _ in blocks:
        d = len(e)
        assert max(np.diff(t.tocsr().indptr)) <= 1
        assert float(abs(v).max()) < 1.01
        gram = v.conj().T @ v - np.eye(d)
        vector_error = max(vector_error, 2 * max_abs_rows(gram) + 64 * d**3 * unit)
        q = t @ v
        residual = h @ q - q * e
        squared_residual += np.linalg.norm(residual, "fro") ** 2
    assert fourier_error < 0.01 and vector_error < 0.01
    unsorted = np.concatenate([b[1] for b in blocks])
    ordered = np.sort(unsorted)
    size_bound = 3 * kappa * n + n / 2
    # Proof.md bounds the residual of the exact product of the two
    # represented binary64 matrices, with an explicit roundoff allowance.
    roundoff = 64 * (n + 4) * unit * (size_bound + max(abs(ordered))) * np.sqrt(dim)
    error = 2 * (2 * np.sqrt(squared_residual) + roundoff + 1e-100)
    ground, rest = blocks[0][1][:2]
    assert ordered[0] == ground and np.count_nonzero(ordered == rest) == 1
    nearest_rest_neighbor = min(abs(ordered[ordered != rest] - rest))
    assert min(ordered[1] - ground, nearest_rest_neighbor) > 4 * error
    rest_gap = rest - ground
    exc = ordered[1:]
    target = exc + rest_gap
    index = np.searchsorted(exc, target)
    # Either neighbor at insertion can be closest. Clamping includes both
    # endpoints when the target lies above the largest energy.
    min_separation = min(
        np.min(abs(exc[np.clip(index + shift, 0, len(exc) - 1)] - target))
        for shift in [-1, 0]
    )
    # Four sorted eigenvalue errors plus subtraction/search arithmetic.
    detuning_bound = float(min_separation - 4 * error - 1e-10)
    assert detuning_bound > 0
    return dict(
        n=n,
        kappa=kappa,
        dimension=dim,
        rest_order_index=int(np.searchsorted(ordered, rest)),
        ground_to_rest=float(rest_gap),
        nearest_level_to_rest=float(nearest_rest_neighbor),
        fourier_gram_bound=fourier_error,
        modes_gram_bound=vector_error,
        observed_frobenius_residual=float(np.sqrt(squared_residual)),
        frobenius_roundoff_allowance=float(roundoff),
        ordered_eigenvalue_error=float(error),
        min_excited_gap_detuning=float(min_separation),
        strict_detuning_lower=detuning_bound,
    )


def run():
    return [certify(n, kappa) for n in (10, 12) for kappa in (8, 4)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output
    output.write_text(json.dumps(run(), indent=1, allow_nan=False) + "\n")
