"""Resolve the energies supporting the two-winding fixed-port quantum channel."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from sympy import Matrix, Poly
from . import measure


def wait_kernel(energies, width):
    delta = energies[:, None] - energies[None, :]
    return np.exp(-1j * (16 + width / 2) * delta) * np.sinc(width * delta / (2 * np.pi))


def acquire(width, rng, shots=2**24):
    words, ix, terms = measure.build(True)
    h = sum(terms) / 8
    levels, v = np.linalg.eigh(h)
    weights = wait_kernel(levels, width)
    initial = ix[0, 0]
    code = [ix[3, a + 8 * b] for a in (0, 7) for b in (0, 7)]
    # Each branch records the full native receiver after the fixed-port herald.
    analyzers = [
        measure.face(words, ix, terms, True, 1, axis, 1).conj().T[:, code]
        for axis in range(3)
    ]
    outputs = []
    counts = []
    for axis in range(3):
        for sign in (1, -1):
            state = measure.face(words, ix, terms, True, 0, axis, sign)[:, initial]
            spectral = v.T @ state
            rho = v @ (weights * np.outer(spectral, spectral.conj())) @ v.T
            received = rho[np.ix_(code, code)]
            bank = []
            for decoder in analyzers:
                final = np.diag(decoder @ received @ decoder.conj().T).real
                plus, minus = [
                    sum(
                        x
                        for x, (site, _) in zip(final, words, strict=True)
                        if site == j
                    )
                    for j in (3, 4)
                ]
                np.testing.assert_allclose(plus + minus, np.trace(received), atol=1e-12)
                prob = np.maximum([plus, minus, 1 - plus - minus], 0)
                bank.append(rng.multinomial(shots, prob / prob.sum()).tolist())
            counts.append(bank)
            accepted = sum(a + b for a, b, _ in bank) / (3 * shots)
            means = [(a - b) / shots for a, b, _ in bank]
            outputs.append(
                (
                    accepted * np.eye(2)
                    + sum(x * p for x, p in zip(means, measure.PAULI, strict=True))
                )
                / 2
            )
    xplus, xminus, yplus, yminus, zplus, zminus = outputs
    cross = (xplus - xminus + 1j * (yplus - yminus)) / 2
    choi = np.block([[zplus, cross], [cross.conj().T, zminus]]) / 2
    # All widths share a 99% error envelope. Each of 18 settings per width
    # supplies one [-1,1] parity mean and one [0,1] acceptance frequency.
    eps = np.sqrt(2 * np.log(2 * 36 * len(WIDTHS) / 0.01) / shots)
    local_frobenius = np.sqrt(((eps / 2) ** 2 + 3 * eps**2) / 2)
    radius = np.sqrt(10) * local_frobenius / 2
    return dict(
        wait_width=width,
        original_starts_per_setting=shots,
        counts=counts,
        acceptance=choi.trace().real,
        pt_minimum=np.linalg.eigvalsh(measure.pt(choi))[0],
        simultaneous_radius=radius,
    )


WIDTHS = (0, 1, 2, 8)


def analyze():
    _, ix, terms = measure.build(True)
    h8 = sum(terms)
    polynomial = Poly(Matrix(h8.tolist()).charpoly().as_expr())
    simple = polynomial.gcd(polynomial.diff()).degree() == 0
    assert simple
    h = h8 / 8
    levels, vectors = np.linalg.eigh(h)
    all_poles = [
        measure.choi(measure.kraus(np.outer(col, col), ix)) for col in vectors.T
    ]
    for pole in all_poles:
        assert min(np.linalg.eigvalsh(measure.pt(pole))) > -2e-14
    fully_resolved = sum(all_poles)
    positive = vectors[:, levels > 0]
    positive_projector = positive @ positive.T
    # This comparison selects the positive spectral subspace mathematically.
    positive_coherent = measure.choi(
        measure.kraus(expm(-16j * h) @ positive_projector, ix)
    )
    rng = np.random.default_rng(492091)
    return dict(
        seed=492091,
        levels=levels.tolist(),
        simple=simple,
        positive_energies=int((levels > 0).sum()),
        coherent_positive_pt_min=np.linalg.eigvalsh(measure.pt(positive_coherent))[0],
        perfect_resolution_acceptance=fully_resolved.trace().real,
        perfect_resolution_pt_min=np.linalg.eigvalsh(measure.pt(fully_resolved))[0],
        arrival_weight_by_pole=[j.trace().real for j in all_poles],
        finite_waits=[acquire(w, rng) for w in WIDTHS],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(analyze(), indent=2) + "\n"
    if parser.parse_args().write:
        Path(__file__).with_name("rest_read.json").write_text(data)
    else:
        print(data, end="")
