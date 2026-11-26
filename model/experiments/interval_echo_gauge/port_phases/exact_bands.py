"""Certify both leading bands and the twelve compressed generators over Q(zeta_20)."""

import json
import numpy as np

from model.executable_model import cell
from model.experiments.interval_echo_angular.charged_orientation.flux_doublet.measure import (
    links,
)

# Basis 1,z,...,z^7; Phi_20(z)=z^8-z^6+z^4-z^2+1.
D = 8


def reduce(p):
    p = p.copy()
    for k in reversed(range(D, len(p))):
        for shift, sign in ((0, -1), (2, 1), (4, -1), (6, 1)):
            p[k - D + shift] += sign * p[k]
    return p[:D]


def product(a, b):
    # Scalars commute with matrices, while matrix order is preserved.
    shape = a.shape[1:] or b.shape[1:]
    p = np.zeros((2 * D - 1, *shape), dtype=object)
    for i in range(D):
        for j in range(D):
            p[i + j] += a[i] @ b[j] if a.ndim == b.ndim == 3 else a[i] * b[j]
    return reduce(p)


def root(k):
    p = np.zeros(20, dtype=object)
    p[k % 20] = 1
    return reduce(p)


def evaluate(a, z, prime):
    x = np.zeros_like(a[0])
    for c in reversed(a):
        x = (x * z + c) % prime
    return x


def gram_determinant(matrix, prime):
    m = matrix.copy() % prime
    det = 1
    for c in range(len(m)):
        pivot = next((j for j in range(c, len(m)) if m[j, c]), None)
        if pivot is None:
            return 0
        if pivot != c:
            m[[pivot, c]] = m[[c, pivot]]
            det = -det
        diag = int(m[c, c])
        det = det * diag % prime
        inverse = pow(diag, -1, prime)
        for j in range(c + 1, len(m)):
            m[j] = (m[j] - (m[j, c] * inverse % prime) * m[c]) % prime
    return det


def certify():
    one = root(0)
    phi = root(2) + root(-2)  # (1+sqrt(5))/2
    beta = root(3) + root(-3)  # 2 cos(3pi/10) > 0
    prime = 41
    z = next(
        a
        for a in range(2, prime)
        if pow(a, 20, prime) == 1 and pow(a, 10, prime) != 1 and pow(a, 4, prime) != 1
    )
    assert (
        pow(z, 8, prime) - pow(z, 6, prime) + pow(z, 4, prime) - pow(z, 2, prime) + 1
    ) % prime == 0
    gram = np.zeros((12, 12), dtype=object)
    outcomes = []
    identity = np.eye(12, dtype=object)
    for k, rank, load, levels in [
        (2, 3, 2, (2 * one + phi, phi - 2 * one, one - 2 * phi)),
        (7, 2, -3, (product(2 * phi - one, beta), beta, -product(phi + one, beta))),
    ]:
        A = np.zeros((D, 12, 12), dtype=object)
        field = links(k)
        for a, b in cell().seams:
            A[:, b, a] = root(field[a, b])
            A[:, a, b] = root(-field[a, b])
        alpha, beta_k, gamma = levels

        def shifted(x):
            return A - x[:, None, None] * identity

        # Three distinct real eigenvalues; alpha is largest for the physical root.
        Q = product(shifted(beta_k), shifted(gamma))
        denominator = product(alpha - beta_k, alpha - gamma)
        assert not np.any(product(shifted(alpha), Q))
        assert np.array_equal(product(Q, Q), product(denominator, Q))
        for p in range(12):
            assert np.array_equal(12 * Q[:, p, p], rank * denominator)
        den = int(evaluate(denominator, z, prime))
        assert den
        P = evaluate(Q, z, prime) * pow(den, -1, prime) % prime
        assert np.array_equal(P @ P % prime, P)
        gram = (gram + load**2 * P * P.T) % prime
        outcomes.append(
            dict(flux=k, leading_rank=rank, charge=load, denominator_mod_prime=den)
        )
    det = gram_determinant(gram, prime)
    assert det
    return dict(
        exact_leading_bands=outcomes,
        prime=prime,
        primitive_root_image=z,
        compressed_current_gram_determinant_mod_prime=det,
        tracefree_current_dimension=12,
    )


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2))
