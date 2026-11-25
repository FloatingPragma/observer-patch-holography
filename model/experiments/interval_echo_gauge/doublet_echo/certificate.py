"""Exact detection of nonzero reciprocal echoes and a relative free-flight phase."""

from math import factorial
import numpy as np
from model.experiments.interval_echo_gauge.doublet_flight.certificate import (
    PRIME,
    modes,
)
from .measure import apparatus


def certificate(kappa, order=60):
    _, sector, addresses, ports, home, _ = apparatus()
    # Integer evaluation of the algebraic initial spectral projector at its
    # finite-field root, in a real unnormalized home-plane basis.
    vectors = np.zeros((sector.P, 2), dtype=np.int64)
    vectors[home] = modes(addresses[0])
    generator = (2 * sector.hamiltonian(kappa)).astype(np.int64)
    coefficients = []
    for n in range(order + 1):
        inv = pow(factorial(n) % PRIME, -1, PRIME)
        coefficients.append((vectors[ports].reshape(2, 2, 2) * inv) % PRIME)
        vectors = (generator @ vectors) % PRIME
    determinants = []
    for n in range(order + 1):
        value = np.zeros(2, dtype=np.int64)
        for a, b in zip(
            coefficients[: n + 1], reversed(coefficients[: n + 1]), strict=True
        ):
            value += a[:, 0, 0] * b[:, 1, 1] - a[:, 0, 1] * b[:, 1, 0]
        determinants.append(value % PRIME)
    skew = []
    for n in range(order + 1):
        value = sum(
            (-1) ** j
            * (
                determinants[j][0] * determinants[n - j][1]
                - determinants[j][1] * determinants[n - j][0]
            )
            for j in range(n + 1)
        )
        skew.append(value % PRIME)

    def first(values):
        return next((i, int(v)) for i, v in enumerate(values) if v)

    return {
        "left_nonzero": first([v[0] for v in determinants]),
        "right_nonzero": first([v[1] for v in determinants]),
        "relative_imaginary_nonzero": first(skew),
    }


if __name__ == "__main__":
    for strength in (8, 4):
        print(strength, certificate(strength))
