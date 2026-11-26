"""Integer enclosure of the two native probabilities at six clock radians."""

from fractions import Fraction
from math import factorial

import numpy as np

# Ports of the two-dimer arms, in the orientation order of the load read.
DIMERS = ((6, 7), (1, 82), (74, 75), (3, 48), (49, 53), (10, 20), (13, 12))
CONTACTS = ((7, 1), (82, 74), (6, 3), (48, 49), (6, 10), (20, 13))
SCALE = 2**104
ORDER = 40
SLICES = 192


def orientation_data():
    # Independent construction on the named native binary endpoint loads.
    addresses = {
        port: (i, bit)
        for i, dimer in enumerate(DIMERS)
        for bit, port in enumerate(dimer)
    }
    diagonal, walls = [], []
    for word in range(128):
        occupied = {p: (word >> i & 1) == bit for p, (i, bit) in addresses.items()}
        energy = 32 * sum(not (occupied[a] and occupied[b]) for a, b in CONTACTS)
        diagonal.append(2 * energy - 192)
        walls.append(
            sum((occupied[a] != occupied[b]) << j for j, (a, b) in enumerate(CONTACTS))
        )
    return np.array(diagonal, dtype=object), np.array(walls)


def remote(phase):
    """Lower and upper bounds on P(load(82) differs from load(74))."""
    if phase not in (1, -1):
        raise ValueError("Use the even or odd source")
    diagonal, walls = orientation_data()
    flips = np.arange(128)[None, :] ^ (1 << np.arange(7))[:, None]
    # M=2H-192I; symmetry bounds its spectral norm by its largest row sum.
    assert max(abs(diagonal)) + 7 <= 256

    def apply(vector):
        return diagonal * vector - vector[flips].sum(axis=0)

    real = np.zeros(128, dtype=object)
    imag = real.copy()
    real[(walls == 4) & (np.arange(128) % 2 == 1)] = SCALE
    real[(walls == 16) & (np.arange(128) % 2 == 1)] = phase * SCALE
    assert sum(real * real) == 2 * SCALE * SCALE
    # One step lasts 1/32. Evaluate the Taylor polynomial of -iM/64
    # by Horner evaluation, rounding each component toward minus infinity.
    for _ in range(SLICES):
        re, im = real, imag
        for j in range(ORDER, 0, -1):
            re, im = real + apply(im) // (64 * j), imag + (-apply(re)) // (64 * j)
        real, imag = re, im
    # A complex vector rounding has norm <=16/SCALE. Its propagation through
    # Horner has gain sum_k 4^k/k! <55. Exact evolution has norm one.
    step_error = 55 * Fraction(16, SCALE) + 2 * Fraction(
        4 ** (ORDER + 1), factorial(ORDER + 1)
    )
    state_error = SLICES * step_error
    assert state_error < Fraction(1, 10**15)
    selected = (walls & 2) != 0
    p = Fraction(sum(real[selected] ** 2 + imag[selected] ** 2), 2 * SCALE**2)
    # Both the unnormalized true vector (norm sqrt(2)) and approximation
    # have norm below two, so their squared norms divided by two differ by <=2e.
    return p - 2 * state_error, p + 2 * state_error


if __name__ == "__main__":
    for name, phase in (("plus", 1), ("minus", -1)):
        lower, upper = remote(phase)
        denominator = 10**12
        print(
            name,
            lower.numerator * denominator // lower.denominator,
            -((-upper.numerator * denominator) // upper.denominator),
            "in units of 10^-12",
        )
