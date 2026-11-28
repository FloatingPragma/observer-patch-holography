"""Exact open-route rest gap and absence of a rest-frequency Bohr pair."""

from fractions import Fraction
from math import isqrt


def below(diagonal, rate, target):
    """Count eigenvalues of a Jacobi matrix strictly below an exact probe."""
    old, prior = Fraction(0), Fraction(1)
    changes = 0
    for entry in diagonal:
        value = (Fraction(entry) - target) * prior - rate**2 * old
        assert value
        changes += int(value * prior < 0)
        old, prior = prior, value
    return changes


def characteristic(diagonal, rate, prime):
    # Ascending coefficient order for det(x I - H).
    old, last = [], [1]
    for d in diagonal:
        out = [0, *last]
        for j, a in enumerate(last):
            out[j] = (out[j] - d * a) % prime
        for j, a in enumerate(old):
            out[j] = (out[j] - rate**2 * a) % prime
        old, last = last, out
    assert last[-1] == 1
    return last


def shifted(f, offset, prime):
    g = []
    for a in reversed(f):
        out = [0, *g]
        for j, b in enumerate(g):
            out[j] = (out[j] + offset * b) % prime
        out[0] = (out[0] + a) % prime
        g = out
    assert len(g) == len(f) and g[-1] == 1
    return g


def modular_gcd_degrees(f, g, prime):
    degrees = []
    while g:
        remainder = f.copy()
        factor = pow(g[-1], -1, prime)
        while remainder and len(remainder) >= len(g):
            shift = len(remainder) - len(g)
            scale = (remainder[-1] * factor) % prime
            for j, a in enumerate(g):
                remainder[shift + j] = (remainder[shift + j] - scale * a) % prime
            while remainder and remainder[-1] == 0:
                remainder.pop()
        f, g = g, remainder
        degrees.append(len(f) - 1)
    return degrees


def certificate(diagonal, rate):
    assert len(diagonal) == 144 and rate == 4
    assert diagonal == [4, 16, *([8, 8, 16] * 47), 4]
    # E*=10-2 sqrt(17), and the first rational cutoff sits at least .0016 higher.
    rest_hi, threshold = Fraction("1.7538"), Fraction("1.7554")
    rest_lo = Fraction("1.7537")
    assert ((10 - rest_hi) / 2) ** 2 < 17 < ((10 - rest_lo) / 2) ** 2
    assert below(diagonal, rate, rest_lo) == 0
    assert below(diagonal, rate, threshold) == 1
    assert threshold - rest_hi == Fraction(1, 625)
    # The coefficient ring is Z[sqrt(17)]. A nonzero
    # resultant modulo the homomorphism sqrt(17) -> s rules out a zero
    # resultant in that ring and hence any real Bohr pair at +/- E*.
    for p in (151, 179, 191, 223, 229, 239, 251, 257, 263):
        assert all(p % d for d in range(2, isqrt(p) + 1))
        poly = characteristic(diagonal, rate, p)
        for s in range(p):
            if s**2 % p != 17:
                continue
            degrees = modular_gcd_degrees(poly, shifted(poly, (10 - 2 * s) % p, p), p)
            if degrees[-1] == 0:
                return {
                    "prime": p,
                    "sqrt17_image": s,
                    "euclidean_divisor_degrees": degrees,
                    "eigenvalues_below_1_7554": 1,
                    "rest_gap_lower": 1 / 625,
                }
    raise ValueError("No nonzero resultant certificate from the tested primes")
