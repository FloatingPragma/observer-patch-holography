"""Exact upper bound for every internal-seam preparation in the M10.3 exchange control."""
from fractions import Fraction
from math import factorial

import numpy as np
from scipy import sparse
from sympy import integer_nthroot

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model import federation


DEGREE = 19


def compressed_polynomial(exchange, indices):
    """Numerator of the compressed Taylor polynomial for exp(-i(K-6I))."""
    size = exchange.shape[0]
    if (exchange.shape != (size, size) or len(indices) == 0 or len(set(indices)) != len(indices)
            or any(not isinstance(j, (int, np.integer)) or not 0 <= j < size for j in indices)):
        raise ValueError("Expected distinct coordinates of a square generator")
    twice = (2 * exchange).tocsr()
    if (twice - twice.T).nnz or not np.array_equal(twice.data, np.rint(twice.data)):
        raise ValueError("Expected symmetric half-integer exchange")
    b = twice.astype(np.int64) - 12 * sparse.eye(size, dtype=np.int64, format="csr")
    if abs(b).sum(axis=1).max() > 12:
        raise ArithmeticError("The exchange spectrum exceeds the Taylor interval")
    power = np.zeros((size, len(indices)), dtype=np.int64)
    power[indices] = np.eye(len(indices), dtype=np.int64)
    denominator = 2**DEGREE * factorial(DEGREE)
    real = np.zeros((len(indices), len(indices)), dtype=object)
    imag = real.copy()
    for order in range(DEGREE + 1):
        weight = (-1)**(order//2) * (denominator // (2**order * factorial(order)))
        part = imag if order % 2 else real
        part += (-weight if order % 2 else weight) * power[indices].astype(object)
        if order < DEGREE:
            # Bound every intermediate sparse dot product before int64 multiplication.
            peak = max(int(power.max()), -int(power.min()))
            if 12*peak > np.iinfo(np.int64).max:
                raise ArithmeticError("The integer recurrence could overflow")
            power = b @ power
    return real, imag, denominator


def home_bound(sector, home):
    indices = [sector.index(a, b) for a, b in sector.car.seams if a//12 == b//12 == home]
    if len(set(indices)) != 30:
        raise ValueError("Incomplete cell-13 detector")
    real, imag, denom = compressed_polynomial(sector.exchange(), indices)
    # Each real and imaginary entry is rounded down by less than one unit.
    # The Frobenius error of the complex return block is < 2*30 / scale.
    scale = 2**20
    r, i = real * scale // denom, imag * scale // denom
    block = np.block([[r, -i], [i, r]])
    gram = block.T @ block
    exponent = 32
    gram_power = np.linalg.matrix_power(gram, exponent)
    norm_bound = max(abs(gram_power).sum(axis=1))
    root, exact = integer_nthroot(int(norm_bound), 2*exponent)
    cap = Fraction(root if exact else root+1, scale)
    # The integral Taylor remainder has a unitary factor for Hermitian K-6I.
    tail = Fraction(6**(DEGREE+1), factorial(DEGREE+1))
    upper = (cap + Fraction(2 * len(indices), scale) + tail)**2
    delta = Fraction(1, 50 if home == 13 else 100)
    # Up to a scalar, H_k-K=k(I-2P) has norm |k|. Duhamel at exposure one
    # adds at most delta to the return amplitude, for every |k| <= delta.
    margin = Fraction(94, 100) - upper - delta**2
    if margin <= 0 or margin**2 <= 4*upper*delta**2:
        raise ArithmeticError("The near-zero strengths cannot be excluded at this precision")
    return {"cell": home, "exposure": 1, "states": len(indices), "probability_upper": str(upper),
            "upper_decimal": float(upper), "excluded_absolute_strength_through": str(delta)}


if __name__ == "__main__":
    import json
    sector = PairSector(federation(0))
    print(json.dumps([home_bound(sector, home) for home in range(sector.car.cells)], indent=2))
