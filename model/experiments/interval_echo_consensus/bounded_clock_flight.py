"""Chebyshev propagation at integer exposures with bounded coefficient and recurrence errors."""
from functools import lru_cache
from math import ceil, factorial, isqrt

import numpy as np
from scipy import sparse


@lru_cache(maxsize=None)
def coefficients(x):
    """Enclose integer-argument Bessel coefficients by fixed-point series from their largest term."""
    if type(x) is not int or x <= 0:
        raise ValueError("Expected a positive integer argument")
    unit = 1 << 96
    result = []
    errors = []
    # Poisson's bound |J_n(x)| <= (x/2)^n/n!; subsequent terms decrease geometrically.
    n = x
    while 4 * x**n * (1 << 55) > (1 << n)*factorial(n):
        n += 1
    for k in range(n):
        peak = (isqrt(k*k + x*x) - k)//2
        order = k + 2*peak
        largest = unit*x**order // ((1 << order)*factorial(peak)*factorial(k+peak))
        total = (-1)**peak * largest
        term = largest
        for j in range(peak, 0, -1):
            term = term*(4*j*(k+j)) // (x*x)
            total += (-1)**(j-1)*term
        term = largest
        j = peak
        while term:
            j += 1
            term = term*(x*x) // (4*j*(k+j))
            total += (-1)**j*term
        # Both outward recurrences have ratios <= 1. Each floor costs < 1
        # fixed-point unit. The alternating omitted tail is at most its first term.
        errors.append((j+3)**2 / unit + 2**-52)
        weight = 1 if k == 0 else 2
        result.append(weight * (total / unit) * (-1j)**(k % 4))
    return np.array(result), 2*sum(errors) + 2**-55


class DetectorEvolution:
    """Only the requested detector coordinates are accumulated, without projecting the flight."""

    def __init__(self, h, times):
        if (h != h.T).nnz or not np.array_equal(2*h.data, np.rint(2*h.data)):
            raise ValueError("Expected a real symmetric half-integer generator")
        diagonal = h.diagonal()
        off = np.asarray(abs(h).sum(axis=1)).ravel() - abs(diagonal)
        low, high = min(diagonal-off), max(diagonal+off)
        center = int((low+high)//2)
        radius = 1 << (max(1, ceil(max(high-center, center-low))) - 1).bit_length()
        self.scaled = ((h - center*sparse.eye(h.shape[0], format="csr")) / radius).tocsr()
        if np.asarray(abs(self.scaled).sum(axis=1)).max() > 1:
            raise ArithmeticError("The scaled absolute matrix is not a contraction")
        series = [coefficients(radius*t) for t in times]
        count = max(len(c) for c, _ in series)
        self.weights = np.zeros((len(times), count), dtype=complex)
        # All matrices are dyadic. A Hermitian contraction has ||U_j|| <= j+1;
        # the same norm bound on |S| also controls sparse dot-product rounding.
        u = 2.**-53
        m = int(np.diff(self.scaled.indptr).max()) + 4
        gamma = m*u/(1-m*u)
        recurrence = 4*gamma*(count+1)**2
        if recurrence >= 1:
            raise ArithmeticError("The rounding induction cannot establish bounded iterates")
        summation = (count+3)*u/(1-(count+3)*u)
        self.probability_error = 0.
        self.amplitude_error = 0.
        for i, (c, coeff_error) in enumerate(series):
            self.weights[i,:len(c)] = c
            # Margin covers the final scalar arithmetic, magnitude squares and
            # summing 30 bins. Twice the home error covers the quiet complement.
            delta = coeff_error + (recurrence + 4*summation)*sum(abs(c)) + 1.e-13
            self.amplitude_error = max(self.amplitude_error, delta)
            self.probability_error = max(self.probability_error, 2*delta*(2+delta) + 1.e-12)
        if self.probability_error >= 1.e-5:
            raise ArithmeticError("The terminal probability error exceeds the registered budget")

    def amplitudes(self, launch, detector):
        """Compressed amplitudes of the centered generator, in one common time gauge."""
        state = np.zeros(self.scaled.shape[0])
        state[launch] = 1.
        selected = np.empty((self.weights.shape[1], len(detector)))
        selected[0] = state[detector]
        previous, state = state, self.scaled @ state
        selected[1] = state[detector]
        for k in range(2, len(selected)):
            previous, state = state, 2*(self.scaled @ state)-previous
            selected[k] = state[detector]
        return self.weights @ selected

    def read(self, launch, detector):
        # The subtracted scalar centre has no effect on any detector probability.
        ports = abs(self.amplitudes(launch, detector))**2
        p = np.column_stack((ports, 1-ports.sum(axis=1)))
        if not np.isfinite(p).all() or p.min() < 0 or p.max() > 1:
            raise ArithmeticError("Invalid terminal distribution")
        return p
