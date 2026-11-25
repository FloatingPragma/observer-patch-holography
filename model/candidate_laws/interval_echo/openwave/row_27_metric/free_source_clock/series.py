"""Exact short-exposure coefficients of the paired-clock discrepancy."""

from fractions import Fraction
from math import comb, factorial

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .measure import ENDPOINTS


def coherence_coefficients(site, gap, order):
    """Coefficients of even powers of t in the coin expectation 1-2 p_flip."""
    board = federation(0)
    ends = set(ENDPOINTS)
    neighbours = [tuple(v for v in row if {u, v} != ends)
                  for u, row in enumerate(board.neighbours)]
    k = int(candidate.KAPPA)
    assert k == candidate.KAPPA
    # Twice the ordinary-contact kinetic energy, contrast and degree cost.
    diagonal = [(1 + 2 * k) * len(row) + 6 * k * (6 - int(board.degree[u]))
                for u, row in enumerate(neighbours)]

    def powers(antisymmetric):
        result = [[int(u == site) for u in range(board.n)]]
        shifts = [2 * gap if antisymmetric or u in ends else 0 for u in range(board.n)]
        for _ in range(order):
            previous = result[-1]
            result.append([(diagonal[u] + shifts[u]) * previous[u] - sum(previous[v] for v in row)
                           for u, row in enumerate(neighbours)])
        return result

    anti, sym = powers(True), powers(False)
    coeff = []
    for n in range(0, order + 1, 2):
        # Expand Re <x|exp(i A t) exp(-i S t)|x>; the two factors need
        # not commute. All powers above are of the doubled generators.
        tot = 0
        for j in range(n + 1):
            inner = sum(v * w for v, w in zip(anti[j], sym[n - j]))
            tot += (-1) ** j * comb(n, j) * inner
        coeff.append(Fraction((-1) ** (n // 2) * tot, 2**n * factorial(n)))
    return coeff


def chebyshev_coefficients(order, cosine):
    """Truncated expansion about cosine[0] = 1, in the variable t squared."""
    assert cosine[0] == 1
    delta = [0, *cosine[1:]]
    power = [1] + [0] * (len(delta) - 1)
    result = power.copy()
    weight = Fraction(1)
    for k in range(1, len(delta)):
        # T_n^(k)(1)/k! follows by differentiating the Chebyshev equation.
        weight *= Fraction(order**2 - (k - 1)**2, (2 * k - 1) * k)
        power = [sum(power[i] * delta[j-i] for i in range(j + 1)) for j in range(len(delta))]
        result = [v + weight * p for v, p in zip(result, power)]
    return result


def first_difference(site, order):
    even = chebyshev_coefficients(33, coherence_coefficients(site, 32, order))
    odd = chebyshev_coefficients(32, coherence_coefficients(site, 33, order))
    return next(((2 * j, a - b) for j, (a, b) in enumerate(zip(even, odd)) if a != b), None)


if __name__ == "__main__":
    for site in (0, 1, 2, 13, 156):
        print(site, first_difference(site, 32))
