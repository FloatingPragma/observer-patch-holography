"""Exact Taylor coefficients of the dimensionless regional modular energy on the cell."""

from collections import defaultdict
from fractions import Fraction
from math import factorial

import numpy as np

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import cell


def scaled_native_contact(strength):
    """Read and exactly integerize four times the trace-centered binary contact."""
    ids = [native_law.SRC.index(w) for w in ((0, 0), (1, 0), (0, 1), (1, 1))]
    h = native_law.seam_term(strength)[np.ix_(ids, ids)]
    centered = 4 * h - np.trace(h) * np.eye(4)
    h4 = centered.astype(int)
    np.testing.assert_array_equal(centered, h4)
    return h4


def integer_contact(strength):
    """Return coefficients in the real basis XX, (XZ)(XZ), ZZ."""
    h4 = scaled_native_contact(strength)
    x = np.array([[0, 1], [1, 0]])
    z = np.diag([1, -1])
    pairs = [np.kron(p, p) for p in (x, x @ z, z)]
    components = [int(np.trace(p.T @ h4)) // 4 for p in pairs]
    np.testing.assert_array_equal(sum(c * p for c, p in zip(components, pairs)), h4)
    return components


def hamiltonian(edges, components):
    terms = defaultdict(int)
    for left, right in edges:
        ports = (1 << left) | (1 << right)
        for key, weight in zip(((ports, 0), (ports, ports), (0, ports)), components):
            terms[key] += weight
    return dict(terms)


def multiply(left, right):
    """Real Pauli basis X^x Z^z. Moving Z^z past X^u gives (-1)^popcount(z & u)."""
    out = defaultdict(int)
    for (x, z), a in left.items():
        for (u, v), b in right.items():
            out[x ^ u, z ^ v] += -a*b if (z & u).bit_count() % 2 else a*b
    return {key: value for key, value in out.items() if value}


def coefficient_of_product(left, right, key):
    tx, tz = key
    total = 0
    for (x, z), a in left.items():
        u = x ^ tx
        b = right.get((u, z ^ tz), 0)
        total += -a*b if (z & u).bit_count() % 2 else a*b
    return total


def conditional_expectation(poly, outside):
    return {(x, z): value for (x, z), value in poly.items() if not (x | z) & outside}


def regional_pair_series(edges, cut, read_pair, strength):
    """Coefficients of beta^1, beta^2, beta^3 in -log E exp(-beta H), in XX,YY,ZZ order."""
    h = hamiltonian(edges, integer_contact(strength))
    h2 = multiply(h, h)
    e1 = conditional_expectation(h, cut)
    e2 = conditional_expectation(h2, cut)
    e1_squared = multiply(e1, e1)
    p = sum(1 << v for v in read_pair)
    if p & cut:
        raise ValueError("the measured pair must lie inside the region")
    series = [[], [], []]
    # YY = -(XZ)(XZ). The stored H is four times the physical centered Hamiltonian.
    for key, sign in (((p, 0), 1), ((p, p), -1), ((0, p), 1)):
        e1e2 = coefficient_of_product(e1, e2, key)
        e2e1 = coefficient_of_product(e2, e1, key)
        e1_cube = coefficient_of_product(e1_squared, e1, key)
        e3 = coefficient_of_product(h2, h, key)
        series[0].append(Fraction(sign * e1.get(key, 0), 4))
        series[1].append(Fraction(sign * (e1_squared.get(key, 0) - e2.get(key, 0)), 2 * 4**2))
        series[2].append(Fraction(sign * (2*e3 - 3*(e1e2 + e2e1) + 4*e1_cube), 12 * 4**3))
    return series


def run(strength=KAPPA):
    board = cell()
    path = [(0, 1), (1, 8)]
    assert all(edge in board.seams for edge in path) and 8 not in board.neighbours[0]
    return {label: regional_pair_series(edges, 1 << 0, (1, 8), strength)
            for label, edges in (("two selected contacts", path), ("all cell contacts", board.seams))}


def path_gap_at_readout(strength=KAPPA):
    """Certified lower bound for exp(beta) f_odd - f_even on the two selected contacts."""
    beta = Fraction(1, 16)
    degree = 12
    # The 4-scaled trace-centered contact has spectral endpoints -(2*k+1), 2*k+3.
    # Subtracting I centers this interval for every nonnegative integer contrast.
    scaled_bond = scaled_native_contact(strength) - np.eye(4, dtype=int)
    scaled_path = (np.kron(scaled_bond, np.eye(2, dtype=int))
                   + np.kron(np.eye(2, dtype=int), scaled_bond)).astype(object)
    radius = beta * (Fraction(strength) + 1)
    if strength < 0 or radius >= degree + 2:
        raise ValueError("the certificate needs nonnegative contrast and a geometric tail ratio below one")
    term = np.eye(8, dtype=int).astype(object)
    polynomial = term.copy()
    for order in range(1, degree + 1):
        term = (term @ scaled_path) * (-beta / (4 * order))
        polynomial += term
    reduced = (polynomial[:4, :4] + polynomial[4:, 4:]) * Fraction(1, 2)
    diagonal = (reduced[1, 1] + reduced[2, 2]) / 2
    exchange = (reduced[1, 2] + reduced[2, 1]) / 2
    # The ratio of successive terms in the exponential tail is <= radius/(degree+2).
    remainder = radius**(degree + 1) / factorial(degree + 1) / (1 - radius / (degree + 2))
    odd_lower = diagonal - exchange - remainder
    if odd_lower <= 0:
        raise ValueError("the odd eigenvalue lower bound must be positive")
    scalar_exp_lower = sum(beta**k / factorial(k) for k in range(5))
    return scalar_exp_lower * odd_lower - (diagonal + exchange + remainder)


if __name__ == "__main__":
    for kappa in (KAPPA, KAPPA / 2):
        for label, series in run(kappa).items():
            print(f"contrast {kappa:g}, {label}")
            for order, projected in enumerate(series, 1):
                print(f"order {order}: " + ", ".join(map(str, projected)))
        print(f"two-contact fixed-temperature gap lower bound {float(path_gap_at_readout(kappa)):.10g}")
