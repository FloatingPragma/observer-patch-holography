"""Certify geometric curvature for every density matrix of the enrolled radar survey."""

from itertools import combinations, product

import numpy as np
import sympy as sp
from mpmath import iv

from .station import ORDER

# Endpoints are in ten-thousandths of one range unit, in directed station order.
BOUNDS = (
    (25004, 27568), (29958, 30000), (29998, 30000), (25017, 27574),
    (28644, 29452), (29962, 30000), (29958, 30000), (28644, 29452),
    (26525, 27574), (29998, 30000), (29962, 30000), (26519, 27568),
)
SPLIT = (0, 3, 4, 7, 8, 11)
PRECISION = 2**48


def positive_definite(matrix):
    return all(sp.expand_complex(matrix[:k, :k].det()).as_real_imag()[0] > 0
               for k in range(1, matrix.rows + 1))


def certify_ranges(effects, error):
    for pair, (lower, upper) in zip(effects, BOUNDS, strict=True):
        for block in pair:
            assert block == block.H
            assert positive_definite(block - (sp.Rational(lower, 10000) + error) * sp.eye(4))
            assert positive_definite((sp.Rational(upper, 10000) - error) * sp.eye(4) - block)


def squared_measures():
    """Heron and Gram polynomials in six squared edges, pulled back to log ranges."""
    edges = list(combinations(range(4), 2))
    variables = sp.symbols("d:6")
    lengths = dict(zip(edges, variables))
    directed = list(zip(*ORDER))
    powers = np.zeros((6, 12), dtype=int)
    for row, (a, b) in zip(powers, edges):
        row[directed.index((a, b))] += 4
        row[directed.index((b, a))] += 4
        for vertex in (a, b):
            for other in range(4):
                if other != vertex:
                    row[directed.index((vertex, other))] -= 1
                    row[directed.index((other, vertex))] += 1

    def gram(a, b, c):
        ab = lengths[tuple(sorted((a, b)))]
        ac = lengths[tuple(sorted((a, c)))]
        return (ab + ac - lengths[tuple(sorted((b, c)))]) / 2

    polynomials = [(lengths[a, b] * lengths[a, c] - gram(a, b, c)**2) / 4
                   for a, b, c in combinations(range(4), 3)]
    matrix = sp.Matrix(3, 3, lambda i, j: lengths[0, i+1] if i == j else gram(0, i+1, j+1))
    polynomials.append(matrix.det() / 36)
    measures = []
    for polynomial in polynomials:
        terms = sp.Poly(polynomial, variables).terms()
        # Each monomial is coefficient * exp(k.log(r)/4).
        k = np.array([np.array(exponents) @ powers for exponents, _ in terms], dtype=object)
        coefficients = [iv.mpf(int(c.p)) / int(c.q) for _, c in terms]
        # Value, log gradient and log Hessian minus diagonal log gradient, all times 16.
        table = np.array([np.r_[16, 4*row, (np.outer(row, row)-4*np.diag(row)).ravel()]
                          for row in k], dtype=object)
        measures.append((coefficients, k, table))
    return measures


def upper_units(value):
    """An integer majorant in fixed point for a nonnegative interval."""
    assert value.a >= 0
    return int((value * PRECISION).b) + 1


def polynomial_jet(measure, log_midpoints, log_widths):
    coefficients, powers, table = measure
    nominal, error = [], []
    for coefficient, row in zip(coefficients, powers):
        value = coefficient * iv.exp(sum(int(k)*x for k, x in zip(row, log_midpoints)) / 4)
        h = int((value * PRECISION).a)
        rounded = iv.mpf(h) / PRECISION
        variation = iv.mpf(int(abs(row) @ log_widths)) / (4 * PRECISION)
        tail = iv.exp(variation) - 1 - variation
        # The centered Taylor bound retains cancellations between monomials in the linear term.
        remainder = abs(value) * abs(tail) + abs(value - rounded) * (1 + variation)
        nominal.append(h)
        error.append(upper_units(remainder))
    nominal = np.array(nominal, dtype=object)
    centers = nominal @ table
    linear = abs((table.T * nominal) @ powers) @ log_widths
    errors = np.array(error, dtype=object) @ abs(table)
    jets = []
    for center, residual, slope in zip(centers, errors, linear):
        interval = iv.mpf([int(center-residual), int(center+residual)]) / (16 * PRECISION)
        linear_radius = iv.mpf(int(slope)) / (64 * PRECISION**2)
        jets.append(interval + iv.mpf([-1, 1]) * linear_radius)
    return jets[0], np.array(jets[1:13], dtype=object), np.array(jets[13:], dtype=object).reshape(12, 12)


def box_curvature(measures, choices):
    # Each endpoint is an exact multiple of 1/20000, including split midpoints.
    bounds = [(2*a, 2*b) for a, b in BOUNDS]
    for i, upper_half in zip(SPLIT, choices):
        a, b = BOUNDS[i]
        bounds[i] = (a+b, 2*b) if upper_half else (2*a, a+b)
    log_midpoints = [iv.ln(iv.mpf(a+b)/40000) for a, b in bounds]
    log_widths = np.array([upper_units(iv.ln(iv.mpf(a+b)/(2*a))) for a, b in bounds], dtype=object)
    # Tr(D)=0 gives |Tr(E_i D)| <= (upper_i-lower_i) ||D||_1 / 2.
    weights = [iv.mpf(hi-lo) / a for (lo, hi), (a, _) in zip(BOUNDS, bounds)]
    gains = [iv.mpf(['23.1042361', '23.1042363'])]*4 + [iv.mpf(['-51.20', '-50.98'])]
    hessian = np.full((12, 12), iv.mpf(0), dtype=object)
    for measure, gain in zip(measures, gains):
        value, gradient, second = polynomial_jet(measure, log_midpoints, log_widths)
        assert value.a > 0
        root = iv.sqrt(value)
        hessian += (second - np.outer(gradient, gradient) * (1 / (2*value))) * (gain / (2*root))
    return sum(abs(hessian[i, j]) * weights[i] * weights[j] for i in range(12) for j in range(12))


def verify(effects, error):
    iv.dps = 30
    certify_ranges(effects, error)
    measures = squared_measures()
    boxes = [box_curvature(measures, choices) for choices in product((False, True), repeat=len(SPLIT))]
    largest = max(bound.b for bound in boxes)
    assert largest < iv.mpf('.72')
    return dict(boxes=len(boxes), curvature_upper_millionths=int((largest * 10**6).b) + 1,
                trace_norm_curvature_limit=0.72, feasible_objective_gap_coefficient=0.14)
