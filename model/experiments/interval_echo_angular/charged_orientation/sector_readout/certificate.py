"""Exact channel identities at repeated roots of integer charge Hamiltonians."""

from itertools import combinations
import numpy as np
from sympy import Matrix, Poly, oo
from sympy.polys.matrices import DomainMatrix


def commutant_dimension(h, locations):
    # Commuting with every local position projector first eliminates all
    # off-block entries. Solve [H,X]=0 for every surviving matrix unit.
    generators = []
    for i, left in enumerate(locations):
        for j, right in enumerate(locations):
            if left == right:
                unit = np.zeros_like(h)
                unit[i, j] = 1
                generators.append((h @ unit - unit @ h).ravel())
    matrix = DomainMatrix.from_Matrix(Matrix(np.column_stack(generators).tolist()))
    return len(generators) - matrix.rank()


def reduce_polynomial(coeff, monic):
    result = list(coeff)
    degree = len(monic) - 1
    for high in reversed(range(degree, len(result))):
        for j in range(degree):
            result[high - degree + j] -= result[high] * monic[j]
    return result[:degree]


def residue_matrices(h, minimal, factor):
    """At a root x of factor, these coefficients evaluate to minimal'(x) P_x."""
    n = len(h)
    ident = np.eye(n, dtype=object)
    # The divided difference (m(y)-m(x))/(y-x) at y=H. The coefficient
    # of x^k is obtained by Horner evaluation on H.
    term = ident
    coeff = [term]
    for entry in minimal.all_coeffs()[1:-1]:
        term = h @ term + int(entry) * ident
        coeff.append(term)
    monic = list(map(int, reversed(factor.all_coeffs())))
    return reduce_polynomial(coeff[::-1], monic), monic


def certify_poles(h, inputs, destinations):
    polynomial = Poly(Matrix(h.tolist()).charpoly().as_expr())
    minimal = polynomial.sqf_part()
    # Check the square-free annihilator against the actual integer H as well.
    obj = h.astype(object)
    ident = np.eye(len(h), dtype=object)
    rem = ident
    for c in minimal.all_coeffs()[1:]:
        rem = obj @ rem + int(c) * ident
    assert not rem.any()
    multiplicities = {}
    checked = 0
    for factor, mult in polynomial.factor_list()[1]:
        # H is integer 8 times the physical formation operator.
        positive = int(factor.count_roots(0, oo))
        multiplicities[mult] = multiplicities.get(mult, 0) + positive
        if mult == 1 or not positive:
            continue
        entries, monic = residue_matrices(obj, minimal, factor)
        # J equals its partial transpose iff every summed 2x2 Kraus minor
        # vanishes. Check every coefficient exactly in Q[x]/factor.
        kraus = np.stack(
            [[mat[np.ix_(idx, inputs)] for mat in entries] for idx in destinations]
        )
        degree = len(entries)
        for r, s in combinations(range(len(destinations[0])), 2):
            for i, j in combinations(range(len(inputs)), 2):
                minor = [0] * (2 * degree - 1)
                for first in range(degree):
                    for second in range(degree):
                        minor[first + second] += sum(
                            k[first, r, i] * k[second, s, j]
                            - k[first, r, j] * k[second, s, i]
                            for k in kraus
                        )
                assert not any(reduce_polynomial(minor, monic))
                checked += 1
    return multiplicities, checked
