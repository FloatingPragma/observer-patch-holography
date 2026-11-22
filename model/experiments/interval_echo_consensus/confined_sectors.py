"""Classify one-cell confinement in the one-load and distinct two-load blank sectors."""
from collections import Counter

import numpy as np
import sympy as sp

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model import cell, federation
from model.executable_model.carrier import automorphisms


def integer_matrix(a):
    if not np.array_equal(a, np.rint(a)):
        raise ArithmeticError("Expected an integer carrier operator")
    return sp.Matrix(a.astype(int))


def boundary_classes():
    actions = automorphisms(cell().seams)
    counts = {}
    for level in range(3):
        b = federation(level)
        gates = [set() for _ in range(b.cells)]
        inside = [[] for _ in range(b.cells)]
        for a, d in b.seams:
            if a//12 == d//12:
                inside[a//12].append((a % 12, d % 12))
        if any(edges != cell().seams for edges in inside):
            raise ValueError("A federation cell differs from the primitive")
        seen = set()
        for edge in b.external:
            for v in edge:
                if v in seen:
                    raise ValueError("The external seams must be a matching")
                seen.add(v)
                gates[v//12].add(v % 12)
        if any(len(g) != 3 for g in gates):
            raise ValueError("Expected three external ports per cell")
        counts[b.cells] = Counter(min(tuple(sorted(a[v] for v in g)) for a in actions) for g in gates)
    return counts


def confined_at(generator, boundary):
    """Unobservable subspace for an exact rational generator and a boundary coordinate read."""
    rows = boundary
    for _ in range(generator.rows):
        vectors = rows.nullspace()
        if not vectors:
            return sp.zeros(generator.rows, 0)
        frame = sp.Matrix.hstack(*vectors)
        if (rows * generator * frame).is_zero_matrix:
            return frame
        # Each pass retains a row basis, avoiding an exponentially growing power stack.
        enlarged = rows.col_join(rows * generator)
        rows = enlarged[list(enlarged.T.rref()[1]), :]
    raise ArithmeticError("The invariant-subspace iteration did not terminate")


def possible_exceptions(kinetic, adjacent, boundary):
    """Every exceptional strength is a common zero of the exhibited polynomial minors."""
    x = sp.Symbol("x")
    frame = sp.Matrix.hstack(*boundary.col_join(boundary*kinetic).nullspace())
    # B P V = B K V = 0. A vector confined under A=K+xP starts in im(V).
    coefficients = [frame]
    rows = []
    for order in range(1, 5):
        terms = [kinetic*v for v in coefficients] + [sp.zeros(*frame.shape)]
        for d, v in enumerate(coefficients):
            terms[d+1] += adjacent*v
        coefficients = terms
        if order >= 2:
            projected = [boundary*v for v in terms]
            for r in range(boundary.rows):
                row = [sum(x**d * v[r, j] for d, v in enumerate(projected) if v[r, j])
                       for j in range(frame.cols)]
                if any(row):
                    rows.append(row)
    obs = sp.Matrix(rows)
    divisor = sp.Poly(0, x)
    # Samples select square minors only. Their exact polynomial determinants give the proof.
    for value in (0, 2, 1, -1, 3):
        for matrix in (obs, obs[::-1, :]):
            pivots = matrix.subs(x, value).T.rref()[1]
            if len(pivots) != frame.cols:
                continue
            minor = matrix[list(pivots), :].to_DM()
            polynomial = sp.Poly(minor.domain.to_sympy(minor.det()), x)
            divisor = sp.gcd(divisor, polynomial).monic()
            if divisor.degree() == 0 or divisor.sqf_part() == sp.Poly(x-2, x, domain=sp.QQ):
                return divisor.as_expr()
    raise ArithmeticError(f"Unresolved roots of the confinement minors: {divisor.as_expr()}")


def classify():
    b = cell()
    sectors = PairSector(b)
    lap = integer_matrix(b.laplacian().toarray())
    kinetic = integer_matrix(2*sectors.exchange().toarray())
    adjacent = sp.diag(*map(int, sectors.adjacent))
    census = boundary_classes()
    results = []
    for gates in sorted(set().union(*census.values())):
        singles = confined_at(lap, sp.eye(b.n)[list(gates), :])
        boundary_ids = np.flatnonzero(np.isin(sectors.aa, gates) | np.isin(sectors.bb, gates))
        boundary = sp.eye(sectors.P)[boundary_ids.tolist(), :]
        divisor = possible_exceptions(kinetic, adjacent, boundary)
        # The only remaining point, if any, is x=2, or kappa=-1/2.
        pairs = confined_at(kinetic+2*adjacent, boundary)
        if pairs.cols:
            if ((kinetic+2*adjacent)*pairs != 12*pairs
                    or 2*pairs.T*adjacent*pairs != pairs.T*pairs):
                raise ArithmeticError("Exceptional energies or adjacent-pair weights differ")
        results.append({"gateways": gates, "cells": {n: counts[gates] for n, counts in census.items()},
                        "one_load_dimension": singles.cols, "necessary_polynomial_in_minus_4k": str(divisor),
                        "pair_dimension_at_minus_half": pairs.cols})
    return results


if __name__ == "__main__":
    import json
    print(json.dumps(classify(), indent=2))
