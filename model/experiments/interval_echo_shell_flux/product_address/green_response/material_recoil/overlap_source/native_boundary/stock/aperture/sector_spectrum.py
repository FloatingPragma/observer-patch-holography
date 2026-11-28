"""Identify the independent native calorimetric data retained by the modular reference."""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache

import numpy as np
from mpmath import iv
from sympy import Matrix, Poly, symbols

from . import peer_fibers as fibers
from . import sector_thermal as thermal


@lru_cache(maxsize=1)
def spectra():
    _, _, _, blocks, ops = thermal.apparatus()
    four_h = ops[0].sum(0)
    shapes = []
    moments = []
    t = symbols("t")
    for group in blocks:
        M = Matrix(four_h[np.ix_(group, group)].tolist()) / 4
        centered = M - M.trace() / len(group) * Matrix.eye(len(group))
        # Each irreducible rational factor fixes its algebraic roots.
        # Exponents per sharp word specify their normalized multiplicities.
        factors = Poly(centered.charpoly(t).as_expr(), t).factor_list()[1]
        shapes.append(
            tuple(
                sorted(
                    (tuple(map(Q, p.all_coeffs())), Q(k, len(group)))
                    for p, k in factors
                )
            )
        )
        moments.append(Q((centered * centered).trace() / len(group)))
    return shapes, moments


def design():
    w, c, ids, cat, _, n, *extra = fibers.apparatus()
    B = Matrix(extra[-1].tolist())
    shapes, moments = spectra()
    _, _, _, _, ops = thermal.apparatus()

    # Within a free category all protected histories and all the previous
    # observables are tested independently of the calorimetric classification.
    def means(feature):
        return [sum(feature[cat == i]) / Q(int(n[i])) for i in range(len(n))]

    strain = ops[1].diagonal(axis1=1, axis2=2)
    for i in range(len(n)):
        assert (
            len(
                {
                    tuple(Q(int(x.sum()), len(g)) for x in strain[:, g])
                    for g in thermal.apparatus()[3]
                    if cat[g[0]] == i
                }
            )
            == 1
        )
    constraints = Matrix(
        [
            means(np.array(x, dtype=int).astype(object))
            for x in [c == 0, w[:, -1] == 0, *strain]
        ]
    )
    V = B.T * Matrix.hstack(*(constraints * B.T).nullspace())
    types = sorted(set(shapes))
    sid = np.array([types.index(k) for k in shapes])[ids]
    var = np.array(moments, object)[ids]
    assert all(len(set(sid[cat == i])) == 1 for i in range(len(n)))
    D = Matrix(
        [means((sid == i).astype(int).astype(object)) for i in range(len(types))]
    )
    variance = Matrix([means(var)])
    # Operating constraints hold the native stiffness response fixed. Work
    # in the full space before the joint-cell and blank readings are fixed.
    strain_space = B.T * Matrix.hstack(*(constraints[2:, :] * B.T).nullspace())
    reads = constraints[:2, :].col_join(variance)
    return V, D, variance, shapes, sid, var, strain_space, reads


def thermal_price(shape, beta):
    """Enclose the centered thermal trace from the exact characteristic factors."""

    def fraction(x):
        return iv.mpf(x.numerator) / x.denominator

    total = iv.mpf(0)
    for coeff, prob in shape:
        c = list(map(fraction, coeff))
        if len(c) == 2:
            roots = [-c[1] / c[0]]
        else:
            assert len(c) == 3
            disc = iv.sqrt(c[1] ** 2 - 4 * c[0] * c[2])
            roots = [(-c[1] + sign * disc) / (2 * c[0]) for sign in (-1, 1)]
        total += fraction(prob) * sum(iv.exp(-beta * r) for r in roots)
    return iv.log(total)


def rational_lines(shape):
    out = {}
    for coeff, p in shape:
        assert len(coeff) == 2
        out[-coeff[1] / coeff[0]] = p
    return out


def convolution(a, b):
    out = Counter()
    for x, u in a.items():
        for y, v in b.items():
            out[x + y] += u * v
    return out


def reconstruction():
    V, D, variance, shapes, _, _, U, reads = design()
    shapes = sorted(set(shapes))
    assert len(shapes) == 8
    narrow, broad, mix1, mix2 = [rational_lines(shapes[i]) for i in (3, 2, 4, 5)]
    assert convolution(narrow, broad) == convolution(mix1, mix2)
    # Check the coefficients of Z_broad-Z_narrow against its exact
    # positive square: (exp(65 beta/4)-exp(-65 beta/4)cosh(beta/2))²/4.
    signed = Counter(
        {Q(-65, 4): Q(1), Q(65, 4) - Q(1, 2): -Q(1, 2), Q(65, 4) + Q(1, 2): -Q(1, 2)}
    )
    sq = convolution(signed, signed)
    assert all(
        broad.get(e, 0) - narrow.get(e, 0) == sq[e] / 4
        for e in set(broad) | set(narrow) | set(sq)
    )
    # Eliminate log Z_5=log Z_2+log Z_3-log Z_4, then determine the
    # coefficients of the three measured changes for every allowed source.
    reduced = D.copy()
    for i, sign in ((2, 1), (3, 1), (4, -1)):
        reduced[i, :] += sign * D[5, :]
    reduced[5, :] = Matrix.zeros(1, D.cols)
    X = reads * U
    assert X.rank() == 3
    C = reduced * U * X.T * (X * X.T).inv()
    assert C * X == reduced * U
    assert list(C[:, 2]) == [0, 0, Q(8, 4225), -Q(8, 4225), 0, 0, 0, 0]
    return C, V, reads
