"""Positive formation levels of the winding preparation with at most one blank contact."""

from collections import defaultdict
from fractions import Fraction

import numpy as np
from sympy import Matrix, Poly, Rational

from model.experiments.interval_echo_angular.charged_orientation.winding_memory import (
    measure as winding,
)

LAW = winding.LAW


def sector(sign, edges, private_port=None):
    """Start the charge with all direction coins zero and follow every native column."""
    terms = winding.blocks()
    source = (0, 0)
    words = [source]
    index = {source: 0}
    columns = []
    for port, heads in words:
        column = defaultdict(int)
        if private_port == port:
            column[port, heads] += 8 * (-8 if sign < 0 else 24)
        for left, right, marked in edges:
            if port not in (left, right):
                continue
            # A directed marked ring edge has its own coin, indexed by left.
            coin = (heads >> left) & 1 if marked else 0
            term = terms[1, 1] if marked else terms["blank"]
            pair = tuple(1 + sign * (p == port) for p in (left, right))
            start = 9 * coin + LAW.SRC.index(pair)
            for end in np.flatnonzero(term[:, start]):
                final = LAW.SRC[end % 9]
                other = (left, right)[final.index(1 + sign)]
                mask = heads ^ ((coin ^ (end // 9)) << left) if marked else heads
                column[other, mask] += int(term[end, start])
            filled = 9 * coin + LAW.SRC.index((1, 1))
            column[port, heads] -= int(term[filled, filled])
        columns.append(column)
        for word, weight in column.items():
            if weight and word not in index:
                index[word] = len(words)
                words.append(word)
    h8 = np.zeros((len(words), len(words)), dtype=int)
    for j, column in enumerate(columns):
        for word, weight in column.items():
            if weight:
                h8[index[word], j] = weight
    np.testing.assert_array_equal(h8, h8.T)
    return h8


def positive_roots(h8, bound):
    """Certify the inertia above zero and exclude the interval (0,bound]."""
    poly = Poly(Matrix(h8.tolist()).charpoly().as_expr())
    x = poly.gens[0]
    while poly.eval(0) == 0:
        poly = poly.exquo(Poly(x, x))

    # All eigenvalues are real. Descartes' rule is exact for a real-rooted
    # polynomial with nonzero constant term, applied to both +/- arguments.
    def signs(coeffs):
        seq = [v > 0 for v in coeffs if v != 0]
        return sum(int(bool(a != b)) for a, b in zip(seq, seq[1:]))

    above_zero = signs(poly.all_coeffs())
    limit = Rational(8 * bound.numerator, bound.denominator)
    shifted = poly.to_field().shift(limit)
    assert shifted.eval(0) != 0
    assert above_zero == signs(shifted.all_coeffs())
    return above_zero


def run():
    board, _, ring, _, _, compensator = winding.flight.layout()
    n = len(ring)
    assert n == 15
    edges = [(i, (i + 1) % n, True) for i in range(n)]
    reference = 8 * sum(winding.sector(n)[1])
    for sign in (-1, 1):
        h = sector(sign, edges)
        np.testing.assert_array_equal(
            h, reference if sign < 0 else reference[::-1, ::-1]
        )
        assert np.linalg.matrix_rank(h) == n
        assert positive_roots(h, Fraction(1, 16)) == 0

    minima = {}
    for mode in ("replace", "side", "ordinary_side"):
        positive = {-1: [], 1: []}
        for i, port in enumerate(ring):
            configuration = edges.copy()
            if mode == "replace":
                a, b, _ = configuration[i]
                configuration[i] = (a, b, False)
            else:
                # Select a native, distinct degree-five branch. Other seams
                # incident to it are sealed and its private contact is declared.
                outer = min(
                    q
                    for q in board.neighbours[port]
                    if q not in ring and q != compensator and board.degree[q] == 5
                )
                assert ring[i] in board.neighbours[outer]
                configuration.append((i, n, False))
            for sign in positive:
                h = sector(sign, configuration, n if mode == "ordinary_side" else None)
                lower = {
                    "replace": Fraction(5),
                    "side": Fraction(7, 2),
                    "ordinary_side": Fraction(9, 20) if sign < 0 else Fraction(7, 2),
                }[mode]
                count = positive_roots(h, lower)
                levels = np.linalg.eigvalsh(h) / 8
                assert count == np.count_nonzero(levels > float(lower))
                positive[sign].append(min(levels[levels > float(lower)]))
        minima[mode] = {
            sign: round(min(levels), 9) for sign, levels in positive.items()
        }
    return minima


if __name__ == "__main__":
    print(run())
