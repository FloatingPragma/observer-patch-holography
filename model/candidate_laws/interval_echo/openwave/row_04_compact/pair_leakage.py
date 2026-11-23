"""Short-time loss of the M10.3 home-and-adjacent pair record.

P projects onto the thirty internal seams of one home in the distinct-port
two-load sector. For a normalized initial state x in P, the probability of
quiet divided by t² tends to ||(I-P) H x||². The agreement interaction is
diagonal. Its contribution to this derivative vanishes for every real κ.
"""
from fractions import Fraction

import numpy as np
import sympy as sp

from model.executable_model import cell, federation
from model.executable_model.carrier import HERE as CARRIER_HOME
from ..submission import row_04_check


def exit_block(board, home_edges):
    """Nonzero rows of (I-P) 2H P, indexed by the actual pair after a swap."""
    patch = set(home_edges)
    exits = {}
    for column, (a, b) in enumerate(home_edges):
        for moving, fixed in ((a, b), (b, a)):
            for target in board.neighbours[moving]:
                pair = tuple(sorted((fixed, target)))
                if target == fixed or pair in patch:
                    continue
                exits.setdefault(pair, [0] * len(home_edges))[column] -= 1
    return exits


def check_internal_spectrum(block):
    """Certify the least squared singular value on the thirty internal exits."""
    if block.shape != (30, 30):
        raise ValueError("Expected thirty internal exit rows and thirty home seams")
    x = sp.Symbol("x")
    expected = (x-16)*(x-4)**15*(x-1)**8*(x**2-12*x+16)**3
    if (block.T*block).charpoly(x).all_coeffs() != sp.Poly(expected, x).all_coeffs():
        raise ArithmeticError("Internal exit Gram spectrum differs from the cell certificate")
    # The quadratic roots are 6 ± 2√5. Both exceed one since 5² > 4·5.
    return Fraction(1, 4)


def certify_pair_leakage():
    local = cell()
    inside = sorted(exit_block(local, local.seams).items())
    floor = check_internal_spectrum(sp.Matrix([row for _, row in inside]))
    # An exterior seam acting on an internal pair puts one load outside.
    # Thus all internal exit rows embed unchanged; extra rows can only
    # increase the squared norm. Check the identical labelling in each home.
    counts = []
    for level, size in enumerate((20, 80, 320)):
        board = federation(level)
        cells = [[] for _ in range(board.cells)]
        for a, b in board.seams:
            if a//12 == b//12:
                cells[a//12].append((a % 12, b % 12))
        if board.cells != size or any(edges != local.seams for edges in cells):
            raise ValueError("The classified internal pair block does not embed in every home")
        counts.append(size)

    # Independent assembly from the hash-checked M10.3 carrier and its
    # published two-load generator. Compare every entry, including zeros,
    # in all thirty detector columns of every home on the scored board.
    board = federation(0)
    edges, index, exchange, adjacent = row_04_check.screen(
        CARRIER_HOME / "carriers" / "FEDERATION_LEVEL0.npz"
    )
    if not np.array_equal(edges, board.seams):
        raise ValueError("M10.3 and the certified board have different seams")
    sharp = (156, 157)
    for home in range(board.cells):
        patch = [(12*home+a, 12*home+b) for a, b in local.seams]
        cols = [index[a, b] for a, b in patch]
        if not np.all(adjacent[cols] == 1):
            raise ArithmeticError("M10.3 fails to recognize a home seam")
        exits = exit_block(board, patch)
        pairs = list(exits)
        rows = [index[pair] for pair in pairs]
        expected = np.array([exits[pair] for pair in pairs])
        actual = (2*exchange[:, cols]).tocsr()
        if (not np.array_equal(actual[rows].toarray(), expected)
                or actual.count_nonzero() - actual[cols].count_nonzero() != np.count_nonzero(expected)):
            raise ArithmeticError("M10.3 quiet rows disagree with the endpoint swaps")
        if home == 13:
            if sharp != next(edge for edge in patch if any(board.degree[p] == 6 for p in edge)):
                raise ValueError("Changed M10.3 sharp launch")
            column = expected[:, patch.index(sharp)]
            sharp_loss = Fraction(int(column @ column), 4)
    return tuple(counts), floor, sharp_loss


if __name__ == "__main__":
    counts, floor, sharp = certify_pair_leakage()
    print(f"{sum(counts)} homes: short-time quiet probability divided by t² tends to at least {floor}")
    print(f"M10.3 launch (156, 157): exact limit {sharp}")
