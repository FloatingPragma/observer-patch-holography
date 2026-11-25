"""Exact signed density-transfer coefficients for the two-path probe."""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import factorial

from model.executable_model import cell
from ..measure import clock_seam
from ..shell_response.measure import PARITY_PORTS as PARITY_WORD


def paths():
    board = cell()
    a, b = clock_seam()
    s, r = sorted(set(board.neighbours[a]) & set(board.neighbours[b]))
    assert s not in board.neighbours[r]
    return a, b, s, r


def parity_terminals():
    board = cell()
    source = set(PARITY_WORD)
    outside = sorted(set(range(board.n)) - source - set(clock_seam()))
    choices = []
    for s, r in combinations(outside, 2):
        near_s, near_r = set(board.neighbours[s]), set(board.neighbours[r])
        if r not in near_s and near_s & near_r:
            score = len(near_s & source) + len(near_r & source)
            choices.append((-score, s, r))
    _, s, r = min(choices)
    return s, r


def brackets(poly, bonds):
    """Apply 2i[H,.] in the Hermitian Pauli basis i^|x&z| X^x Z^z."""
    result = defaultdict(int)
    for (u, v), entry in poly.items():
        for x, z, weight in bonds:
            if ((x & v).bit_count() + (z & u).bit_count()) % 2 == 0:
                continue
            new = (x ^ u, z ^ v)
            phase = ((x & z).bit_count() + (u & v).bit_count()
                     - (new[0] & new[1]).bit_count() + 2 * (z & u).bit_count()) % 4
            # 4H = -sum(XX + YY + (2k+1)ZZ) after dropping terms commuting with n_r(t).
            result[new] += (1 if phase == 1 else -1) * weight * entry
    return {pair: value for pair, value in result.items() if value}


def short_coefficients(kappa, half_order=3, terminals=None):
    a, b, sender, receiver = paths()
    if terminals is not None:
        sender, receiver = terminals
    if len({a, b, sender, receiver}) != 4:
        raise ValueError("source and terminals must be disjoint")
    bonds = []
    for i, j in cell().seams:
        ij = (1 << i) | (1 << j)
        bonds.extend([(ij, 0, 1), (ij, ij, 1), (0, ij, 2 * kappa + 1)])
    ab, s = (1 << a) | (1 << b), 1 << sender
    # The receiver is n_r=(I-Z_r)/2; the changed sender is -Z_s.
    # Tr(-Z_s Z_s)=-2 cancels the receiver's nonconstant coefficient.
    receiver_terms = [{(0, 1 << receiver): 1}]
    for _ in range(half_order):
        receiver_terms.append(brackets(receiver_terms[-1], bonds))
    signed = {}
    # Coefficients are traces of Pauli words against the source tangent.
    launches = {"exchange": {(ab, s): -2, (ab, ab | s): -2},
                "aligned": {(ab, s): 1, (ab, ab | s): 1, (0, ab | s): -2}}
    launches["parity"] = {} if sender in PARITY_WORD else {
        (0, s | sum(1 << p for p in PARITY_WORD)): 1}
    for name, launch in launches.items():
        terms = [launch]
        for _ in range(half_order):
            terms.append(brackets(terms[-1], bonds))
        slopes = []
        for n in range(half_order + 1):
            # 2i[H,.] is skew with respect to the real Hilbert-Schmidt product.
            dot = sum(c * terms[n].get(p, 0) for p, c in receiver_terms[n].items())
            slopes.append(Fraction((-1)**n * dot, 2**(2*n) * factorial(2*n)))
        signed[name] = slopes
    return signed


def shell_averages(kappa):
    """Equal weights on every ordered spectator pair at the indicated graph distance."""
    board = cell()
    a, b, _, _ = paths()
    shells = {key: [] for key in (1, 2, 3)}
    for s in range(board.n):
        for r in range(board.n):
            if len({a, b, s, r}) != 4:
                continue
            shared = set(board.neighbours[s]) & set(board.neighbours[r])
            distance = 1 if s in board.neighbours[r] else 2 if shared else 3
            coefficients = short_coefficients(kappa, 2, (s, r))
            shells[distance].append([coefficients[key][2] for key in ("exchange", "aligned")])
    return {distance: (len(rows), *(sum(x)/len(rows) for x in zip(*rows, strict=True)))
            for distance, rows in shells.items()}


if __name__ == "__main__":
    for kappa in (8, 0):
        print(kappa, short_coefficients(kappa))
