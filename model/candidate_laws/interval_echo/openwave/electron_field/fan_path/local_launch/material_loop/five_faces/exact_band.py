"""Exact wheel reduction and rational stationary current enclosure."""

from dataclasses import dataclass
from fractions import Fraction
from itertools import pairwise, product
from math import isqrt

from model.executable_model import federation

CENTER = 1047
RIM = (1046, 1048, 1053, 1052, 1050)
N = len(RIM)


@dataclass(frozen=True)
class Bounds:
    lo: Fraction
    hi: Fraction

    def __add__(self, other):
        return Bounds(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other):
        return Bounds(self.lo - other.hi, self.hi - other.lo)

    def __mul__(self, other):
        corners = [a * b for a, b in product((self.lo, self.hi), (other.lo, other.hi))]
        return Bounds(min(corners), max(corners))

    def __truediv__(self, other):
        assert other.lo * other.hi > 0
        return self * Bounds(1 / other.hi, 1 / other.lo)

    def sign(self):
        assert self.lo * self.hi > 0, self
        return 1 if self.lo > 0 else -1


def point(value):
    q = Fraction(value)
    return Bounds(q, q)


def square_root(x):
    scale = 10**26
    assert x.lo > 0
    a, b = (isqrt(q.numerator * scale**2 // q.denominator) for q in (x.lo, x.hi))
    out = Bounds(Fraction(a, scale), Fraction(b + 1, scale))
    assert out.lo**2 <= x.lo <= x.hi <= out.hi**2
    return out


def native_transitions():
    board = federation(2)
    faces = {frozenset(face) for face in board.faces()}
    assert all(frozenset((CENTER, a, b)) in faces for a, b in pairwise((*RIM, RIM[0])))
    assert all(board.degree[p] == 5 for p in RIM)
    # A spoke bit 1 points to the rim. A rim bit 1 points to its clockwise port.
    basis = {}
    for y, bits, x in product(range(N), product((0, 1), repeat=N), range(N)):
        if any(bits[j - 1] + 1 - bits[j] + (j == y) != 1 + (j == x) for j in range(N)):
            continue
        ell = sum(bits)
        assert x == (y + ell) % N
        assert bits == tuple(int((j - y) % N < ell) for j in range(N))
        assert (x, ell) not in basis
        spokes = tuple(int(j == y) for j in range(N))
        basis[x, ell] = (spokes, bits)
    assert set(basis) == set(product(range(N), range(N + 1)))
    inverse = {(x, s, b): (x, ell) for (x, ell), (s, b) in basis.items()}
    for (x, ell), (s, b) in basis.items():
        ring, material = set(), set()
        for j in range(N):
            k = (j + 1) % N
            updated = tuple(bit ^ (i == j) for i, bit in enumerate(b))
            # A circulating triangle reverses all three heads.
            if (s[j], b[j], s[k]) in ((1, 1, 0), (0, 0, 1)):
                next_spokes = tuple(bit ^ (i in (j, k)) for i, bit in enumerate(s))
                ring.add(inverse[x, next_spokes, updated])
            # The positive load leaves a contact that points to its port.
            if x == (k if b[j] else j):
                destination = j if b[j] else k
                material.add(inverse[destination, s, updated])
        steps = [d for d in (-1, 1) if 0 <= ell + d <= N]
        assert ring == {(x, ell + d) for d in steps}
        assert material == {((x + d) % N, ell + d) for d in steps}


DIAGONALS = (Fraction(25, 2), 17, 17, 17, 17, Fraction(25, 2))


def continuant(energy, beta):
    """Enclose det(EI-J), its E and beta derivatives, and all leading minors."""
    older, prior = (point(0),) * 3, (point(1), point(0), point(0))
    minors = [prior[0]]
    for diagonal in DIAGONALS:
        x = energy - point(diagonal)
        value, by_e, by_beta = prior
        ov, oe, ob = older
        row = (
            x * value - beta * ov,
            value + x * by_e - beta * oe,
            x * by_beta - ov - beta * ob,
        )
        older, prior = prior, row
        minors.append(row[0])
    return prior, minors


def negative_pivots(energy, beta):
    # Consecutive leading-minor ratios are the exact no-pivot LDL pivots.
    signs = [minor.sign() for minor in continuant(point(energy), beta)[1]]
    return sum(a != b for a, b in pairwise(signs))


def certify():
    native_transitions()
    # For the k=+1 Fourier convention exp(-2pi ix/5), the upper-ell bond is
    # -(g + t exp(i(2pi+phi)/5)). Its phase can be gauged away on the open leg.
    # The degree offsets are 8 + (g+t) at endpoints, 8 + 2(g+t) inside.
    # beta is the bond modulus squared, with cos(2pi/5) = (sqrt(5)-1)/4.
    sqrt5 = square_root(point(5))
    beta = point("15.25") + sqrt5
    sin_phase = square_root((point(5) + sqrt5) / point(8))
    lo, hi = Fraction("8.52"), Fraction("8.54")
    assert negative_pivots(lo, beta) == 6 and negative_pivots(hi, beta) == 5
    for _ in range(50):
        center = (lo + hi) / 2
        if negative_pivots(center, beta) == 6:
            lo = center
        else:
            hi = center
    energy = Bounds(lo, hi)
    (_, by_e, by_beta), _ = continuant(energy, beta)
    # beta'(0) = -(4/5) sin(2pi/5), and simple-root differentiation gives
    # E'(0) = -P_beta beta'(0)/P_E. Rotation gives ∂H/∂phi = -mean bond J.
    derivative = point("0.8") * sin_phase * by_beta / by_e
    assert (
        Fraction("8.5260398792981")
        < energy.lo
        < energy.hi
        < Fraction("8.5260398792982")
    )
    assert (
        Fraction("0.1488899031607")
        < derivative.lo
        < derivative.hi
        < Fraction("0.1488899031608")
    )
    return energy, derivative


if __name__ == "__main__":
    energy, derivative = certify()
    print("Certified E in (8.5260398792981, 8.5260398792982).")
    print("Certified dE/dphi in (0.1488899031607, 0.1488899031608) for k=+1.")
