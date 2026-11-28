"""Independent integer eight-site construction in spin-separated Fock order."""

from fractions import Fraction
from functools import cache
from itertools import combinations

import numpy as np
from scipy import sparse

L = 8
MASK = (1 << L) - 1
ROOT_HALF = float.fromhex("0x1.6a09e667f3bcdp-1")
EPS = 2**-53
PHASES = np.array(
    [
        1,
        ROOT_HALF * (1 + 1j),
        1j,
        ROOT_HALF * (-1 + 1j),
        -1,
        ROOT_HALF * (-1 - 1j),
        -1j,
        ROOT_HALF * (1 - 1j),
    ],
    dtype=complex,
)
DWT = np.array([0, 1, ROOT_HALF, 0, 0.5, 0, 0, 0, ROOT_HALF / 2])
assert (Fraction(ROOT_HALF) - Fraction(EPS)) ** 2 < Fraction(1, 2)
assert (Fraction(ROOT_HALF) + Fraction(EPS)) ** 2 > Fraction(1, 2)


@cache
def modes(n):
    return tuple(sum(1 << i for i in indices) for indices in combinations(range(L), n))


@cache
def stock(up, down):
    words = tuple(a | (b << L) for a in modes(up) for b in modes(down))
    where = {word: j for j, word in enumerate(words)}
    return words, where


def fermion(word, orbital, create):
    if bool(word & (1 << orbital)) == create:
        return None
    sign = (-1) ** ((word & ((1 << orbital) - 1)).bit_count() & 1)
    return word ^ (1 << orbital), sign


def csr(transitions, shape):
    rows, cols, values = zip(*transitions, strict=True) if transitions else ((), (), ())
    return sparse.csr_matrix((values, (rows, cols)), shape=shape, dtype=np.int64)


@cache
def one_spin(n):
    words = modes(n)
    lookup = {w: j for j, w in enumerate(words)}
    items = []
    for col, w in enumerate(words):
        for i in range(L):
            out = fermion(w, i, False)
            if out is None:
                continue
            word, sign = out
            for distance in (-1, 1):
                move = fermion(word, (i + distance) % L, True)
                if move is not None:
                    items.append((lookup[move[0]], col, -sign * move[1]))
    return csr(items, (len(words), len(words)))


@cache
def hamiltonian(up, down):
    a, b = one_spin(up), one_spin(down)
    identity = lambda d: sparse.eye(d, dtype=np.int64)  # noqa: E731
    double = [16 * (x & y).bit_count() for x in modes(up) for y in modes(down)]
    h = sparse.kron(a, identity(b.shape[0]), format="csr")
    h += sparse.kron(identity(a.shape[0]), b, format="csr") + sparse.diags(
        double, dtype=np.int64
    )
    assert (h != h.T).nnz == 0
    return h


def rotate(word, sites=1):
    up, down = word & MASK, word >> L
    stay = (1 << (L - sites)) - 1

    def move(w):
        wrap = w >> (L - sites)
        shift = (w & stay) << sites | wrap
        phase = (-1) ** (wrap.bit_count() * (w.bit_count() - wrap.bit_count()))
        return shift, phase

    a, phase_a = move(up)
    b, phase_b = move(down)
    return a | (b << L), phase_a * phase_b


@cache
def translation(up, down, sites=1):
    words, where = stock(up, down)
    data = []
    for j, word in enumerate(words):
        out, sign = rotate(word, sites)
        data.append((where[out], j, sign))
    return csr(data, (len(words), len(words)))


@cache
def field(up, down, site, spin, create):
    old, _ = stock(up, down)
    change = 1 if create else -1
    new, where = stock(up + change * (spin == 0), down + change * (spin == 1))
    entries = []
    for j, word in enumerate(old):
        out = fermion(word, site + L * spin, create)
        if out is not None:
            entries.append((where[out[0]], j, out[1]))
    return csr(entries, (len(new), len(old)))


@cache
def raising(up, down):
    return sum(
        field(up, down - 1, site, 0, True) @ field(up, down, site, 1, False)
        for site in range(L)
    )


def check_symmetries(up, down):
    h = hamiltonian(up, down)
    t = translation(up, down)
    assert (t @ h != h @ t).nnz == 0
    assert (t**4 != translation(up, down, 4)).nnz == 0
    assert (t**8 != sparse.eye(h.shape[0], dtype=np.int64)).nnz == 0
    if up < L and down > 0:
        lowering_chart = raising(up, down)
        hh = hamiltonian(up + 1, down - 1)
        tt = translation(up + 1, down - 1)
        assert (lowering_chart @ h != hh @ lowering_chart).nnz == 0
        assert (lowering_chart @ t != tt @ lowering_chart).nnz == 0


@cache
def orbits(up, down):
    words, where = stock(up, down)
    seen = set()
    result = []
    for word in words:
        if word in seen:
            continue
        to = word
        sign = 1
        orbit = []
        while to not in seen:
            seen.add(to)
            orbit.append((where[to], sign))
            to, factor = rotate(to)
            sign *= factor
        assert to == word
        result.append((sign, orbit))
    return result


def block(up, down, k):
    """Form the exact integer coefficients for each power of exp(i pi/4)."""
    n = len(stock(up, down)[0])
    kept = [
        orbit
        for sign, orbit in orbits(up, down)
        if (k * len(orbit)) % L == (0 if sign == 1 else L // 2)
    ]
    lengths = [len(o) for o in kept]
    owner = np.full(n, -1, dtype=np.int64)
    pos = np.zeros(n, dtype=np.int64)
    sign = np.zeros(n, dtype=np.int64)
    basis_entries = []
    for i, orbit in enumerate(kept):
        for j, (index, factor) in enumerate(orbit):
            owner[index] = i
            pos[index] = j
            sign[index] = factor
            basis_entries.append(
                (index, i, factor * PHASES[-k * j % L] * DWT[len(orbit)])
            )
    shape = (n, len(kept))
    if basis_entries:
        ri, ci, vv = zip(*basis_entries, strict=True)
        basis = sparse.csr_matrix((vv, (ri, ci)), shape=shape)
    else:
        basis = sparse.csr_matrix(shape, dtype=complex)
    h = hamiltonian(up, down).tocoo()
    retain = (owner[h.row] >= 0) & (owner[h.col] >= 0)
    i, j, amplitude = h.row[retain], h.col[retain], h.data[retain]
    coeff = np.zeros((L, len(kept), len(kept)), dtype=np.int64)
    phase = (k * (pos[i] - pos[j])) % L
    np.add.at(coeff, (phase, owner[i], owner[j]), amplitude * sign[i] * sign[j])
    # Hermiticity can be checked already at the polynomial coefficient level.
    np.testing.assert_array_equal(coeff[(-np.arange(L)) % L], coeff.transpose(0, 2, 1))
    grams = np.outer(DWT[lengths], DWT[lengths])
    base = np.tensordot(PHASES, coeff, axes=1) * grams / 16
    abs_scale = float(np.linalg.norm(np.sum(abs(coeff), axis=0) * grams / 16))
    representation_error = 64 * EPS * max(1, abs_scale)
    return (base + base.conj().T) / 2, basis, representation_error
