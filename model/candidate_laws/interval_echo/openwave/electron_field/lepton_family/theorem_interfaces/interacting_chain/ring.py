"""The fixed repulsive contact chain in the archived onsite spin frame."""

from functools import cache
from hashlib import sha256
from itertools import combinations
import json

import numpy as np
from scipy import sparse

from ..charged_endpoints.measure import read_archive


ANCESTRY = {
    "model/mass_generation_audit/chain.py": "49c1d92d2dde9993069eb0779d30188f71a98691a77f32d13b2e00e67251125a",
    "model/mass_generation_audit/results.json": "712f46a8e2bac224e93577313aa9dc452c2df4d4f5a5d69d1b6bd88da84c397e",
    "model/periodic_momentum/momentum.py": "f94e72cd4a49dc6d621adfa6ac6f4cdd2afa1999183ea766a2fbd2688db019b9",
}


def source():
    for name, digest in ANCESTRY.items():
        assert sha256(read_archive(name)).hexdigest() == digest
    data = json.loads(read_archive("model/mass_generation_audit/results.json"))[
        "unchanged_source_join"
    ]
    assert data["source"] == "periodic_momentum.PeriodicGeometry((L,1,1),c=1/8,U=1)"
    assert (
        data["complete_Fock_witness_dimension"] == 256 and not data["contact_changed"]
    )
    assert data["hopping_after_exact_frame"] == "t=c/2=1/16"
    sigma_x = np.array([[0, 1], [1, 0]])
    for x in range(4):
        g = np.linalg.matrix_power(-1j * sigma_x, x)
        next_g = np.linalg.matrix_power(-1j * sigma_x, x + 1)
        np.testing.assert_array_equal(g @ (1j * sigma_x) @ next_g.conj().T, -np.eye(2))
        assert np.linalg.det(g) == 1


def car(word, mode):
    return word ^ (1 << mode), (-1) ** ((word & ((1 << mode) - 1)).bit_count())


@cache
def sector(length, up, down):
    def patterns(count, spin):
        return [
            sum(1 << (2 * x + spin) for x in sites)
            for sites in combinations(range(length), count)
        ]

    words = sorted(u + d for u in patterns(up, 0) for d in patterns(down, 1))
    return words, {word: j for j, word in enumerate(words)}


def operator(rows, columns, values, shape):
    return sparse.csr_matrix((values, (rows, columns)), shape=shape)


@cache
def hamiltonian(length, up, down):
    """Build 16H, with integer entries, in one complete spin-count block."""
    words, index = sector(length, up, down)
    rows, cols, values = [], [], []
    for col, word in enumerate(words):
        doubles = sum((word >> (2 * x)) & 3 == 3 for x in range(length))
        if doubles:
            rows.append(col)
            cols.append(col)
            values.append(16 * doubles)
        for source_mode in range(2 * length):
            if not word & (1 << source_mode):
                continue
            intermediate, first = car(word, source_mode)
            for step in (-2, 2):
                target_mode = (source_mode + step) % (2 * length)
                if intermediate & (1 << target_mode):
                    continue
                target, second = car(intermediate, target_mode)
                rows.append(index[target])
                cols.append(col)
                values.append(-first * second)
    h = operator(rows, cols, values, (len(words), len(words)))
    assert (h - h.T).nnz == 0
    return h


def shifted(word, length):
    wrap = word >> (2 * (length - 1))
    shifted_word = ((word << 2) % (1 << (2 * length))) + wrap
    return shifted_word, (-1) ** (
        wrap.bit_count() * (word.bit_count() - wrap.bit_count())
    )


@cache
def orbits(length, up, down, k):
    words, index = sector(length, up, down)
    visited = set()
    columns = []
    for seed in words:
        if seed in visited:
            continue
        orbit, j, sign = [], seed, 1
        while j not in visited:
            visited.add(j)
            orbit.append((index[j], sign))
            j, perm_sign = shifted(j, length)
            sign *= perm_sign
        assert j == seed
        if (k * len(orbit)) % length == (0 if sign == 1 else length // 2):
            columns.append(orbit)
    return columns


def fourier_basis(length, up, down, k):
    blocks = orbits(length, up, down, k)
    rows, cols, data = [], [], []
    for col, block in enumerate(blocks):
        for n, (row, sign) in enumerate(block):
            rows.append(row)
            cols.append(col)
            data.append(
                sign
                * np.exp(-2j * np.pi * (n * k % length) / length)
                / np.sqrt(len(block))
            )
    return operator(rows, cols, data, (len(sector(length, up, down)[0]), len(blocks)))


def field(length, up, down, site, spin, add):
    words = sector(length, up, down)[0]
    step = 1 if add else -1
    target_words, index = sector(
        length, up + step * (spin == 0), down + step * (spin == 1)
    )
    rows, cols, values = [], [], []
    mode = 2 * site + spin
    for col, word in enumerate(words):
        if bool(word & (1 << mode)) == add:
            continue
        changed, phase = car(word, mode)
        rows.append(index[changed])
        cols.append(col)
        values.append(phase)
    return operator(rows, cols, values, (len(target_words), len(words)))


def spin_plus(length, up, down):
    """The calibrated global spin raising operator in the onsite frame."""
    return sum(
        field(length, up, down - 1, x, 0, True) @ field(length, up, down, x, 1, False)
        for x in range(length)
    )
