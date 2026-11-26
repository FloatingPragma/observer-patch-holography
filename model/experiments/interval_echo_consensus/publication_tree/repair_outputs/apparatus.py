"""Two local ternary echo instruments, a native repair, and an optional third rail."""

from functools import cache

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

law = candidate.native_law
INPUTS = tuple(law.SRC)
WORDS = ("00", "01", "10", "11")


@cache
def geometry(level):
    board = federation(level)
    a, b = min(tuple(sorted(edge)) for edge in board.external)
    neighbors = [
        sorted(i for i in board.neighbours[p] if i // 12 == p // 12) for p in (a, b)
    ]
    return {
        "cells": board.cells,
        "gateway": [int(a), int(b)],
        "owners": [int(a // 12), int(b // 12)],
        "standards_zero_two": [[int(i) for i in n[:2]] for n in neighbors],
        "vacant_arm": int(neighbors[0][2]),
    }


@cache
def detector_words():
    """Verify the project-and-unwrite instrument on each sharp standard and input."""
    u = law.echo()
    balanced = np.arange(law.DIM) % 9 == INPUTS.index((1, 1))
    ket = np.eye(law.DIM)
    words = []
    for measured in range(3):
        bits = []
        for reference in (0, 2):
            source = ket[:, INPUTS.index((measured, reference))]
            positive = u.T @ (balanced * (u @ source))
            bit = int(measured + reference == 2)
            np.testing.assert_allclose(positive, bit * source, atol=1e-13)
            bits.append(str(bit))
        words.append("".join(bits))
    return tuple(words)


@cache
def writers():
    # The comparator stores the same signed contrast and uniform coin at the same head.
    # Its classical load pair is preserved. Distinct inputs still have orthogonal outputs.
    native = law.echo()[:, :9]
    changed = np.zeros_like(native)
    for j, (a, b) in enumerate(INPUTS):
        for coin in (0, 1):
            head = 1 + law.LABELS.index((coin, a - b))
            changed[9 * head + j, j] = 1 / np.sqrt(2)
    np.testing.assert_allclose(changed.T @ changed, np.eye(9), atol=1e-13)
    return {
        "repair": native,
        "unchanged_loads": changed,
        "parked": np.eye(law.DIM)[:, :9],
    }


def level_probabilities(code):
    """Columns are preparations. The inaccessible head is traced out."""
    return (abs(code.reshape(11, 3, 3, -1)) ** 2).sum(axis=0)


def read_distribution(levels):
    table = np.zeros((4, 4))
    where = [WORDS.index(w) for w in detector_words()]
    table[np.ix_(where, where)] = levels
    return table


def possible(pair):
    s = sum(pair)
    return (s // 2, (s + 1) // 2)


def verify_marginal_channels():
    """Compare both local density matrices on all 81 (including non-Hermitian) matrix units."""
    v = writers()["repair"].reshape(11, 3, 3, 9)
    left = np.einsum("mabx,mcby->xyac", v, v.conj())
    right = np.einsum("mbax,mbcy->xyac", v, v.conj())
    expected = np.zeros_like(left)
    for x, pair in enumerate(INPUTS):
        for r in possible(pair):
            expected[x, x, r, r] += 0.5
    for actual in (left, right):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)


@cache
def free_wait(kappa, arm):
    exchange, contrast = law.local_term_parts()
    h = exchange + kappa * contrast
    v = writers()["repair"]
    np.testing.assert_allclose(h @ v, 0, atol=1e-12)
    # The independent exchange arm stays on its blank head throughout for any ABC state.
    small = h[:9, :9]
    np.testing.assert_allclose(h[9:, :9], 0, atol=1e-12)
    blank_arm = np.kron(small, np.eye(3)).reshape((3,) * 6)
    ac = blank_arm.transpose(0, 2, 1, 3, 5, 4).reshape(27, 27)
    both = np.kron(h, np.eye(3))
    if arm:
        both += np.kron(np.eye(11), ac)
    loads = np.indices((11, 3, 3, 3))[1:].sum(axis=0).ravel()
    full = np.kron(v, np.eye(3)[:, :1])
    result = {0: np.zeros((3, 3, 9)), 1: np.zeros((3, 3, 9))}
    for n in range(5):
        members = np.flatnonzero(loads == n)
        columns = [i for i, pair in enumerate(INPUTS) if sum(pair) == n]
        block = both[np.ix_(members, members)]
        outside = both[np.ix_(loads != n, members)]
        np.testing.assert_allclose(outside, 0, atol=1e-12)
        energy, vectors = np.linalg.eigh(block)
        prepared = vectors.T @ full[np.ix_(members, columns)]
        for i, t in enumerate((np.pi / 2, np.pi)):
            propagated = np.zeros((len(loads), len(columns)), dtype=complex)
            propagated[members] = vectors @ (
                np.exp(-1j * energy * t)[:, None] * prepared
            )
            p = abs(propagated.reshape(11, 3, 3, 3, -1)) ** 2
            result[i][:, :, columns] = p.sum(axis=(0, 3))
    return result


def phase_return(exposure, resolve_sign):
    exchange, contrast = law.local_term_parts()
    g = (exchange + candidate.KAPPA * contrast)[:9, :9]
    pulse = expm(-1j * g * exposure)
    x = INPUTS.index((0, 2))
    v = writers()["repair"] @ pulse[:, x]
    branches = []
    for a, b in INPUTS:
        ends = np.arange(law.DIM) % 9 == INPUTS.index((a, b))
        if resolve_sign:
            for d in (-2, 2):
                signs = np.repeat([False] + [label[1] == d for label in law.LABELS], 9)
                branches.append(v * ends * signs)
        else:
            branches.append(v * ends)
    u = law.echo()
    # Branches of all registered endpoint reads remain in the total, including failures.
    probabilities = sum(abs(pulse.conj().T @ (u @ w)[:9]) ** 2 for w in branches)
    np.testing.assert_allclose(sum(probabilities), 1, atol=5e-13)
    return probabilities.reshape(3, 3)


def clock():
    exchange, contrast = law.local_term_parts()
    h = (exchange + candidate.KAPPA * contrast)[:9, :9]
    index = INPUTS.index((0, 2))
    for t in (0, np.pi / 2, np.pi):
        yield abs(expm(-1j * t * h)[:, index].reshape(3, 3)) ** 2
