"""One charge, twelve physical ports, and a shared exposure on a three-address return."""

from fractions import Fraction
from functools import cache
from math import asin, pi, sqrt

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.executable_model.symmetry import rotations

from .. import measure as marked


def chart():
    a, b, plus, minus = marked.addresses(marked.CELL_SEAMS[0])
    start = [a, plus, minus, b]
    for g in rotations():
        if any(g[g[g[i]]] != i or g[i] == i for i in range(12)):
            continue
        rows = [start, [g[i] for i in start], [g[g[i]] for i in start]]
        if len(set(sum(rows, []))) == 12:
            return np.array(rows), list(g)
    raise ArithmeticError("Three rotated copies do not tile this cell.")


def preparation(ports, weights):
    remaining = 1.0
    ops = []
    for b, w in zip(ports[1:], weights[1:], strict=True):
        if w:
            ops.append(
                ("rotate", int(ports[0]), int(b), asin(min(1.0, sqrt(w / remaining))))
            )
            remaining -= w
    return ops


def invert(word):
    return [(kind, a, b, -phi) for kind, a, b, phi in reversed(word)]


def compile_word(word):
    U = np.eye(12, dtype=complex)
    for kind, a, b, phi in word:
        if kind == "rotate":
            marked.rotate(U, a, b, phi, 1)
        elif kind == "phase":
            marked.phase(U, a, phi, 1)
        else:
            trail = marked.route(a, b)
            forward = list(zip(trail, trail[1:]))
            for left, right in forward + forward[:-1][::-1]:
                marked.contact(U, left, right, pi, 1)
    return U


@cache
def build():
    sites, rotation = chart()
    a, b, plus, minus = marked.addresses(marked.CELL_SEAMS[0])
    uniform = abs(marked.prepare(1)) ** 2
    central = float(uniform[a] + uniform[b])
    tip = float(uniform[plus])
    other = float(uniform[minus])
    theta = tip * (1 - central) + (1 - tip) * central
    np.testing.assert_allclose(
        [central, tip, other, theta], [1 / 6, 1 / 12, 1 / 12, 2 / 9], atol=1e-13
    )
    # The receiver and comparison pulses only use native vacancy contacts and
    # private dwells. A one-radian-per-probability convention supplies theta.
    input_word = []
    output_word = []
    select = []
    for row in sites:
        input_word.extend(
            preparation(row, [central, tip, other, 1 - central - tip - other])
        )
        output_word.extend(preparation(row, [1 / 3, 1 / 3, 1 / 3, 0]))
        select.extend(
            [("phase", int(row[1]), 0, theta), ("phase", int(row[2]), 0, -theta)]
        )
    for i, j in [(1, 1), (2, -1)]:
        cycle = list(map(int, sites[:, i]))
        order = (
            [(cycle[1], cycle[2]), (cycle[0], cycle[1])]
            if j == 1
            else [(cycle[0], cycle[1]), (cycle[1], cycle[2])]
        )
        select.extend(("swap", a, b, 0.0) for a, b in order)
    word = input_word + select + invert(output_word)
    U = compile_word(word)
    back = compile_word(invert(word))
    np.testing.assert_allclose(back @ U, np.eye(12), atol=2e-12)
    ports = sites[:, 0]
    A = U[np.ix_(ports, ports)]
    shift = np.roll(np.eye(3), 1, axis=0)
    target = (
        sqrt(2) * np.eye(3) + np.exp(1j * theta) * shift + np.exp(-1j * theta) * shift.T
    ) / 6
    np.testing.assert_allclose(A, target, atol=1e-12)
    F = np.exp(2j * np.pi * np.outer(np.arange(3), np.arange(3)) / 3) / sqrt(3)
    p = abs(np.diag(F.conj().T @ A @ F)) ** 2
    # The distinct historical 3,600-address analyzer has coherent flag
    # cardinalities 500,250,250 and common denominator 60 sqrt(3). Its three
    # coefficients have precisely the same shape, up to reversal of the
    # cyclic convention. This equality compares effects, not time evolutions.
    historical_scale = Fraction(250, 60**2 * 3)
    assert Fraction(1, 36) / historical_scale == Fraction(6, 5)
    old = (
        np.sqrt(500) + 2 * np.sqrt(250) * np.cos(theta + 2 * np.pi * np.arange(3) / 3)
    ) ** 2 / (60**2 * 3)
    np.testing.assert_allclose(np.sort(p), np.sort(old) * 6 / 5, atol=2e-13)
    return sites, rotation, theta, U, back, F, p, word


def tick(phase):
    sites, _, _, U, back, F, p, _ = build()
    ports = sites[:, 0]
    pulse = compile_word([("phase", int(a), 0, phase) for a in ports])
    raw = (back @ pulse @ U)[np.ix_(ports, ports)]
    measured = F.conj().T @ raw @ F
    np.testing.assert_allclose(
        measured, np.diag(1 + np.expm1(1j * phase) * p), atol=3e-12
    )
    # A finite inverse pulse has roundoff. Reuse the exactly proved dilation
    # identity for the long product after checking its direct 12-rail columns.
    return 1 + np.expm1(1j * phase) * p


def release():
    sites, _, _, _, _, F, _, _ = build()
    board = federation(0)
    h = board.laplacian() / 2 + diags(marked.KAPPA * (2 * board.degree - 6))
    basis = np.eye(board.n)[:, sites[:, 0]]
    rest = h @ basis
    leak = rest - basis @ (basis.T @ rest)
    gram = leak.T @ leak
    np.testing.assert_array_equal(4 * gram, 3 * np.eye(3) + 2 * np.ones((3, 3)))
    check = expm_multiply(-1j * h, basis @ F)
    return dict(
        grams_times_four=np.rint(4 * gram).astype(int).tolist(),
        minimum_leakage_squared=float(np.linalg.eigvalsh(gram)[0]),
        survival_on_three_ports_at_one=np.sum(
            abs(check[sites[:, 0]]) ** 2, axis=0
        ).tolist(),
    )
