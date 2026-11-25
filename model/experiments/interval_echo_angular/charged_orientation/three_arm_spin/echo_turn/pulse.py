"""Physical two-endpoint unwrite and rewrite on every native head/load column."""

from functools import cache
from math import pi

import numpy as np
from scipy.linalg import eigh

from ...live_trimer.measure import LAW


@cache
def operators():
    exchange, contrast = LAW.local_term_parts()
    private = [LAW.SRC.index((n, 0)) for n in range(3)]
    np.testing.assert_allclose(contrast.diagonal()[private], [0, 1, 4], atol=2e-14)
    h = exchange + 8 * contrast
    ev, cols = eigh(h)
    np.testing.assert_allclose(ev, np.rint(ev), atol=3e-12)
    echo = LAW.echo()
    np.testing.assert_allclose(echo @ echo, np.eye(99), atol=3e-14)
    return ev, cols, echo


@cache
def wait(t):
    ev, cols, _ = operators()
    return (cols * np.exp(-1j * ev * t)) @ cols.T


def col(blank, end, hole):
    """The blank or a fixed even sign parity direction, and the vacancy address."""
    p = LAW.SRC.index((0, 1) if hole == 0 else (1, 0))
    z = np.zeros(99)
    if blank:
        z[p] = 1.0
    else:
        for sign in (-1, 1):
            z[9 + 9 * LAW.LABELS.index((end, sign)) + p] = 2**-0.5
    return z


def phase(end, t):
    # A physical blank private zero reference provides contrast cost 8*n².
    diag = np.array([p[end] ** 2 for p in LAW.SRC] * 11)
    return np.exp(-8j * t * diag)


@cache
def primitive(a):
    """Orient the written seam from physical endpoint a to the opposite g."""
    g = 1 - a
    echo = operators()[2]
    # Full directed pulse: the written head is fixed at g on both outputs.
    # Probe position distinguishes the incoming bit; its quantum state stays material.
    carry = wait(pi / 8)
    effects = []
    for bit in (0, 1):
        input_col = col(False, a if bit == 0 else g, a)
        output_col = col(False, g, g if bit == 0 else a)
        z = np.vdot(output_col, carry @ input_col)
        np.testing.assert_allclose(carry @ input_col, z * output_col, atol=5e-14)
        np.testing.assert_allclose(z, 1 if bit == 0 else -1, atol=5e-14)
        effects.append(z)
    # Make the two written branches of W in phase. Echo takes them to the
    # symmetric blank vacancy. Reverse the blank splitter to restore the stock.
    before = phase(a, 3 * pi / 16)[:, None] * wait(pi / 16)
    after = wait(3 * pi / 2) * phase(a, pi / 16)[None, :]
    reset = after @ echo @ before
    fixed = col(False, g, g)
    blank = col(True, 0, g)
    amp = np.vdot(blank, reset @ fixed)
    np.testing.assert_allclose(reset @ fixed, amp * blank, atol=2e-13)
    np.testing.assert_allclose(
        reset.conj().T @ blank, amp.conjugate() * fixed, atol=2e-13
    )
    assert abs(abs(amp) - 1) < 1e-14
    null = after @ before @ fixed
    assert np.linalg.norm(null[:9]) < 1e-13
    # For equal loads the returned stocked word is exact on real physical
    # (1,1) closures and real private (0,0) closures.
    for pair in ((0, 0), (1, 1)):
        k = LAW.SRC.index(pair)
        assert np.linalg.norm(echo[:9, k]) < 1e-14
        assert np.linalg.norm(echo @ echo[:, k] - np.eye(99)[:, k]) < 1e-14
    full = wait(pi)[:9, :9]
    swap = np.eye(9)[:, [LAW.SRC.index((y, x)) for x, y in LAW.SRC]]
    np.testing.assert_allclose(full, swap, atol=2e-13)
    return tuple(effects), amp
