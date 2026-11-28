"""Sharp native transition tables for a binary rail and its in-place collar."""

from functools import lru_cache
import numpy as np
from model.candidate_laws.interval_echo.candidate import native_law as law

MEM = [None, *[(c, d) for c in (0, 1) for d in (-1, 0, 1)]]
SIGNED = np.array([0, *[d for c, d in MEM[1:]]], dtype=np.int16)
FLAG = np.abs(SIGNED)
CURRENT = np.array(
    [0, *[((d + 1) // 2 - c) if d else 0 for c, d in MEM[1:]]], dtype=np.int16
)
LEFT = np.tile([0, 0, 1, 1], len(MEM))
RIGHT = np.tile([0, 1, 0, 1], len(MEM))
MEMORY = np.arange(28) // 4
LEFT_CHARGE = LEFT + CURRENT[MEMORY]
RIGHT_CHARGE = RIGHT - CURRENT[MEMORY]
SETTINGS = ("retained", "retained_exchange", "fresh")


def exact_kernels():
    # F is sqrt(2) times the native write on four load words. The other
    # 24 labels span six written memories times the same four words.
    f = np.zeros((24, 4), dtype=np.int64)
    r = np.zeros((24, 24), dtype=np.int64)
    for a in (0, 1):
        for b in (0, 1):
            i = 2 * a + b
            for c in (0, 1):
                y = law.native_event(a, b, c)
                j = 4 * (MEM.index((c, a - b)) - 1) + 2 * y[0] + y[1]
                f[j, i] = 1
    for j, m in enumerate(MEM[1:]):
        c, d = m
        k = MEM.index((c, -d)) - 1
        r[4 * k : 4 * k + 4, 4 * j : 4 * j + 4] = np.eye(4, dtype=np.int64)
    f2 = f @ f.T
    eye = np.eye(24, dtype=np.int64)
    write = np.zeros((28, 28), dtype=np.int64)
    write[4:, :4] = 32 * f**2
    write[:4, 4:] = 32 * f.T**2
    write[4:, 4:] = 16 * (2 * eye - f2) ** 2
    # At a native half contact only the reversal parity carries phase.
    # The echo-conjugated reversal on blank loads is SWAP. Its written
    # block is R + WW^T(I-R); there are no cross terms.
    swap = np.eye(4, dtype=np.int64)[[0, 2, 1, 3]]
    np.testing.assert_array_equal(r @ f, f @ swap)
    twice_s = np.zeros((28, 28), dtype=np.int64)
    twice_s[:4, :4] = 2 * swap
    twice_s[4:, 4:] = 2 * r + f2 @ (eye - r)
    np.testing.assert_array_equal(twice_s @ twice_s, 4 * np.eye(28, dtype=np.int64))
    half = 32 * np.eye(28, dtype=np.int64) + 8 * twice_s**2
    return write, half


@lru_cache(maxsize=1)
def operators():
    basis = [
        (0 if m is None else 1 + law.LABELS.index(m)) * 9 + law.SRC.index((a, b))
        for m in MEM
        for a in (0, 1)
        for b in (0, 1)
    ]
    U = law.echo()
    h0, h1 = law.local_term_parts()
    e, v = np.linalg.eigh(h0 + 8 * h1)
    V = (v * np.exp(-1j * np.pi / 2 * e)) @ v.T
    # Half contact on the written reversal-even/odd sectors at kappa=8.
    # The integer-contrast phase is one. U maps this spectral representation
    # to the native external port and collar basis.
    rr, _ = law._label_ops()
    R = np.zeros_like(U)
    R[:9, :9] = np.eye(9)
    R[9:, 9:] = np.kron(rr, np.eye(9))
    with_reversal = U @ (((np.eye(99) + R) - 1j * (np.eye(99) - R)) / 2) @ U
    np.testing.assert_allclose(with_reversal, V, rtol=0, atol=4e-13)
    outside = np.setdiff1d(np.arange(99), basis)
    ps = []
    for op in (U, V):
        assert np.max(abs(op[np.ix_(outside, basis)])) < 4e-13
        p = abs(op[np.ix_(basis, basis)].T) ** 2
        ints = np.rint(p * 64).astype(int)
        np.testing.assert_allclose(ints / 64, p, rtol=0, atol=1e-12)
        np.testing.assert_array_equal(ints.sum(axis=1), 64)
        np.testing.assert_array_equal(ints, ints.T)
        ps.append(ints)
    write, exchange = exact_kernels()
    for expected, actual in zip((write, exchange), ps, strict=True):
        np.testing.assert_array_equal(expected, actual)
    x, y = np.nonzero(write)
    for label in (LEFT_CHARGE, RIGHT_CHARGE):
        np.testing.assert_array_equal(label[x], label[y])
    x, y = np.nonzero(exchange)
    assert np.any(LEFT_CHARGE[x] != LEFT_CHARGE[y])
    np.testing.assert_array_equal(FLAG[MEMORY[x]], FLAG[MEMORY[y]])
    matrices = {
        "retained": write * 64,
        "retained_exchange": write @ exchange,
        "fresh": write[np.arange(28) % 4] * 64,
    }
    return matrices, write, exchange


@lru_cache(maxsize=1)
def samplers():
    result = {}
    for name, mat in operators()[0].items():
        g = int(np.gcd.reduce(mat.ravel()))
        table = mat // g
        rows = np.array(
            [np.repeat(np.arange(28, dtype=np.uint8), row) for row in table]
        )
        assert rows.shape[1] & (rows.shape[1] - 1) == 0
        result[name] = rows
    return result
