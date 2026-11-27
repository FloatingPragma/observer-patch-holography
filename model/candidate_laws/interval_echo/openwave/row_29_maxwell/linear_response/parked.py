"""Exact fixed-load reduction of the 99-state native seam and a local echo read."""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law

HERE = Path(__file__).resolve().parent


def integer_blocks():
    """Return four times both coefficients, exact on blank and written memories.

    The write is Wp / sqrt(2). Since the pre-echo read preserves the image of W,
    the written-to-written compression of U G U* equals G (I - W W*).
    """
    w = law._write_pattern()
    rev, contrast = law._label_ops()
    swap, square = law._load_ops()
    i9 = np.eye(9, dtype=int)
    r = np.kron(rev, i9)
    c = np.kron(contrast, i9)
    np.testing.assert_array_equal(r @ w, w @ swap)
    np.testing.assert_array_equal(c @ w, w @ square)
    np.testing.assert_array_equal(w.T @ w, 2 * i9)
    i90 = np.eye(90, dtype=int)
    complement2 = 2 * i90 - w @ w.T
    exch4 = np.zeros((99, 99), dtype=int)
    cost4 = np.zeros_like(exch4)
    exch4[:9, :9], cost4[:9, :9] = 2 * (i9 - swap), 4 * square
    exch4[9:, 9:], cost4[9:, 9:] = (i90 - r) @ complement2, 2 * c @ complement2
    for exact, native in zip((exch4, cost4), law.local_term_parts(), strict=True):
        np.testing.assert_allclose(exact / 4, native, rtol=0, atol=2e-14)
    return exch4, cost4


def certificate():
    exchange, cost = integer_blocks()
    rows = []
    for load in (0, 1, 2):
        ix = np.arange(11) * 9 + law.SRC.index((load, load))
        rest = np.setdiff1d(np.arange(99), ix)
        for coefficient in (exchange, cost):
            assert not coefficient[np.ix_(rest, ix)].any()
        for k in (0, 8):
            block4 = (exchange + k * cost)[np.ix_(ix, ix)]
            # Both terms come from the native coefficients. The confined
            # memory's characteristic energy levels are integers on both laws.
            ev = np.linalg.eigvalsh(block4 / 4)
            levels, mult = np.unique(ev.round().astype(int), return_counts=True)
            np.testing.assert_allclose(ev, np.repeat(levels, mult), atol=3e-13)
            annihilator = np.eye(11, dtype=object)
            for value in levels:
                annihilator = annihilator @ (
                    block4.astype(object) - 4 * value * np.eye(11, dtype=int)
                )
            assert not annihilator.any()
            rows.append(
                dict(
                    load=load,
                    contrast=k,
                    levels=levels.tolist(),
                    multiplicities=mult.tolist(),
                )
            )
    # Zero-reference load stays fixed and its initial blank is invariant.
    for load in (0, 1, 2):
        ix = law.SRC.index((load, 0))
        target = np.zeros(99, dtype=int)
        target[ix] = 4 * load**2
        np.testing.assert_array_equal(cost[:, ix], target)
    return rows


def witness():
    """Write local sign at total load two, let the material fly, undo the echo.

    A remote owner holds its own independent unbiased sign at unit/unit.
    Local material and the unbiased sign can also be read on entangled copies.
    """
    shots = 2**32
    rng = np.random.default_rng(290773)
    h0, h1 = integer_blocks()
    U = law.echo()
    indices = [law.SRC.index(pair) for pair in ((2, 0), (0, 2))]
    code = U[:, indices]
    material = np.arange(11) * 9 + law.SRC.index((1, 1))
    np.testing.assert_allclose(np.sum(np.abs(code[material]) ** 2, axis=0), 1)
    records = []
    # Local zero-energy written columns remain eigenstates even when all
    # contacts are on. Other memories contribute only their own phases.
    for kappa in (0, 8):
        H = (h0 + kappa * h1) / 4
        np.testing.assert_allclose(H @ code, 0, atol=3e-14)
        for t in (0, 4, 16, 32, 64):
            decoded = U @ expm(-1j * H * t) @ code
            np.testing.assert_allclose(
                abs(decoded[indices]) ** 2, np.eye(2), atol=1e-12
            )
            for bit in range(2):
                source = rng.multinomial(shots, np.eye(2)[bit]).tolist()
                # Sum over the independently prepared remote sign. Both paths
                # use the same measured local address read.
                remote_probs = np.mean(abs(decoded[indices]) ** 2, axis=1)
                remote = rng.multinomial(
                    shots, remote_probs / remote_probs.sum()
                ).tolist()
                records.append(
                    dict(contrast=kappa, wait=t, bit=bit, source=source, remote=remote)
                )
    return records


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = (
        json.dumps({"blocks": certificate(), "local_echo": witness()}, indent=2) + "\n"
    )
    path = HERE / "collar.json"
    if args.write:
        path.write_text(text)
    else:
        assert text == path.read_text()
