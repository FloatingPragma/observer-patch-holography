"""Current as an exact coin difference on a written one-defect seam."""

from itertools import product

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate

LAW = candidate.native_law


def exact_columns():
    pattern = LAW._write_pattern()
    projection = 2 * np.eye(90, dtype=int) - pattern @ pattern.T
    reverse, contrast = LAW._label_ops()
    eye = np.eye(10, dtype=int)
    fiber = np.eye(9, dtype=int)
    # Eight times the written block Q G Q with Q = I - W W*, W = pattern / sqrt(2).
    doubled = [np.kron(eye - reverse, fiber), np.kron(2 * contrast, fiber)]
    charges = np.array([coin if abs(d) == 1 else 0 for coin, d in LAW.LABELS])
    catalog = []
    for name, operator in zip(("exchange", "contrast"), doubled, strict=True):
        numerator = projection @ operator @ projection
        assert not np.any(projection @ operator @ pattern)
        actual = LAW.local_term_parts()[0 if name == "exchange" else 1]
        np.testing.assert_allclose(numerator / 8, actual[9:, 9:], atol=5e-14)
        for sign in (-1, 1):
            transitions = 0
            for initial in ((1, 1), (1 + sign, 1), (1, 1 + sign)):
                word = LAW.SRC.index(initial)
                for head in range(10):
                    for row in np.flatnonzero(numerator[:, 9 * head + word]):
                        final = LAW.SRC[row % 9]
                        assert final in (initial, initial[::-1])
                        motion = int(final[1] == 1 + sign) - int(initial[1] == 1 + sign)
                        assert charges[row // 9] - charges[head] == -sign * motion
                        transitions += 1
            catalog.append(dict(term=name, sign=sign, checked_nonzero=transitions))
    return catalog


def triangle(sign, written, angle):
    """All hops allowed on a guarded triangle, zero or three active |d|=1 heads."""
    edges = [(0, 1), (1, 2), (2, 0)]
    words = list(product(range(2), repeat=3)) if written else [()]
    keys = list(product(range(3), words))
    index = {word: i for i, word in enumerate(keys)}
    h = np.zeros((len(keys), len(keys)), complex)
    local = sum(
        a * b for a, b in zip((1, candidate.KAPPA), LAW.local_term_parts(), strict=True)
    )
    # Sign-even |d|=1 fiber; its coin is recorded separately on each seam.
    heads = np.zeros((11, len(range(2)) if written else 1))
    if written:
        for c in (0, 1):
            for d in (-1, 1):
                heads[1 + LAW.LABELS.index((c, d)), c] = 1 / np.sqrt(2)
    else:
        heads[0] = 1
    dress = []
    for col, (i, w) in enumerate(keys):
        dress.append(np.exp(-1j * sign * angle * w[0]) if written else 1)
        for e, (a, b) in enumerate(edges):
            initial = (1 + sign * (i == a), 1 + sign * (i == b))
            initial_index = LAW.SRC.index(initial)
            hi = w[e] if written else 0
            source = np.zeros(99)
            source[initial_index::9] = heads[:, hi]
            column = local @ source
            retained = 0.0
            for j in sorted({i, b if i == a else a if i == b else i}):
                final = (1 + sign * (j == a), 1 + sign * (j == b))
                component = heads.T @ column[LAW.SRC.index(final) :: 9]
                for hj, value in enumerate(component):
                    dest = w[:e] + (hj,) + w[e + 1 :] if written else w
                    movement = int(j == b) - int(i == b)
                    h[index[j, dest], col] += value * np.exp(
                        1j * angle * movement * (e == 0)
                    )
                    retained += abs(value) ** 2
            np.testing.assert_allclose(
                retained, np.vdot(column, column).real, atol=2e-12
            )
    return h, np.array(dress), keys


def native_flux_test(rng, shots):
    records = []
    for sign, written in product((-1, 1), (False, True)):
        h0, _, keys = triangle(sign, written, 0)
        h, gauge, _ = triangle(sign, written, 2 * np.pi / 20)
        diff = np.linalg.norm(h - (gauge[:, None] * h0 * gauge.conj()[None, :]), 2)
        if written:
            assert diff < 1e-11
        u0, u = expm(-1j * h0), expm(-1j * h)
        # One localized charge and a sharp coin background provide independent
        # positional records. Maximal shift over all preparations is also retained.
        probs0 = abs(u0) ** 2
        probs = abs(u) ** 2
        n = len(keys) // 3

        def read(p):
            return p.reshape(3, n, -1).sum(axis=1)

        if written:
            # Among src states at port 0, choose the first with transport on both paths.
            col = next(c for c in range(n) if min(read(probs0)[:, c]) > 1e-5)
        else:
            col = 0
        archives = []
        for p in (read(probs0)[:, col], read(probs)[:, col]):
            stable = np.maximum(0, np.round(p, 12))
            archives.append(rng.multinomial(shots, stable / stable.sum()).tolist())
        records.append(
            dict(
                sign=sign,
                written=written,
                source=keys[col],
                gauge_error=diff,
                spectrum_shift=float(
                    max(abs(np.linalg.eigvalsh(h) - np.linalg.eigvalsh(h0)))
                ),
                max_port_probability_shift=float(abs(read(probs) - read(probs0)).max()),
                port_probabilities=[
                    read(probs0)[:, col].tolist(),
                    read(probs)[:, col].tolist(),
                ],
                records=archives,
            )
        )
    return dict(full_seam=exact_columns(), triangle_reads=records)
