"""Exact native sign-basis and sharp readback checks."""

import numpy as np
from scipy.linalg import expm
from .exchange import LAW


def local_certificate():
    w = LAW._write_pattern().astype(np.int64)
    r, d2 = LAW._label_ops()
    swap, contrast = LAW._load_ops()
    mem = np.eye(10, dtype=np.int64)
    material = np.eye(9, dtype=np.int64)
    complement = 2 * np.eye(90, dtype=np.int64) - w @ w.T
    checks = 0
    for k in (0, 8):
        g = np.kron(mem - r + 2 * k * d2, material)
        h = complement @ g @ complement
        if np.any(complement @ g @ w):
            raise ArithmeticError("Blank and written blocks couple")
        if not np.array_equal(
            2 * w.T @ g @ w, 4 * (material - swap) + 8 * k * contrast
        ):
            raise ArithmeticError("Wrong blank block")
        for magnitude in (0, 1, 2):
            for sign in (1, -1) if magnitude else (1,):
                e = np.zeros((90, 18), dtype=np.int64)
                for c in range(2):
                    for index in range(9):
                        column = 9 * c + index
                        e[9 * LAW.LABELS.index((c, magnitude)) + index, column] = 1
                        if magnitude:
                            e[9 * LAW.LABELS.index((c, -magnitude)) + index, column] = (
                                sign
                            )
                mu = k * magnitude**2 + int(sign < 0)
                expected = 8 * mu * np.eye(18, dtype=np.int64)
                totals = sorted({a + b for a, b in LAW.SRC if abs(a - b) == magnitude})
                for total in totals:
                    pattern = np.zeros(18, dtype=np.int64)
                    for c in range(2):
                        # Output word depends only on total and the coin.
                        word = LAW.native_event(
                            *next((a, b) for a, b in LAW.SRC if a + b == total), c
                        )
                        pattern[9 * c + LAW.SRC.index(word)] = 1
                    expected -= 4 * mu * np.outer(pattern, pattern)
                if not np.array_equal(h @ e, e @ expected):
                    raise ArithmeticError("Written reduction differs")
                off = expected.copy()
                np.fill_diagonal(off, 0)
                if np.max(off) > 0:
                    raise ArithmeticError("Positive off-diagonal entry")
                checks += 18
    head = LAW.echo()
    exchange, cost = LAW.local_term_parts()
    swap_gate = expm(-1j * np.pi * (exchange + 8 * cost))
    np.testing.assert_allclose(swap_gate[:, :9], np.eye(99)[:, :9] @ swap, atol=4e-13)
    for load in range(3):
        sharp = np.eye(99)[:, LAW.SRC.index((load, 0))]
        np.testing.assert_allclose(cost @ sharp, load**2 * sharp, atol=4e-13)
    for col, (a, b) in enumerate(LAW.SRC):
        for delta in range(-2, 3):
            diag = np.r_[np.zeros(9), np.repeat([d == delta for c, d in LAW.LABELS], 9)]
            back = head.T @ (head[:, col] * diag)
            if np.linalg.norm(back - np.eye(99)[:, col] * int(a - b == delta)) > 1e-12:
                raise ArithmeticError("Sharp difference head failed")
    return dict(
        written_integer_columns_checked=checks,
        blank_dimensions=9,
        parity_blocks_per_contact=[
            "blank",
            "zero",
            "unit_even",
            "unit_odd",
            "double_even",
            "double_odd",
        ],
        native_sharp_read_columns=9,
        difference_outcomes=5,
    )


def reflection_sectors():
    """Reflection intertwines the two arms. Distinguish the source and target parities."""
    from .exchange import generator, source, target

    data = []
    for n in (1, 2):
        for kind in ("native", "exchange_minus", "flux_pi"):
            words, _, h = generator(n, 8, kind)
            rot = np.zeros_like(h)
            for col, w in enumerate(words):
                m = [-x % 8 for x in w]
                inversion = sum(a > b for j, a in enumerate(m) for b in m[j + 1 :])
                s = (-1) ** inversion if kind == "exchange_minus" else 1
                if kind == "flux_pi" and 0 in w:
                    s = -s
                rot[words.index(tuple(sorted(m))), col] = s
            if not np.array_equal(rot @ h, h @ rot):
                raise ArithmeticError("Guide reflection broken")
            parities = [
                int(rot[words.index(x), words.index(x)]) for x in (source(n), target(n))
            ]
            if not all(parities):
                raise ArithmeticError("Event has no reflection parity")
            data.append(
                {
                    "holes": n,
                    "law": kind,
                    "source_target_reflection_parities": parities,
                    "sharp_target_arm_sign": int(np.prod(parities)),
                }
            )
        hn = generator(n, 8, "native" if n == 1 else "flux_pi")[2]
        if not np.array_equal(hn, generator(n, 8, "exchange_minus")[2]):
            raise ArithmeticError("Fixed-particle ring statistics identity failed")
    return data
