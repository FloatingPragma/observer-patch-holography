"""Native energy capture from an independently flown load on the 80-cell federation."""

from fractions import Fraction
from math import factorial, isqrt

import numpy as np
from scipy import sparse

from model.experiments.interval_echo_einstein_station.mobile_collar_source import (
    measure as formation,
)
from matched_einstein_limit.source_implication.full_geometry.station_balls.reciprocal_screen.live_reference import (
    apparatus,
)

TICKS = (0, 32, 128)
NAMES = ("sharp", "forward", "reverse")
PRECISION = 2**54


def local_proof():
    law = formation.emission.LAW
    w = law._write_pattern()
    reversal, _ = law._label_ops()
    exchange, _ = law._load_ops()
    full_reversal = np.kron(reversal, np.eye(9, dtype=int))
    np.testing.assert_array_equal(w.T @ w, 2 * np.eye(9, dtype=int))
    np.testing.assert_array_equal(full_reversal @ w, w @ exchange)
    np.testing.assert_array_equal(
        w.T @ (w - full_reversal @ w), 2 * (np.eye(9) - exchange)
    )
    # A passive contact cannot write the parked sign, for *any* of the nine load words.
    basis = np.kron(
        np.diag([int(abs(d) == 1) for _, d in law.LABELS]), np.eye(9, dtype=int)
    )
    twice_odd = np.zeros((99, 99), dtype=int)
    twice_odd[9:, 9:] = basis @ (np.eye(90, dtype=int) - full_reversal)
    doubled = formation.emission.four_h(8)
    np.testing.assert_array_equal(doubled @ twice_odd, twice_odd @ doubled)
    np.testing.assert_array_equal(twice_odd @ twice_odd, 2 * twice_odd)
    # On capture the coin is maximally mixed and factors from the signed orientation.
    # Partial trace is checked on all four input matrix units, including coherences.
    cols = [law.SRC.index(word) for word in ((2, 1), (1, 2))]
    written = np.stack(
        [
            w.reshape(10, 9, 9)[
                [law.LABELS.index((c, d)) for c in range(2) for d in (1, -1)], :, j
            ].reshape(2, 2, 9)
            for j in cols
        ]
    )
    for a in range(2):
        for b in range(2):
            reduced = np.einsum("ijv,klv->ijkl", written[a], written[b])
            expected = np.zeros((2, 2, 2, 2), dtype=int)
            expected[:, a, :, b] = np.eye(2, dtype=int)
            np.testing.assert_array_equal(reduced, expected)
    # Independent discard-free energy assay on the blank load sector: Hadamard then position.
    _, hadamard, _ = apparatus.blank()
    np.testing.assert_allclose(
        hadamard.conj().T @ np.diag([0, 1]) @ hadamard,
        (np.eye(2) - np.fliplr(np.eye(2))) / 2,
        rtol=0,
        atol=2e-13,
    )


def apparatus_map():
    board = formation.federation(1)
    route = formation.region(board)
    # First native edge in the one-step ring of the left regular endpoint, outside capture.
    near = set(board.neighbours[route[0]]) - set(route)
    edge = min((a, b) for a, b in board.seams if a in near and b in near)
    assert edge == (5, 11) and board.degree[route[0]] == 5
    assert board.degree[edge[0]] == board.degree[edge[1]] == 5
    # Twice the one-excess blank-board Hamiltonian up to a common integer scalar.
    diagonal = 17 * board.degree + 48 * (6 - board.degree)
    center = (int(max(diagonal)) + int(min(diagonal))) // 2
    h = sparse.diags(diagonal - center, dtype=np.int64) - board.adjacency().astype(
        np.int64
    )
    h = h.tocsr()
    assert (h - h.T).count_nonzero() == 0
    return board, route, edge, h


def measure_vector(real, imag, route, mobile_norm2):
    a, b = route[:2]
    denom = int(mobile_norm2) * PRECISION**2

    def norm2(x, y):
        return int(x) ** 2 + int(y) ** 2

    capture = Fraction(norm2(real[a], imag[a]) + norm2(real[b], imag[b]), denom)
    energy = Fraction(norm2(real[a] - real[b], imag[a] - imag[b]), 2 * denom)
    return float(capture), float(energy)


def flights():
    """Fixed-point Taylor propagation with a whole-state norm error for each column."""
    board, route, edge, h = apparatus_map()
    real = np.zeros((board.n, 3), dtype=np.int64)
    imag = np.zeros_like(real)
    real[edge[0]] = PRECISION
    imag[edge[1]] = [0, PRECISION, -PRECISION]
    row = int(abs(h).sum(axis=1).max())
    order = 18
    assert row < 64
    per_step = Fraction(3 * order * (1 + isqrt(2 * board.n)), PRECISION)
    per_step += 2 * Fraction(row, 64) ** (order + 1) / factorial(order + 1)
    assert max(TICKS) * per_step < Fraction(1, 10**9)
    early_arrivals = {}
    out = {}
    for tick in range(max(TICKS) + 1):
        if tick == 1:
            for j, name in enumerate(NAMES):
                norm2 = 1 if j == 0 else 2
                early_arrivals[name] = float(
                    Fraction(
                        int(real[edge[1], j]) ** 2 + int(imag[edge[1], j]) ** 2,
                        norm2 * PRECISION**2,
                    )
                )
        if tick in TICKS:
            for j, (name, norm2) in enumerate(zip(NAMES, (1, 2, 2), strict=True)):
                cap, energy = measure_vector(real[:, j], imag[:, j], route, norm2)
                out[f"{name}_{tick // 32}"] = {"capture": cap, "exchange": energy}
        if tick == max(TICKS):
            break
        # Exact truncation bound uses a vector norm at most two before every step.
        squared = real.astype(object) ** 2 + imag.astype(object) ** 2
        assert max(squared.sum(axis=0)) < (2 * PRECISION) ** 2
        x, y = real.copy(), imag.copy()
        for k in range(1, order + 1):
            assert int(max(abs(x).max(), abs(y).max())) * row < 2**62
            x, y = (h @ y) // (64 * k), -(h @ x) // (64 * k)
            assert int(abs(real).max()) + int(abs(x).max()) < 2**62
            assert int(abs(imag).max()) + int(abs(y).max()) < 2**62
            real += x
            imag += y
    err = max(TICKS) * per_step
    # |<P>_hat-<P>_true| <= (2+delta)delta. Includes binary64 conversion.
    assert (2 + err) * err + Fraction(1, 10**15) < Fraction(1, 10**10)
    return out, {
        "max_row_twice_H": row,
        "vector_norm_error_upper": float(err),
        "capture_or_energy_probability_error_upper": 1e-10,
        "port_11_at_one_over_32": early_arrivals,
    }
