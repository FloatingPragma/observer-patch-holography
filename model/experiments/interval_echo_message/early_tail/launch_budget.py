"""Certify receiver click probabilities and repeat counts with integer arithmetic."""
from fractions import Fraction
import math

import numpy as np
from scipy.sparse import eye
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .. import channel


ORDER = 22


def amplitudes(centered_twice, ports, reciprocal_time):
    """Real and imaginary Taylor numerators over one common denominator."""
    scale = 2 * reciprocal_time
    denominator = scale**ORDER * math.factorial(ORDER)
    n = centered_twice.shape[0]
    term = np.zeros((n, 2), dtype=np.int64)
    term[ports, [0, 1]] = 1
    real = np.zeros_like(term, dtype=object)
    imag = np.zeros_like(real)
    # Integer powers are safe since every absolute row sum is at most six.
    if 6**ORDER > np.iinfo(np.int64).max:
        raise OverflowError("Taylor powers need wider integers")
    for order in range(ORDER + 1):
        coefficient = denominator // (scale**order * math.factorial(order))
        sign = (1, 1, -1, -1)[order % 4]
        if order % 2:
            imag -= sign * coefficient * term.astype(object)
        else:
            real += sign * coefficient * term.astype(object)
        if order < ORDER:
            term = centered_twice @ term
    error = Fraction(3, reciprocal_time)**(ORDER + 1) / math.factorial(ORDER + 1)
    return real, imag, denominator, 2 * error + error**2


def cell_intervals(data, cell):
    real, imag, denominator, error = data
    rows = slice(12 * cell, 12 * (cell + 1))
    norm_squares = (real[rows]**2 + imag[rows]**2).sum(axis=0)
    result = []
    for numerator in norm_squares:
        p = Fraction(numerator, denominator**2)
        result.append((max(Fraction(0), p - error), min(Fraction(1), p + error)))
    return result


def check(station):
    board = federation(station["level"])
    ports = station["ports"]
    encoded = channel.encode(board, station)
    np.testing.assert_allclose(np.abs(encoded[ports, [0, 1]]), 1, atol=1e-12, rtol=0)
    h = candidate.generator(board).tocsr()
    centered = 2 * (h - (6 * candidate.KAPPA + 3) * eye(board.n, format="csr"))
    integer = (board.laplacian() - 6 * eye(board.n, format="csr")).astype(np.int64)
    if (centered - integer).count_nonzero() or np.max(abs(integer).sum(axis=1)) > 6:
        raise ValueError("native one-load generator differs from the certified integer matrix")
    distance = shortest_path(board.adjacency(), directed=False, unweighted=True, indices=ports)
    minima = distance.reshape(2, board.cells, 12).min(axis=2).astype(int)
    native = amplitudes(integer, ports, reciprocal_time=100)
    half = amplitudes(integer, ports, reciprocal_time=200)
    target = Fraction(str(channel.SPEC["accuracy_target"]))
    for depth in range(1, 7):
        cells = np.flatnonzero((minima.min(axis=0) == depth) & (minima[0] != minima[1]))
        if not len(cells):
            print(f"{board.cells} cells, port distance {depth}: no eligible receiver")
            continue
        cell = int(cells[0])
        near = int(minima[:, cell].argmin())
        bands = cell_intervals(native, cell)
        low, high = bands[near]
        far_high = bands[1 - near][1]
        required = math.ceil((2 * target - 1) / max(high, far_high))
        shots = math.ceil(Fraction(1, 20) / low)
        lower_accuracy = (1 + shots * low / (1 + shots * low) - shots * far_high) / 2
        half_upper = (1 + shots * max(pair[1] for pair in cell_intervals(half, cell))) / 2
        if not half_upper < target < lower_accuracy or shots < required:
            raise ValueError("launch budget or half-exchange control missed the gate")
        print(f"{board.cells} cells, port distance {depth}, receiver {cell}, near bit {near}: "
              f"at least {required:,} independent launches, {shots:,} suffice; "
              f"native accuracy >= {float(lower_accuracy):.6f}, half exchange <= {float(half_upper):.6f}")


if __name__ == "__main__":
    for source in channel.SPEC["stations"]:
        check(source)
