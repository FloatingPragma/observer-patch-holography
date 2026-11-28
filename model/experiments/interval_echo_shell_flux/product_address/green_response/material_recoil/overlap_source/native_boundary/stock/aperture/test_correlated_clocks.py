from fractions import Fraction

from mpmath import mp
import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    clock_covariance as clocks,
    clock_flight as flight,
    positive_resolvent as proof,
)


def test_clock_pair_inversion_for_distinct_correlated_lines():
    # A correlated law whose two off-diagonal high-frequency cells differ.
    # Combining the two positive lines into one would lose this test.
    matrix = clocks.source_receiver(2)
    inv = np.linalg.pinv(matrix)
    freq = np.array([[2, 1, 0], [0, 2, 1], [1, 3, 0]], dtype=float) / 10
    records = matrix @ freq @ matrix.T
    np.testing.assert_allclose(inv @ records @ inv.T, freq, atol=6e-13, rtol=0)


def test_exact_mixed_resolvent_against_dense_high_precision_inverse():
    twice = np.array(
        [
            [6, -1, 0, 0, 0],
            [-1, 6, -1, 0, 0],
            [0, -1, 6, -1, 0],
            [0, 0, -1, 6, -1],
            [0, 0, 0, -1, 6],
        ]
    )
    marks = np.array([[1, 1, 0, 0, 0], [0, 0, 0, 1, 1]])
    omega, gap = Fraction(3, 2), Fraction(1, 2)
    mixed, _ = proof.mixed_columns(twice, marks, omega, gap)
    with mp.workdps(120):

        def diagonal(a, b):
            h = mp.matrix(twice.tolist()) / 2 + mp.diag(
                (a * marks[0] + b * marks[1]).tolist()
            )
            return (h - mp.mpf(3) / 2 * mp.eye(5)) ** -1

        direct = mp.zeros(5)
        for a in (32, 33):
            for b in (32, 33):
                direct += (
                    diagonal(0, 0) - diagonal(a, 0) - diagonal(0, b) + diagonal(a, b)
                )
        for i in range(5):
            assert (
                0
                < mp.mpf(mixed[0][i]) / proof.S**3
                <= direct[i, i]
                <= mp.mpf(mixed[1][i]) / proof.S**3
            )


def test_four_port_clock_against_full_propagation_on_irregular_chart():
    h = np.diag([4.0, 5.0, 6.0, 4.0, 3.0])
    for i, j in ((0, 1), (0, 2), (2, 3), (3, 4)):
        h[i, j] = h[j, i] = -0.5
    marks = np.array([[1, 1, 0, 0, 0], [0, 0, 0, 1, 1]])
    psi = np.ones(5) / np.sqrt(5)
    maps = np.eye(5)
    result, _ = flight.boundary_response(h, marks, psi, maps, 4.375, -1.25)
    for i, a in enumerate(clocks.LINE):
        for j, b in enumerate(clocks.LINE):
            full = h - 1.25 * np.diag(a * marks[0] + b * marks[1])
            for n, age in enumerate(flight.MEAN_AGES):
                target = np.linalg.solve(
                    (1 - 1j * age * 4.375) * np.eye(5) + 1j * age * full, psi
                )
                np.testing.assert_allclose(result[i, j, n], target, atol=3e-14, rtol=0)
