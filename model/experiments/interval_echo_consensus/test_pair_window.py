from fractions import Fraction

import numpy as np
import pytest
from scipy.linalg import expm

from .pair_window import TIMES, gram_ceiling, window_ceiling


@pytest.mark.parametrize("strength,limit", [(0, .934), (Fraction(1, 2), .939)])
def test_every_preparation_has_an_appointment_below_the_window(strength, limit):
    assert window_ceiling(strength) < limit < .94


def test_compression_keeps_complex_amplitudes_and_scalar_gauge():
    # A complex Hermitian generator on four sites with a one-site receiver.
    h = np.array([[0, 1j, 0, 0], [-1j, 0, 1, 0], [0, 1, 0, .5], [0, 0, .5, 1]])
    for shift in (0., .32):
        u = np.stack([expm(-1j*t*(h+shift*np.eye(4))) for t in TIMES])
        rows = u[:, :1, :].reshape(4, 4)
        actual = np.linalg.eigvalsh(rows.conj().T @ rows / 4)[-1]
        ceiling = gram_ceiling(u[1:, :1, :1], 1e-14, actual + .001)
        assert actual <= float(ceiling) < actual+.002
    with pytest.raises(ArithmeticError, match="does not enclose"):
        gram_ceiling(np.ones((3, 1, 1)), 0, .94)
    with pytest.raises(ValueError):
        gram_ceiling(np.full((3, 1, 1), np.nan), 0, .94)
