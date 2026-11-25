from fractions import Fraction
from itertools import combinations

import numpy as np

from model.executable_model import rotations
from model.experiments.seam_program.collective_rate.shell_response.measure import (
    PARITY_PORTS, finite_measurements, log_series, parities, region,
)


def test_every_three_port_read_and_complete_shell():
    ports = region()
    names, signs, norms = parities(list(range(12)))
    assert norms == [1, 30, 495]
    assert names == ["single", "rotated", "stationary"]
    words = np.arange(4096)
    for observed in combinations(range(12), 3):
        marginal = sum(((words >> p) & 1) << j for j, p in enumerate(observed))
        for sign in signs:
            assert np.all(np.bincount(marginal, weights=sign) == 0)
    # Check every rotated shell directly by summing the actual sharp probabilities.
    for rotation in rotations():
        shell = [rotation[i] for i in ports[2:]]
        marginal = sum(((words >> p) & 1) << j for j, p in enumerate(shell))
        assert np.all(np.bincount(marginal, weights=signs[1]) == 0)
    shell = sum(((words >> p) & 1) << j for j, p in enumerate(ports[2:]))
    assert np.all(np.bincount(shell, weights=signs[0]) == 0)
    for number in range(13):
        select = [i for i in words if int(i).bit_count() == number]
        assert np.ptp(signs[2, select]) == 0
    assert len({tuple(sorted(rotation[i] for i in PARITY_PORTS)) for rotation in rotations()}) == 30


def test_exact_parity_response_and_diagonal_control():
    leading = [Fraction(22897, 3840), Fraction(22897, 19200), Fraction(21309491, 950400)]
    controls = [0, 0, Fraction(8192, 495)]
    for exchange, fifth in ((True, leading), (False, controls)):
        results = log_series(exchange=exchange)
        assert list(results.values()) == [[0]*5 + [c] for c in fifth]
    slopes = finite_measurements()[1]["modular_slopes"]
    np.testing.assert_allclose(list(slopes.values()), [3.7207594e-6, 7.58593075e-7, 8.21019164e-5],
                               rtol=1e-7, atol=0)
