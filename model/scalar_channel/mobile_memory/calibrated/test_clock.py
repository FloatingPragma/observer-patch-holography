import math

import numpy as np
from scipy import sparse

from .clock import audit


def _dense_grade_integral(source, beta, initial, wait):
    energies, modes = np.linalg.eigh(source)
    preparation = modes.conj().T @ initial
    transformed_grade = modes.conj().T @ np.diag(beta) @ modes
    frequencies = energies[:, None] - energies[None, :]
    phase = frequencies * wait / 2
    integrated_phases = wait * np.exp(1j * phase) * np.sinc(phase / np.pi)
    return float(np.vdot(preparation, (transformed_grade * integrated_phases) @ preparation).real)


def test_first_read_clock_has_exact_rest_phase_with_active_spectator():
    spectator = sparse.csr_matrix([[0.5, 0.25], [0.25, -0.5]])
    wait, gap = 2, 1 / 8
    result = audit(spectator, [1, 1], [1, 0], wait, gap=gap, groups=[7, 9])
    np.testing.assert_allclose(result["ramsey_overlap"], (math.cos(wait * gap), -math.sin(wait * gap)),
                               atol=2e-12, rtol=0)
    np.testing.assert_allclose(result["phase_clock_rate"], 1, atol=2e-12, rtol=0)
    np.testing.assert_allclose(result["exact_integrated_grade"], wait, atol=2e-12, rtol=0)
    assert result["terminal_groups"] == [7, 9]
    assert np.asarray(result["joint_records"]).shape == (2, 2, 2, 2)


def test_interacting_clock_matches_independent_dense_derivative():
    first = np.array([[1, 1], [1, -1]], dtype=float)
    second = np.array([[2, 1], [1, -2]], dtype=float)
    initial = np.kron(np.linalg.eigh(first)[1][:, -1], np.linalg.eigh(second)[1][:, -1])
    free = np.kron(first, np.eye(2)) + np.kron(np.eye(2), second)
    interaction = np.diag([0, 0, 0, 2])
    assert np.linalg.norm(free @ interaction - interaction @ free) > 1
    coupled = free + interaction
    beta = [1, 1, -1, -1]
    wait = 1.5
    result = audit(sparse.csr_matrix(coupled), beta, initial, wait, groups=[0, 1, 2, 3])
    control = audit(sparse.csr_matrix(free), beta, initial, wait)
    expected = _dense_grade_integral(coupled, beta, initial, wait)
    np.testing.assert_allclose(result["exact_integrated_grade"], expected, atol=2e-12, rtol=0)
    np.testing.assert_allclose(control["exact_grade_rate"], 1 / math.sqrt(2), atol=2e-12, rtol=0)
    assert abs(result["phase_clock_rate"] - control["phase_clock_rate"]) > 1e-3
    assert abs(result["phase_clock_rate"] - expected / wait) <= result["normalized_phase_error_bound"]
    assert abs(result["sine_clock_time"] - expected) <= result["sine_clock_time_error_bound"]
