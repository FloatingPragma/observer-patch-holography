"""Independent graph, Fourier, record and pulse checks."""

import json
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.interval_echo.openwave.row_31_orbitals.written_transitions import measure as m


def independent(board, seam):
    mat = np.zeros((board.n, board.n))
    for a, b in board.seams:
        if {a, b} == set(seam):
            continue
        mat[a, a] += 0.5
        mat[b, b] += 0.5
        mat[a, b] -= 0.5
        mat[b, a] -= 0.5
    for port in seam:
        mat[port, port] -= 8
    return mat


def test_native_matrix_and_two_sided_minmax_census_on_all_sizes():
    for level in range(3):
        board = m.federation(level)
        seam = m.EMITTER_SEAM[level]
        record = m.census(level)
        assert record["negative_level_count"] == 2
        assert record["negative_band"] == [-8., -5.5]
        np.testing.assert_allclose(record["endpoint_compression"], [[-6, 0], [0, -5.5]])
        if level < 2:
            direct = independent(board, seam)
            np.testing.assert_allclose(direct, m.generator(board, seam, 8), atol=1e-13)
            e = np.linalg.eigvalsh(direct)
            assert sum(e < 0) == 2 and e[2] > 0
            assert e[0] >= -8 and e[1] <= -5.5


def test_phase_kick_has_a_native_reference_operator():
    from model.candidate_laws.contrast_memory_spherical_v1 import law
    _, cost = law.local_term_parts()
    # Native diagonal-only reference contact: 8 D^2 for duration pi/8.
    pulse = expm(-1j * np.pi * cost)
    columns = np.eye(law.DIM)[:, [law.SRC.index((n, 0)) for n in (0, 1, 2)]]
    np.testing.assert_allclose(pulse @ columns, columns * [1, -1, 1], atol=3e-13)


def test_clock_filter_by_explicit_native_conjugated_pulses(monkeypatch):
    from model.executable_model import cell
    board = cell()
    h = independent(board, tuple(board.seams[0]))
    e, v = np.linalg.eigh(h)
    rng = np.random.default_rng(174)
    state = rng.normal(size=board.n) + 1j * rng.normal(size=board.n)
    state /= np.linalg.norm(state)
    n, step, target = 19, 0.031, 0.483
    monkeypatch.setattr(m, 'BINS', n)
    monkeypatch.setattr(m, 'STEP', step)
    pulse = np.eye(board.n)
    pulse[0, 0] = -1
    direct = np.zeros(board.n, complex)
    for j in range(n):
        u = expm(-1j * h * step * j)
        direct += np.exp(-1j * target * step * j) * u.conj().T @ pulse @ u @ state / n
    np.testing.assert_allclose(m.filtered(e, v, state, [0], target), direct, atol=2e-13)
    g = v[:, 0]
    np.testing.assert_allclose(m.stationary_filtered(e, v, g, [0], target),
                               m.filtered(e, v, g, [0], target), atol=2e-13)


def test_complete_counts_selection_hydrogen_refusal_and_held_profile():
    data = json.loads((m.HERE / 'results.json').read_text())
    for panel in data['panels']:
        assert len(panel['scan']) == 41
        for reading in panel['scan'] + [panel[k] for k in (
            'selected_native_line', 'common_phase_control', 'hydrogen_target',
            'planted_spectrum_detector_control', 'zero_contrast_same_input')]:
            assert len(reading['counts_by_port_then_clock_refusal']) == 12 * panel['cells'] + 1
            assert sum(reading['counts_by_port_then_clock_refusal']) == m.SHOTS
        assert panel['fitted_line']['held_probability_rms'] < 5e-5
        assert abs(panel['fitted_line']['frequency'] - 0.47993) < 2e-5
        assert panel['hydrogen_target']['counts_by_port_then_clock_refusal'][-1] == m.SHOTS
        assert panel['hydrogen_target']['acceptance'] < panel['hydrogen_target_acceptance_upper_bound'] < 2e-6
        assert panel['planted_spectrum_detector_control']['acceptance'] > 0.08
        assert panel['selected_native_line']['acceptance'] > 0.08
        assert panel['common_phase_control']['acceptance'] < 0.0004
        assert panel['zero_contrast_same_input']['acceptance'] < 4e-5
        assert panel['maximum_target_detuning'] < panel['alias_half_period']


def test_receipt_mutation_fails_with_optimization_and_invalid_norm_fails():
    with pytest.raises(ArithmeticError):
        m.record(np.array([1.001]), np.random.default_rng(1))
    code = f"from {m.__package__} import measure; measure.compare({{'count':2}},{{'count':1}})"
    result = subprocess.run([sys.executable, '-O', '-B', '-c', code], capture_output=True, text=True)
    assert result.returncode != 0 and 'ArithmeticError' in result.stderr


def test_empty_native_memory_suppresses_every_one_load_contact_column():
    from model.candidate_laws.contrast_memory_spherical_v1 import law
    x, c = law.local_term_parts()
    full = x + 8 * c
    zero = law.SRC.index((0, 0))
    empty_record = law.echo()[:, zero].reshape(-1, law.N_SRC)[:, zero]
    for occupation in ((0, 0), (1, 0), (0, 1)):
        written = np.kron(empty_record, np.eye(law.N_SRC)[:, law.SRC.index(occupation)])
        np.testing.assert_allclose(full @ written, 0, atol=2e-13)
