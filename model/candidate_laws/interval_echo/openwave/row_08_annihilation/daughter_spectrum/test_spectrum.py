"""Independent checks of native gaps, finite clocks, and source/readout controls."""

from itertools import product
import json
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.interval_echo.openwave.row_08_annihilation.daughter_spectrum import measure as m


def independent_wing(boundary):
    words = list(product((0, 1), repeat=3))
    h = np.zeros((8, 8))
    for i, bits in enumerate(words):
        loads = [x for bit in bits for x in (2 * bit, 2 - 2 * bit)]
        # One attachment and the two contacts between, rather than inside, dimers.
        h[i, i] = 3 * 32.5 + 32.5 * (
            (boundary != loads[0]) + (loads[1] != loads[2]) + (loads[3] != loads[4]))
        for j in range(3):
            changed = tuple(bit ^ (k == j) for k, bit in enumerate(bits))
            h[words.index(changed), i] = -0.5
    return h


def test_native_gap_and_far_dimer_matrix_elements_from_loads():
    for boundary in (0, 2):
        h = independent_wing(boundary)
        np.testing.assert_allclose(h, m.ex.wing(8, 3, boundary), atol=1e-13)
        levels, vectors = np.linalg.eigh(h)
        np.testing.assert_allclose(levels[1:4] - levels[0],
                                   [31.8007556392, 32.5153809758, 33.2146344369], atol=1e-9)
        x = np.kron(np.eye(4), [[0, 1], [1, 0]])
        response = abs((vectors.T @ x @ vectors)[0, 1:4])**2
        np.testing.assert_allclose(response, [0.24976331, 0.49952663, 0.24976331], atol=1e-8)


def test_finite_time_bin_formula_by_explicit_clock_sum():
    h = independent_wing(0)
    n, step, target = 43, 0.031, 32.5153809758
    x = np.kron(np.eye(4), [[0, 1], [1, 0]])
    direct = np.zeros((8, 8), complex)
    for j in range(n):
        u = expm(-1j * h * j * step)
        direct += np.exp(1j * target * j * step) * u.conj().T @ x @ u / n
    effect, kraus = m.local_probe(h, target, n, step)
    np.testing.assert_allclose(kraus, direct, atol=2e-12)
    np.testing.assert_allclose(effect, direct.conj().T @ direct, atol=2e-12)


def bank():
    return json.loads((m.HERE / 'results.json').read_text())


def test_complete_two_phase_counts_and_resolved_late_response():
    data = bank()
    for arm in data['arms']:
        for panel in arm['readings']:
            for source in (0, 1):
                spectrum = panel['spectra'][source]
                assert {x['ring_phase'] for x in spectrum['spectrum_outcome_labels']} == {0, 1}
                assert len(spectrum['joint_spectrum']['counts']) == len(spectrum['spectrum_outcome_labels']) + 1
                assert sum(spectrum['joint_spectrum']['counts']) == m.SHOTS
                for record in panel['detector_records'][source]:
                    assert len(record['counts']) == 9
                    assert sum(record['counts']) == m.SHOTS
                    assert all(x >= 0 for x in record['probabilities'])
    late = data['arms'][0]['readings'][1]
    records = late['detector_records']
    mixed = next(r for r in late['source_excess'] if (r['left'], r['right']) == (1, 0))
    assert mixed['resolved_positive']
    assert mixed['exact_wall_minus_sham_double_click'] > 0.001
    assert 0.0005 < mixed['measured_wall_minus_sham_double_click'] < 0.0015
    assert records[1][3]['probabilities'][3] + records[1][3]['probabilities'][7] > 0.001


def test_control_backgrounds_and_degenerate_spectral_grouping():
    data = bank()
    for panel in data['arms'][1]['readings']:
        assert panel['cut_conditional_detector_factorization_residual'] < 1e-10
    # Zero contrast retains the unequal-load diagonal exchange terms; it is
    # not a sum of three independent flips. Audit labels without assuming degeneracy.
    for source in data['arms'][2]['readings'][0]['spectra']:
        labels = source['spectrum_outcome_labels']
        assert len({(x['ring_phase'], x['left_energy'], x['right_energy']) for x in labels}) == len(labels)
    for arm in data['arms'][1:]:
        assert all(not r['resolved_positive'] for p in arm['readings'] for r in p['source_excess'])


def test_invalid_probabilities_and_optimized_receipt_mismatch():
    rng = np.random.default_rng(33)
    with pytest.raises(ArithmeticError):
        m.sample([-1e-15, 1 + 1e-15], rng)
    with pytest.raises(ArithmeticError):
        m.sample([0.1, 0.2], rng)
    code = f"from {m.__package__} import measure; measure.compare({{'count':1}},{{'count':2}})"
    result = subprocess.run([sys.executable, '-O', '-B', '-c', code], capture_output=True, text=True)
    assert result.returncode != 0
    assert 'ArithmeticError' in result.stderr


def test_exactly_degenerate_energy_group_is_basis_independent():
    levels = np.array([0., 1., 1., 1., 2., 2., 2., 3.])
    local = [(levels, None), (levels, None)]
    rng = np.random.default_rng(111)
    state = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    state /= np.linalg.norm(state)
    weights, labels = m.group_spectrum(local, abs(state)**2, 0)
    rotation = np.eye(8)
    for indices in ([1, 2, 3], [4, 5, 6]):
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        rotation[np.ix_(indices, indices)] = q
    transformed = rotation @ state @ rotation.T
    repeated, same_labels = m.group_spectrum(local, abs(transformed)**2, 0)
    assert len(labels) == 16
    assert same_labels == labels
    np.testing.assert_allclose(weights, repeated, atol=2e-15)
